      subroutine frmemtri(iout, iprint, nbasis, ngrid, jsym2e, neq,
     &     neqshl, neqsh2, ir12, tolinv, accdes, irwfac, irwe, v, lenv)

      implicit none
      
C     This subroutine performs contractions of the RI decomposed two body integrals with the 
C     quadrature matrices. The contraction is:
C     
C     E(P,Q) = sum_{i,j} sum_{k,l} X(i,P) X(j,P) (ij | L) S(L,M)^-1 (M, kl) X(k,Q) X(l,Q) 
C
C     The algorithm is built to fulfill two goals: maximum utilization of RAM
C     and the use of BLAS routines. This algorithm may have to be merged into full
C     integral contraction code.
C
C     Code flow:  
C                    estimate memory requirements and sizes
C                    
C                    form RI 2e integrals and overlap S in memory or on disk
C     
C                    for each block of L DBF indicies do in parallel 
C                            get (kl | L)
C                            sum_{j} (ij | L) R(j,P)
C                            O(iP | L) = (iP | L) * R(i,P)
C                            F(P, L) = sum_{i} O(iP | L)
C
C                            E(P,Q)  = sum_{L,M} F(P,L) S(L,M)^-1 F(M,Q)
C                                 
C                            update E(P,Q) += E(P,Q)
C
C     
C     
C     Arguments:
C     
C     [in]  nbasis             - number of basis functions
C     [in]  ngrid              - number of grid points in the quadrature
C     ----passed to frm2eri---
C     [in]  jsym2e             - controls symmetry.           
C     [in]  neq, neqshl, neqs2 - symmetry related RWFs from Gaussian
C     [in]  ir12               - metric to use: 1 - 1/r, 2 - overlap
C     [in]  tolinv             - tolerance for the inverse of the overlap
C     [in]  accdes             - desired accuracy from PRISM
C     ----------------------
C     [in]  irwfac             - RWF with quadratures
C     [out] irwe               - RWF with E matrices 
C

#include "commomob2.inc"
      
C     Arguments 

      integer    iout, iprint, lenv
      real*8     v(*)
      
      integer    jsym2e, neq, neqshl, neqsh2
      integer    ir12
      real*8     tolinv, accdes
      integer    irwfac, irwe

C     Pointers and sizes

      integer    jend, next
      integer    jov, j2eint,jwbuf,jovbuf, jfbuf, jcbuf, jxbufp, jxbufq
      integer    jscr, jebuf
      integer    lenwbuf, lenfbuf, lenxbufp, lenxbufq, lenscr, lenebuf

C     other

      integer    nbasis, nbsq, ngrid, ndbf
      integer    numkl, nump, numq
      integer    jnumkl, jnuml, jnumk, jnump, jnumq
      integer    jk, jl, jp, jq 
      integer    jkl, jklx 

      integer    ikl, ip, iq, ik, il, ij
      integer    ioff, ldlen

      logical    keep, dbg

      integer    irwria, irwrib, irwric, irwrid
      integer    dims(4)
      
C     functions
      
      integer    igetrw
      real*8     ddot, dsdot

C     Load /B2/ with the density basis.

      call rdwrb(-10,0)
      
C     Figure out NDBF and NDBF6D.

      call getmbn(I5DB2,I7FB2,NShelB,ShelTB,ShelCB,ndbf)
      call getmbn(1,1,NShelB,ShelTB,ShelCB,ndbf6d)

      jend   = 1 
      nbsq   = nbasis*nbasis
      ndb6sq = ndbf6d*ndbf6d
      
      call ctrmemest2(iout, iprint, nbasis, ngrid, ndbf,
     &     keepov, keepints, numkl, nump, numq, lenwbuf, lenfbuf,
     &     lenxbufp, lenxbufq, lenscr, lenebuf, v(jend), lenv - jend) 

C---------------------------------------------------+      
C     Compute integrals                             |
C---------------------------------------------------+

      call frm2eri (iout, iprint, jsym2e, neq, neqshl, 
     &     neqsh2, ir12, keepov, keepints, tolinv, accdes, jovoff,
     &     next, irwria, irwrib, v(jend), lenv - jend)

C     Allocate files for half projected integrals and half projected integrals times overlap 
      irwric = igetrw(ngrid*ndbf6d)
      irwrid = igetrw(ngrid*ndbf6d)
      
      if (keepov) then
         jov    = jend
         jend   = jend + jovoff
      endif

      if (keepints) then
         j2eint = jend
         jend   = jend + next
      endif
      
C     Memory allocation
      
      jwbuf  = jend
      jovbuf = jwbuf  + lenwbuf
      jfbuf  = jovbuf + lenfbuf
      jcbuf  = jfbuf  + lenfbuf
      jxbufp = jcbuf  + lenfbuf
      jxbufq = jxbufp + lenxbufp
      jscr   = jxbufq + lenxbufq
      jebuf  = jscr   + lenscr
      jend   = jebuf  + lenebuf

      write (iout, *) 'blah from frmemt'
      call tstcor (jend-1, lenv, 'frmemt: memory allocation')

C     Allocate RWF for the E matrix
      
      call conddf (irwe, ngrid*ngrid)

C     Top of the contraction loop
      jkl = 0
      jk  = 0
      jl  = 0
 100  continue

C     Ugly stepsize computation. We always round to whole counts 
C     of L index if possible
      
        if (jkl + numkl .gt. ndbf6d) then
           jnumkl = (ndbf6d - jkl)
        else
           jnumkl = numkl
        endif
        
C     Load KL block from memory or file 

        ioff  = jkl*nbsq 
        ldlen = jnumkl*nbsq

        if (keepints) then
C           call amove (ldlen, v(j2eint+ioff), v(jwbuf))
           jwbuf = j2eint + ioff
        else
           call fileio(2,-irwrib,ldlen,v(jwbuf),ioff )
        endif

C     Load overlap block
        
C     Loop over left quadrature

        do 200 jp = 0, ngrid - 1, nump
           if (jp + nump .gt. ngrid) then
              jnump = ngrid - jp
           else
              jnump = nump
           endif
           
C     Load letf quadrature block from file

           ioff  = jp*nbasis
           ldlen = jnump*nbasis 
           
           call fileio(2,-irwfac,ldlen,v(jxbufp),ioff)
C
C     Do F(kl, P) = sum_{i,j} X(i,P) X(j,P) (i,j| kl)
C 
           do 210 ikl = 0, jnumkl - 1
              call dgemm ('n', 'n', nbasis, jnump, nbasis, 1.0d0,
     & v(jwbuf + ikl*nbsq), nbasis, v(jxbufp), nbasis, 0.0d0,
     & v(jscr), nbasis)

              do 220 ip = 0, jnump - 1 
                 v(jfbufp + ikl + ip*jnumkl) = ddot (nbasis,
     & v(jscr + ip*nbasis), 1, v(jxbufp + ip*nbasis), 1)
 220          continue
 210       continue
   
           call updemt(C)
           call updemt(D)
           
 200    continue
           
        jkl = jkl + jnumkl
        if (jkl .lt. ndbf6d) goto 100
C     end of the contraction loop

C     Block matrix multiplication
        do jp = 0, ngrid - 1
           do jq = 0, ngrid - 1 
C     use full ndbf6d blocks
              call updemt(E)
           enddo
        enddo
        
        dbg = .true.
        if (dbg .and. keepints) then
           dims(1) = nbasis
           dims(2) = nbasis
           dims(3) = ndbf6d
           call outoctfmt(iout, 'ints.mat', 'ints', 3, dims, v(j2eint))
        endif
        
        end
      
      
