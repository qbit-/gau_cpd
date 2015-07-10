      subroutine frmemtri(iout, iprint, nbasis, ngrid, jsym2e, neq,
     &     neqshl, neqsh2, ir12, tolinv, accdes, ndbf,
     &     irwfac, irwria, irwrib,
     &     irwe, v, lenv)

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

#include "commonb2.inc"
C     Needed for commonb2 to work with implicit none

      integer    nshelb, janb, maxtyb 
C     Arguments 

      integer    iout, iprint, lenv
      real*8     v(*)
      
      integer    jsym2e, neq, neqshl, neqsh2
      integer    ir12
      real*8     tolinv, accdes
      integer    irwfac, irwe

C     Pointers and sizes

      integer    jend, next, jovoff
      integer    jov, j2eint,jwbuf, jfbuf, jxbufp, jxbufq
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

      logical    keepints, keepov, dbg

      integer    irwria, irwrib, irwric, irwrid
      integer    dims(4)
      
C     functions
      
      integer    igetrw
      real*8     ddot, dsdot

C     Load /B2/ with the density basis.

      call rdwrb(-10,0)
      
C     Figure out NDBF

      call getmbn(I5DB2,I7FB2,NShelB,ShelTB,ShelCB,ndbf)

      jend   = 1 
      nbsq   = nbasis*nbasis
          
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
      irwric = igetrw (ngrid*ndbf, 0)
      irwrid = igetrw (ngrid*ndbf, 0)
      
      if (keepov) then
         jov    = jend
         jend   = jend + jovoff
      else
         jov    = jend
         jend   = jov  + ndbf*ndbf
         call   tstcor (jend-1, lenv, 'frmemtri: overlap')
         call   fileio (2, -irwria, ndbf*ndbf, v(jov), 0)
      endif

      if (keepints) then
         j2eint = jend
         jend   = jend + next
      endif
      
C     Memory allocation
      
      jwbuf  = jend
      jfbuf  = jwbuf  + lenwbuf
      jxbufp = jfbuf  + lenfbuf
      jxbufq = jxbufp + lenxbufp
      jscr   = jxbufq + lenxbufq
      jebuf  = jscr   + lenscr
      jend   = jebuf  + lenebuf

      call tstcor (jend-1, lenv, 'frmemtri: memory allocation')

C     Allocate RWF for the E matrix
      
      call conddf (irwe, ngrid*ngrid)

C     Top of the contraction loop
      jkl = 0
      jk  = 0
      jl  = 0
 100  continue

C     Ugly stepsize computation. We always round to whole counts 
C     of L index if possible
      
        if (jkl + numkl .gt. ndbf) then
           jnumkl = (ndbf - jkl)
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
                 v(jfbuf + ikl + ip*jnumkl) = ddot (nbasis,
     & v(jscr + ip*nbasis), 1, v(jxbufp + ip*nbasis), 1)
 220          continue
 210       continue
   
           call updmatf (iout, iprint, ndbf, jkl, jnumkl, ngrid, 
     &      jp, jnump, v(jfbuf), irwric, 0, v(jscr), lenscr)
           
 200    continue
           
        jkl = jkl + jnumkl
        if (jkl .lt. ndbf) goto 100
C     end of the contraction loop

C     Form A(L, Q) = S(L, M)^(-1) * F(M, Q)
        
        do 300 jp = 0, ngrid - 1, nump
           if (jp + nump .gt. ngrid) then
              jnump = ngrid - jp
           else
              jnump = nump
           endif
           
           ioff  = jp*ndbf
           ldlen = jnump*ndbf
           
           call fileio (2, -irwric, ldlen, v(jxbufp), ioff)  
           call dgemm ('n', 'n', ndbf, jnump, ndbf, 1.0d0,
     & v(jov), ndbf, v(jfbuf), ndbf, 0.0d0, v(jscr), ndbf)
           call fileio (1, -irwrid, ldlen, v(jscr), ioff)

 300    continue
        
C     Block matrix multiplication to form E(P, Q) =  F(P, L) * A(L, Q) 
        do 400 jp = 0, ngrid - 1, nump
           if (jp + nump .gt. ngrid) then
              jnump = ngrid - jp
           else
              jnump = nump
           endif
           
           ioff  = jp*ndbf
           ldlen = jnump*ndbf
           
           call fileio (2, -irwric, ldlen, v(jxbufp), ioff)

           do 500 jq = 0, ngrid - 1, numq 
              if (jq + numq .gt. ngrid) then
                 jnumq = ngrid - jq
              else
                 jnumq = numq
              endif
              
              ioff  = jq*ndbf
              ldlen = jnumq*ndbf

              call fileio (2, -irwrid, ldlen, v(jxbufq), ioff)

              call dgemm ('t', 'n', jnump, jnumq, ndbf, 1.0d0,
     &    v(jxbufp), ndbf, v(jxbufq), ndbf, 0.0d0, v(jebuf), jnump)  

              call updmatf (iout, iprint, ngrid, jp, jnump, ngrid, jq,
     &    jnumq, v(jebuf), irwe, 1, v(jscr), lenscr)
 500       continue
 400    continue
        
        end
      
      
