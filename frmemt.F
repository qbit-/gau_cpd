      subroutine frmemt(iout, iprint, nbasis, ngrid, jsym2e, neq,
     &     neqshl, neqsh2, accdes, irwfac, irwe, v, lenv)

      implicit none

C     This subroutine performs contractions of the two body integrals with the 
C     quadrature matrices. The contraction is:
C     
C     E(P,Q) = sum_{i,j} sum_{k,l} X(i,P) X(j,P) (ij | kl) X(k,Q) X(l,Q) 
C
C     The algorithm is built to fulfill two goals: maximum utilization of RAM
C     and the use of BLAS routines. It may not necessarily be optimal in parallel
C     case nor require less calls to integral generation/loading code
C
C     Code flow:  
C                    estimate memory requirements and sizes
C                    
C                    form 2e integrals in memory or on disk
C     
C                    for each block of kl indicies do in parallel 
C                        for each block of P 
C                            get (ij | kl)
C                            sum_{j} (ij | kl) R(j,P)
C                            L(iP | kl) = (iP | kl) * R(i,P)
C                            F(P, kl) = sum_{i} L(iP | kl)
C
C                            transpose F(P,kl)
C
C                            for each block of Q
C                                 sum_{l} F(P,kl) R(l,Q)
C                                 G(P,kQ) = F(P,kQ) * R(k,Q)
C                                 E(P,Q)  = sum_{k} F(Q,kQ)
C                                 
C                                 update E(P,Q) += E(P,Q)
C
C     
C     
C     Arguments:
C     
C     [in]  nbasis             - number of basis functions
C     [in]  ngrid              - number of grid points in the quadrature
C     ----passed to frm2eints---
C     [in]  jsym2e             - controls symmetry.           
C     [in]  neq, neqshl, neqs2 - symmetry related RWFs from Gaussian
C     [in]  accdes             - desired accuracy from PRISM
C     ----------------------
C     [in]  irwfac             - RWF with quadratures
C     [out] irwe               - RWF with E matrices 
C

C     Arguments 

      integer    iout, iprint, lenv
      real*8     v(*)
      
      integer    jsym2e, neq, neqshl, neqsh2
      real*8     accdes
      integer    irwfac, irwe

C     Pointers and sizes

      integer    jend, next
      integer    j2eint, jwbuf, jfbuf, jxbufp, jxbufq, jscr, jebuf
      integer    lenwbuf, lenfbuf, lenxbufp, lenxbufq, lenscr, lenebuf

C     other

      integer    nbasis, nbsq, ngrid
      integer    numkl, nump, numq
      integer    jnumkl, jnuml, jnumk, jnump, jnumq
      integer    jk, jl, jp, jq 
      integer    jkl, jklx 

      integer    ikl, ip, iq, ij
      integer    ioff, ldlen

      logical    incorepos, dbg

      integer    irw2efull
      integer    dims(4)
      
C     functions
      
      integer    igetrw
      real*8     asum


      jend = 1 
      nbsq = nbasis*nbasis

      call ctrmemest(iout, iprint, nbasis, ngrid, incorepos, 
     &     numkl, nump, numq, lenwbuf, lenfbuf, lenxbufp, lenxbufq,
     &     lenscr, lenebuf, v(jend), lenv) 

C---------------------------------------------------+      
C     Compute integrals                             |
C---------------------------------------------------+

      irw2efull = igetrw(nbsq*nbsq, 0)

      call frm2eints (iout, iprint, jsym2e, neq, neqshl, 
     &     neqsh2, incorepos, accdes, next, irw2efull, v(jend),
     & lenv-jend)

      if (incorepos) then
         j2eint = jend
         jend = next
      endif
      
C     Memory allocation
      
      jwbuf  = jend
      jfbuf  = jwbuf  + lenwbuf
      jxbufp = jfbuf  + lenfbuf
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
      
        if (jkl + numkl .gt. nbsq) then
           jnumkl = (nbsq - jkl)
        else
           jnumkl = numkl
        endif
        
        jnuml = jnumkl / nbasis
        if (jnuml .ne. 0) then
           jnumk  = nbasis
           jnumkl = jnuml*nbasis
        else
           jnuml  = 1
           jnumkl = jnumk
        endif
        
C     Load KL block from memory or file 

        ioff  = jkl*nbsq 
        ldlen = jnumkl*nbsq

        if (incorepos) then
C           call amove (ldlen, v(j2eint+ioff), v(jwbuf))
           jwbuf = j2eint + ioff
        else
           call fileio(2,-irw2efull,ldlen,v(jwbuf),ioff )
        endif

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

              do 215 ij = 0, jnump*nbasis - 1 
                 v(jscr + ij) = v(jscr + ij)*v(jxbufp + ij)
 215       continue

              do 220 ip = 0, jnump - 1 
                 v(jfbuf + ikl + ip*jnumkl) = asum (nbasis,
     & v(jscr + ip*nbasis), 1)
 220          continue

 210       continue
           
           do 300 jq = 0, ngrid - 1, numq
              if (jq + numq .gt. ngrid) then
                 jnumq = ngrid - jq
              else
                 jnumq = numq
              endif
              
C     Load right quadrature from file 

              ioff  = jq*nbasis
              ldlen = jnumq*nbasis 
              
              call fileio(2,-irwfac,ldlen,v(jxbufq),ioff)
C     
C     Do E(Q, P) = sum_{k,l} X(k, Q) X(l, Q) F(kl, P) and transpose it
C     
C     HERE WE NEED SUBMATRIX OF X in jxbufq . WILL FAIL IF L=1, K<NBASIS
              do 310 ip = 0, jnump - 1
               call dgemm ('n', 'n', jnumk, jnumq, jnuml, 1.0d0,  
     & v(jfbuf + ip*jnumkl), jnumk, v(jxbufq + jl), nbasis, 0.0d0, 
     & v(jscr), jnumk)
               do 320 iq = 0, jnumq - 1

                  do 325 ij = 0, jnumk - 1
                     v(jscr + ij + iq*jnumk) =
     & v(jscr + ij + iq*jnumk)*v(jxbufq + jk + iq*nbasis + ij)
 325              continue

                  v(jebuf + ip + iq*jnump) =
     &                 asum (jnumk, v(jscr + iq*jnumk), 1)
                  
 320           continue
 310        continue
            
            call updemt (iout, iprint, ngrid, jp, jnump, jq, jnumq, 
     &             v(jebuf), irwe, v(jend), lenv)
 300       continue
 200    continue

        jk = mod(jk + jnumk, nbasis)
        jl = mod(jl + jnumkl/nbasis, nbasis)
        
        jkl = jkl + jnumkl
        if (jkl .lt. nbsq - 1) goto 100
C     end of the contraction loop

        dbg = .true.
        if (dbg .and. incorepos) then
           dims(1) = nbasis
           dims(2) = nbasis
           dims(3) = nbasis
           dims(4) = nbasis
           call outoctfmt(iout, 'ints.mat', 'ints', 4, dims, v(j2eint))
        endif
        
        end
      
