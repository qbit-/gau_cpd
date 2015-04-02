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
      integer    jwbuf, jxbufp, jxbufq, jscr, jebuf
      integer    lenwbuf, lenxbufp, lenxbufq, lenscr, lenebuf

C     other

      integer    nbasis, nbsq, ngrid
      integer    numk, numl, nump, numq
      integer    jnumk, jnuml, jnump, jnumq
      integer    jk, jl, jp, jq 
      integer    ikl

      logical    incorepos

      integer    irw2efull

C     functions
      
      integer    igetrw

      jend = 1 
      nbsq = nbasis*nbasis

      call ctrmemest(iout, iprint, nbasis, ngrid, incorepos, 
     &     numk, numl, nump, numq, lenwbuf, lenxbufp, lenxbufq,
     &     lenscr, v(jend), lenv) 
        
C---------------------------------------------------+      
C     Compute integrals                             |
C---------------------------------------------------+

      irw2efull = igetrw(nbsq*nbsq, 0)

      call frm2eints (iout, iprint, jsym2e, neq, neqshl, 
     &     neqsh2, incorepos, accdes, next, irw2efull, v(jend),
     & lenv-jend)

      jend = next

C     Memory allocation
      
      jwbuf  = jend
      jxbufp = jwbuf  + lenwbuf
      jxbufq = jxbufp + lenxbufp
      jscr   = jxbufq + lenxbufq
      jend   = jscr   + lenscr
      
      call tstcor (jend-1, lenv, 'frmemt: memory allocation')

      do 100 jk = 1, nbasis, numk
       jnumk = numk
       if (jk + numk - 1 .gt. nbasis) then
        jnumk = nbasis - jk*numk
       endif

       do 200 jl = 1, nbasis, numl
        jnuml = numl
        if (jl + numl - 1 .gt. nbasis) then
         jnuml = nbasis - jl*numl
        endif
            
        do 300 jp = 1, ngrid, nump
         jnump = nup
         if (jp + nup - 1 .gt. ngrid) then
          jnump = ngrid - jp*nump
         endif

         do 400 jq = 1, ngrid, numq
          jnumq = numq
          if (jq + numq - 1 .gt. ngrid) then
           jnumq = ngrid - jq*numq
          endif
            
C            call getklblock (iout, iprint, v(jend), lenv)
C            call getxblock ()

          do 500 icol = 1, jnumk*jnuml
                          
 500      continue

C            call transpose() 
C            call getxblock ()

          do 600 icol = 1, jnump

 600      continue
 400     continue
 300    continue
 200   continue
 100  continue

      end
