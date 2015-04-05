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
      integer    numk, numl, nump, numq
      integer    jnuml, jnumk, jnump, jnumq
      integer    jk, jl, jp, jq 
      integer    jkl, jklx 

      logical    incorepos

      integer    irw2efull

C     functions
      
      integer    igetrw

      jend = 1 
      nbsq = nbasis*nbasis

      call ctrmemest(iout, iprint, nbasis, ngrid, incorepos, 
     &     numk, numl, nump, numq, lenwbuf, lenfbuf, lenxbufp, lenxbufq,
     &     lenscr, lenebuf, v(jend), lenv) 
        
C---------------------------------------------------+      
C     Compute integrals                             |
C---------------------------------------------------+

      irw2efull = igetrw(nbsq*nbsq, 0)

      call frm2eints (iout, iprint, jsym2e, neq, neqshl, 
     &     neqsh2, incorepos, accdes, next, irw2efull, v(jend),
     & lenv-jend)

      if (incorepos) j2eint = jend
      jend = next

C     Memory allocation
      
      jwbuf  = jend
      jfbuf  = jwbuf  + lenwbuf
      jxbufp = jfbuf  + lenfbuf
      jxbufq = jxbufp + lenxbufp
      jscr   = jxbufq + lenxbufq
      jebuf  = jscr   + lenscr
      jend   = jebuf  + lenebuf
      
      call tstcor (jend-1, lenv, 'frmemt: memory allocation')

      jkl = 1
 100  continue

C     Ugly stepsize computation. We try to load by whole counts 
C     of L index, e.g. nbasis*nbsq

        jklx = jkl + (numk + numl*nbasis)    
        if (jklx / nbasis .gt. jkl / nbasis) then
           if (jklx .gt. nbsq) then
              jnuml = (nbsq - jkl) / nbasis
              jnumk = mod((nbsq - jkl), nbasis)
           else
              jnuml = (jklx - jkl) / nbasis
              jnumk = mod((jklx - jkl), nbasis)
           endif
        endif
        jl = jkl / nbasis + 1
        jk = mod(jkl, nbasis) + 1
        
        call gtklblock (iout, iprint, jk, jl, jnumk, jnuml, 
     &       incorepos, v(jwbuf), v(j2eint), irw2efull, v(jend), lenv)
        
        
        jkl = jklx
        if (jkl .lt. nbsq) goto 100

      end
