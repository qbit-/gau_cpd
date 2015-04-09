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
      integer    jnuml, jnumk, jnump, jnumq, jnumkl
      integer    jk, jl, jp, jq 
      integer    jkl, jklx 

      integer    ikl, ip, iq, ij
      integer    ioff, ldlen

      logical    incorepos

      integer    irw2efull

C     functions
      
      integer    igetrw
      real*8     asum


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

C     Allocate RWF for the E matrix
      
      call conddf (irwe, ngrid*ngrid)

C     Top of the contraction loop
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
        jnumkl = numk + numl*nbasis

C     Load KL block from memory or file 

        ioff  = (jk-1)*nbsq + (jl-1)*nbasis*nbsq
        ldlen = numk*nbsq + numl*nbsq*nbasis

        if (incorepos) then
           call amove (ldlen, v(j2eint+ioff), v(jwbuf))
        else
           call fileio(2,-irw2efull,ldlen,v(jwbuf),ioff )
        endif

C     Loop over left quadrature

        do 200 jp = 1, ngrid, nump
           if (jp + nump .gt. ngrid) then
              jnump = ngrid - jp
           else
              jnump = nump
           endif
           
C     Load letf quadrature block from file

           ioff  = (jp-1)*nbasis
           ldlen = jnump*nbasis 
           
           call fileio(2,-irwfac,ldlen,v(jxbufp),ioff)
C
C     Do F(kl, P) = sum_{i,j} X(i,P) X(j,P) (i,j| kl)
C 
           do 210 ikl = 1, jnumkl
              call dgemm ('n', 'n', nbasis, jnump, nbasis, 1.0d0,
     & v(jwbuf + (ikl-1)*nbsq), nbasis, v(jxbufp), nbasis, 0.0d0,
     & v(jscr), nbasis)

              do 215 ij = 1, nbsq
                 v(jscr + ij - 1) = v(jscr + ij - 1)*v(jxbufp + ij - 1)
 215       continue

              do 220 ip = 1, jnump 
                 v(jfbuf + (ikl-1) + (ip-1)*nbasis) = asum (nbasis,
     & v(jscr + (ip-1)*nbasis), 1)
 220          continue

 210       continue
           
           do 300 jq = 1, ngrid, numq
              if (jq + numq .gt. ngrid) then
                 jnumq = ngrid - jq
              else
                 jnumq = numq
              endif
              
C     Load right quadrature from file 

              ioff  = (jq-1)*nbasis
              ldlen = jnumq*nbasis 
              
              call fileio(2,-irwfac,ldlen,v(jxbufq),ioff)
C     
C     Do E(Q, P) = sum_{k,l} X(k, Q) X(l, Q) F(kl, P) and transpose it
C     
C     HERE WE NEED SUBMATRIX OF X in jxbufq . WILL FAIL IF L=1, K<NBASIS
              do 310 ip = 1, jnump
               call dgemm ('n', 'n', jnumk, jnumq, jnuml, 1.0d0,  
     & v(jfbuf + (ip-1)*jnumkl), jnumkl, v(jxbufq + jk-1),nbasis, 0.0d0, 
     & v(jscr), jnumk)
               do 320 iq = 1, jnumq

                  do 325 ij = 1, jnumk
                     v(jscr + (iq-1)*jnumk + ij-1) =
     & v(jscr + (iq-1)*jnumk + ij-1)*v(jxbufq+jk-1+(iq-1)*nbasis+ij-1)
 325              continue

                  v(jebuf + (ip-1) + (iq-1)*jnump) = asum (jnumk, 
     & v(jscr + (iq-1)*jnumk), 1)

 320           continue
 310          continue

              call updemt (iout, iprint, ngrid, jp, jnump, jq, jnumq, 
     &             v(jebuf), irwe, v(jend), lenv)
 300       continue
 200  continue

        jkl = jklx
        if (jkl .lt. nbsq) goto 100
C     end of the contraction loop

      end
