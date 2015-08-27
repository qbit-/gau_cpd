      subroutine frmemt_nd(iout, iprint, nbasis, ngrid, jsym2e, neq,
     &     neqshl, neqsh2, accdes, irwfac, irw2efull, irwe, v, lenv)

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
C     estimate memory requirements and sizes
C     
C     for each block of P 
C       for each block of kl indicies do in parallel 
C       get (ij | kl)
C          L(iP | kl) = sum_{j} (ij | kl) R(j,P)
C          F(P, kl) = sum_{i} R(i,P) * L(iP | kl)
C     
C          for each block of Q
C          G(l,k,Q) = R(l,Q)*R(k,Q) (Direct product)
C          E(P,Q)  = sum_{k,l} F(P,kl) * G(lk,Q)
C     
C          update E(P,Q) += E(P,Q)
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
C     [out] irw2efull          - RWF with 2e integrals
C     [out] irwe               - RWF with E matrices 
C     

C     Arguments 

      integer    iout, iprint, lenv
      double precision     v(*)
      
      integer    jsym2e, neq, neqshl, neqsh2
      double precision     accdes
      integer    irwfac, irwe, irw2efull

C     Pointers and sizes

      integer    jend, j2eint, jx, next
      integer    jwbuf, jfbuf, jxbufp, jxbufq, jscr, jebuf
      integer    lenwbuf, lenfbuf, lenxbufp, lenxbufq, lenscr, lenebuf

C     other

      integer    nbasis, nbsq, ngrid
      integer    numkl, nump, numq
      integer    jnumkl, jnuml, jnumk, jnump, jnumq
      integer    jk, jl, jp, jq 
      integer    jkl
      integer    iupd
      
      integer    ikl, ip, iq
      integer    ioff, ldlen

      logical    keepints, keepx

C     integer    irw2efull

      
C     functions

      double precision     ddot
      

      jend = 1 
      nbsq = nbasis*nbasis

      call ctrmemest_nd(iout, iprint, nbasis, ngrid, keepints, 
     &     keepx, numkl, nump, numq, lenwbuf, lenfbuf, lenxbufp,
     &     lenxbufq, lenscr, lenebuf, v(jend), lenv - jend) 

C---------------------------------------------------+      
C     Compute integrals                             |
C---------------------------------------------------+

C     Load required quantities if we keep them in memory

      if (keepints) then
         j2eint = jend
         call frm2eints (iout, iprint, jsym2e, neq, neqshl, 
     &     neqsh2, accdes, keepints,.false., next, irw2efull, v(j2eint),
     & lenv-jend)
         jend = jend + next
      endif

      if (keepx) then
         jx   = jend
         call fileio(2,-irwfac,ngrid*nbasis,v(jx),0)
         jend = jx + ngrid*nbasis
      endif

C     Allocate RWF for the E matrix
      
      call conddf (irwe, ngrid*ngrid)

C     Fill RWF with zeros if we are not going to write
C     the full size matrix in the first pass

      if ((nump .ne. ngrid) .or. (numq .ne. ngrid)) then
         call aclearf (ngrid*ngrid, irwe, v(jend), lenv - jend)
      endif
      
C     Memory allocation

      jwbuf  = jend
      jfbuf  = jwbuf  + lenwbuf
      jxbufp = jfbuf  + lenfbuf
      jxbufq = jxbufp + lenxbufp
      jscr   = jxbufq + lenxbufq
      jebuf  = jscr   + lenscr
      jend   = jebuf  + lenebuf

      call tstcor (jend-1, lenv, 'frmemt: memory allocation')

C     Top of the contraction loop
      iupd = 0
      
C     Loop over left quadrature

      do 100 jp = 0, ngrid - 1, nump
         if (jp + nump .gt. ngrid) then
            jnump = ngrid - jp
         else
            jnump = nump
         endif
         
C     Load letf quadrature block from file

         ioff  = jp*nbasis
         ldlen = jnump*nbasis 

         if (keepx) then
            jxbufp = jx + ioff
         else
            call fileio(2,-irwfac,ldlen,v(jxbufp),ioff)
         endif
         
         jkl = 0
         jk  = 0
         jl  = 0
 200     continue

C     Ugly stepsize computation. We always round to whole counts 
C     of L index if possible
         
         if (jkl + numkl .gt. nbsq) then
            jnumkl = (nbsq - jkl)
         else
            jnumkl = numkl
         endif
         
         jnuml = jnumkl / nbasis
         if (jnuml .gt. 0) then
            jnumk  = nbasis
            jnumkl = jnuml*nbasis
         else
            jnuml  = 1
            jnumk  = min(jnumkl, nbasis - jk)
            jnumkl = jnumk
         endif
         
C     Load KL block from memory or file 

         ioff  = jkl*nbsq 
         ldlen = jnumkl*nbsq

         if (keepints) then
C     call amove (ldlen, v(j2eint+ioff), v(jwbuf))
            jwbuf = j2eint + ioff
         else
            call fileio(2,-irw2efull,ldlen,v(jwbuf),ioff )
         endif

C     
C     Do F(kl, P) = sum_{i,j} X(i,P) X(j,P) (i,j| kl)
C     
         do 210 ikl = 0, jnumkl - 1
            call dgemm ('n', 'n', nbasis, jnump, nbasis, 1.0d0,
     &           v(jwbuf + ikl*nbsq), nbasis, v(jxbufp), nbasis, 0.0d0,
     &           v(jscr), nbasis)

            do 220 ip = 0, jnump - 1 
               v(jfbuf + ip + ikl*jnump) = ddot (nbasis,
     &              v(jxbufp + ip*nbasis), 1, v(jscr + ip*nbasis), 1)
 220        continue
 210     continue
         
         do 300 jq = 0, ngrid - 1, numq
            if (jq + numq .gt. ngrid) then
               jnumq = ngrid - jq
            else
               jnumq = numq
            endif
            
C     Load right quadrature from file 

            ioff  = jq*nbasis
            ldlen = jnumq*nbasis 

            if (keepx) then
               jxbufq = jx + ioff
            else
               call fileio(2,-irwfac,ldlen,v(jxbufq),ioff)
            endif
C     
C     Do E(P, Q) = sum_{k,l} F(P, kl) X(k, Q) X(l, Q)
C
            call aclear(jnumq*jnumkl, v(jscr))
            do 310 iq = 0, jnumq - 1
               call dger (jnuml, jnumk, 1.0d0,
     &          v(jxbufq+jl+iq*nbasis), 1, v(jxbufq+jk+iq*nbasis), 1,  
     &              v(jscr+iq*jnumkl), jnuml)
 310        continue
            
            call dgemm ('n', 'n', jnump, jnumq, jnumkl, 1.0d0,
     &           v(jfbuf), jnump, v(jscr), jnumkl, 0.0d0,
     &           v(jebuf), jnump)

            call updmatf (iout, iprint, ngrid, jp, jnump, ngrid, jq,
     &           jnumq, v(jebuf), irwe, iupd, v(jend), lenv - jend)
            if (iupd .eq. 0) iupd = 1
 300     continue

         jkl = jkl + jnumkl
         jk = mod(jk + jnumk, nbasis)
         jl = jkl / nbasis
         
         if (jkl .lt. nbsq) goto 200
 100     continue

C     end of the contraction loop

         end
      
