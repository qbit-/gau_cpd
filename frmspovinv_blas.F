      subroutine frmspovinv (iout, iprint, nbas6d, ngrid, tol,
     &     irwfac, irwspo, v, lenv)

      implicit none

C     Computes the Moore-Penrose pseudoinverse of the matrix
C     S(P,Q) = ( X(mu,P)^t X(mu,Q) ) o ( X(nu,P)^t X(nu,Q) )
C     X is stored in the irwfac, s is put into irwspo
C
C     Options:
C     [in] nbas6d - cartesian basis size (dimension of mu)
C     [in] ngrid  - number of grid points (dimension of P)
C     [in] tol    - numeric tolerance for the pseudoinverse
C     [in] irwfac - RWF with quadrature
C     [out] irwspo - RWF with S^-1 

      integer    iout, iprint, nbas6d, ngrid
      integer    irwfac, irwspo, lenv
      real*8     tol, v
      dimension  v(*)

      real*8     one, zero
      parameter  ( one = 1.0d0, zero = 0.0d0)
      
C     other variables
      integer    jend, jx, js, jsv, ju, jvt, jscr, ii
      integer    jscopy, junit
      integer    lenx, lens, lensv, lenu, lenvt, lenscr
      integer    info
      real*8     l(ngrid,ngrid)
      real*8     r(ngrid,ngrid)
C     define sizes
      lenx = nbas6d * ngrid
      lens = ngrid*ngrid
      lenu = lens
      lenvt = lens
      lensv = ngrid
      lenscr = 8*ngrid
      
C     allocate memory
      jend = 1
      jx   = jend
      js   = jx + lenx
      jsv  = js + lens
      ju   = jsv + lensv
      jvt  = ju + lenu
      jscr = jvt + lenvt
      jend = jscr + lenscr
      
      call tstcor (jend-1, lenv, 'frmspovinv')

C     allocate RWF
      
      call conddf (irwspo, ngrid*ngrid)

C     form S(P,Q)
      
      call fileio (2, -irwfac, lenx, v(jx), 0)

      if (iprint .gt. 5)  call outmts (iout,'X:', 0,
     & 0, v(jx), nbas6d, ngrid, nbas6d, ngrid)
      
      call dgemm ( 't', 'n', ngrid, ngrid, nbas6d, one, v(jx), nbas6d,
     &     v(jx), nbas6d, zero, v(js), ngrid )

C     Hadamard square of R^t * R
      do 10 ii = 1, lens
         v(js - 1 + ii) = v(js - 1 + ii)*v(js - 1 + ii)
 10   continue

      if (iprint .gt. 5) call outmts (iout,'S:',0,0,
     &     v(js),ngrid,ngrid,ngrid,ngrid)
      
C     form S(P,Q)^(-1)
C     svd
      call dgesdd( 'A', ngrid, ngrid, v(js), ngrid,
     &     v(jsv), v(ju), ngrid, v(jvt), ngrid, v(jend), lenv-jend,
     &     v(jscr), info )

C      if (info .ne. 0) call gauerr('LAPACK SVD failed in frmspov')

      do 20 ii = 1, lensv
          if ( v(jsv - 1 + ii) .gt. tol ) then
             v(jsv - 1 + ii) = one / v(jsv - 1 + ii)
          else
             v(jsv - 1 + ii) = zero
          endif
  20   continue

C     s*V^t
      call ascaln(ngrid,ngrid, v(jsv), v(ju))

C     assemble back
      call dgemm ( 't', 't', ngrid, ngrid, ngrid, one, v(jvt), ngrid,
     &     v(ju), ngrid, zero, v(js), ngrid )
      
      if (iprint .gt. 5)  call outmts (iout,'S^-1:', 0,
     & 0, v(js), ngrid, ngrid, ngrid, ngrid)

C     save S^-1
      call fileio (1, -irwspo, lens, v(js), 0) 
      
      return
      end
