      subroutine frmspovinv (iout, iprint, nbas6d, ngrid, tol,
     &     irwfac, irwspo, v, lenv)

      implicit none

C     Computes the Moore-Penrose pseudoinverse of the matrix
C     S(P,Q) = ( X(mu,P)^t X(mu,Q) ) o ( X(nu,P)^t X(nu,Q) )
C     X is stored in the irwfac, s is put into irwspo
C
C     Options:
C     nbas6d - cartesian basis size (dimension of mu)
C     ngrid  - number of grid points (dimension of P)
C     tol    - numeric tolerance for the pseudoinverse

      integer    iout, iprint, nbas6d, ngrid
      integer    irwfac, irwspo, lenv
      real*8     tol, v
      dimension  v(*)

      real*8     one, zero
      parameter  ( one = 1.0d0, zero = 0.0d0)
      
C     other variables
      integer    jend, jx, js, jsv, ju, jvt, jscr, ii
      integer    lenx, lens, lensv, lenu, lenvt, lenscr
      integer    info
      
C     define sizes
      lenx = nbas6d * ngrid
      lens = ngrid*ngrid
      lenu = lenx
      lenvt = lenx
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

C     form S(P,Q)
      
      call fileio (2, -irwfac, lenx, v(jx), 0)

      if (iprint .gt. 5)  call outmts (iout,'X:', 0,
     & 0, v(jx), nbas6d, ngrid, nbas6d, ngrid)
      call aclear (v(js), lens)
      
      call dgemm ( 't', 'n', ngrid, ngrid, nbas6d, one, v(jx), nbas6d,
     &     v(jx), nbas6d, zero, v(js), ngrid )

      do 10 ii = 1, lens
         v(js - 1 + ii) = v(js - 1 + ii)*v(js - 1 + ii)
 10   continue
      
C     form S(P,Q)^(-1)
C     svd
      call dgesdd ( 'S', ngrid, ngrid, v(js), ngrid, v(jsv), v(ju),
     &     ngrid, v(jvt), ngrid, v(jend), lenv - jend, v(jscr), info)

      do 20 ii = 1, lensv
         if ( v(jsv - 1 + ii) .gt. tol ) then
            v(jsv - 1 + ii) = one / v(jsv - 1 + ii)
         else
            v(jsv - 1 + ii) = zero
         endif
 20   continue

C     s*V^t
      do 30 ii = 1, lensv
         call dscal(ngrid,v(js-1+ii), v(jvt - 1 + (ii-1)*ngrid), 1)
 30   continue

C     assemble back
      call dgemm ( 'n', 'n', ngrid, ngrid, ngrid, one, v(ju), ngrid,
     &     v(jvt), ngrid, zero, v(js), ngrid )
               
      
      return
      end
