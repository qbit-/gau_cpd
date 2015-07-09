      subroutine frmspovinv (iout, iprint, nbasis, ngrid, tol,
     &     irwfac, nrank, irwspo, v, lenv)

      implicit none

C     Computes the Moore-Penrose pseudoinverse of the matrix
C     S(P,Q) = ( X(mu,P)^t X(mu,Q) ) o ( X(nu,P)^t X(nu,Q) )
C     X is stored in the irwfac, s is put into irwspo
C
C     Options:
C     [in] nbasis - cartesian basis size (dimension of mu)
C     [in] ngrid  - number of grid points (dimension of P)
C     [in] tol    - numeric tolerance for the pseudoinverse
C     [in] irwfac - RWF with quadrature
C     [out] irwspo - RWF with S^-1 

      integer    iout, iprint, nbasis, ngrid
      integer    irwfac, irwspo, lenv
      real*8     tol, v, xx
      dimension  v(*)

      real*8     one, zero, half
      parameter  ( one = 1.0d0, zero = 0.0d0 )
      
C     other variables
      integer    jend, jx, js,  jspk, jsv, ju, jscr, ii
      integer    lenx, lens, lensv, lenu, lenspk, lenscr
      integer    nrank
      
C     define sizes
      
      lenx   = nbasis * ngrid
      lens   = ngrid * ngrid 
      lenspk = ngrid * (ngrid + 1) / 2
      lenu   = lens
      lensv  = ngrid
      
C     allocate memory
      jend = 1
      jx   = jend
      js   = jx + lenx
      ju   = js + lens
      jsv  = ju + lenu
      jend = jsv + lensv
      
      call tstcor (jend-1, lenv, 'frmspovinv')

C     allocate RWF
      
      call conddf (irwspo, ngrid*ngrid)

C     form S(P,Q)
      
      call fileio (2, -irwfac, lenx, v(jx), 0)

      if (iprint .gt. 5)  call outmts (iout,'X:', 0,
     & 0, v(jx), nbasis, ngrid, nbasis, ngrid)
      
      call aclear (lens, v(js))
      call dsyrk ('u', 't', ngrid, nbasis, one, v(jx), nbasis,
     &     zero, v(js), ngrid)

      call linear (v(js), v(js), ngrid, ngrid)
      
C     Hadamard square of R^t * R
      do 10 ii = 0, lenspk - 1
         v(js + ii) = v(js + ii)*v(js + ii)
 10   continue

      if (iprint .gt. 5) call ltoutl (iout,'S:',0, ngrid, 
     &     v(js))
      
C     form S(P,Q)^(-1)

      call geninv(iout,iprint,'S:', 0, tol, ngrid, v(js), v(ju),
     &     v(jsv), .false., xx, nrank, v(jend), lenv - jend)
      
C     save S^-1
      call fileio (1, -irwspo, lenspk, v(js), 0) 
      
      return
      end
