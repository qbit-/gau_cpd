      subroutine frmspovinv_blas (iout, iprint, nbasis, ngrid, itol,
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
C     [out] nrank - rank of the pseudoinverse
C     [out] irwspo - RWF with S^-1 

      integer    iout, iprint, nbasis, ngrid, itol, nrank
      integer    irwfac, irwspo, lenv
      real*8     v
      dimension  v(*)

      real*8     one, zero
      data       one/1.0d0/,zero/0.0d0/
      
C     other variables
      integer    jend, jx, js, jsv, ju, jvt, jscr, ii
      integer    lenspk
      integer    lenx, lens, lensv, lenu, lenvt, lenscr
      integer    info
      real*8     tol

C     Functions
      real*8     gfloat
      
 1000 format(' Eigenvalues of S:',/,10(1PD10.2))
 1020 format(' SDV returned info=',I12,' to frmspovinv.')
 1100  format (' Overlap: min. sing. value = ', E20.12,
     & '  Сond. number = ', E20.12)
 1110  format ('          rank = ', I8)

      
C     define sizes
      lenx = nbasis * ngrid
      lens = ngrid*ngrid
      lenu = lens
      lenvt = lens
      lensv = ngrid
      lenscr = 9*ngrid
      lenspk = ngrid*(ngrid+1)/2
      
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
     & 0, v(jx), nbasis, ngrid, nbasis, ngrid)
      
      call dgemm ( 't', 'n', ngrid, ngrid, nbasis, one, v(jx), nbasis,
     &     v(jx), nbasis, zero, v(js), ngrid )

C     Hadamard square of R^t * R
      do 10 ii = 0, lens - 1
         v(js + ii) = v(js + ii)*v(js + ii)
 10   continue

      if (iprint .gt. 5) call outmts (iout,'S:',0,0,
     &     v(js),ngrid,ngrid,ngrid,ngrid)
      
C     form S(P,Q)^(-1)
C     svd
      info = 0
      call dgesvd( 'A','A', ngrid, ngrid, v(js), ngrid,
     &     v(jsv), v(ju), ngrid, v(jvt), ngrid, v(jend),
     & lenv - jend, info )

      if ( info .ne. 0 ) then
         write (iout,1020) info
         call gauerr('SVD failed in frmspovinv.')
      endIf
      
      if (iprint .ge. 2) then
         write (iout, 1000) (v(jsv+ii), ii=0,lensv-1) 
      endif
      
      nrank = 0
      tol   = Gfloat(10)**(-itol)
      do 20 ii = 0, lensv - 1 
          if ( abs(v(jsv + ii)) .gt. tol ) then
             v(jsv + ii) = one / v(jsv + ii)
             nrank = nrank + 1
          else
             v(jsv + ii) = zero
          endif
  20   continue

C     s*V^t
      call ascaln(ngrid,ngrid, v(jsv), v(ju))

C     assemble back
      call dgemm ( 't', 't', ngrid, ngrid, ngrid, one, v(jvt), ngrid,
     &     v(ju), ngrid, zero, v(js), ngrid )
      
      if (iprint .gt. 5)  call outmts (iout,'S^-1:', 0,
     & 0, v(js), ngrid, ngrid, ngrid, ngrid)

      call linear (v(js),v(js), ngrid, ngrid)
C     print summary
      write (iout, 1100) one / v(jsv+nrank-1), v(jsv+nrank-1) / v(jsv)
      write (iout, 1110) nrank
      
C     save S^-1
      call fileio (1, -irwspo, lenspk, v(js), 0) 
      
      return
      end
