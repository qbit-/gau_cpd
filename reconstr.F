      subroutine reconstr (iout, iprint, nbasis, ngrid,
     &     irwfac, irwz, irwrec, v, lenv)
C     Reconstruct the tensor from it's decomposition
      implicit none

      integer iout, iprint, lenv
      real*8  v(*)
      
      integer nbasis, ngrid
      integer irwfac, irwfact, irwz, irwrec

      integer jend
      integer jveci, jvecj, jveck, jvecl
      integer jz, jvecij, jveckl, jscr
      integer i, j, k, l, ii

      integer igetrw
      real*8  ddot
      real*8  elem

      integer jx, jxt
      real*8  xx
C
      
      irwfact = igetrw (ngrid*nbasis, 0)
      call mattrn (1,1,ngrid,nbasis,3,irwfac,irwfact,lenv/2,v(1))  
      
C       
      jend   = 1
      jz     = jend
      jveci  = jz + ngrid*ngrid
      jvecj  = jveci  + ngrid
      jveck  = jvecj  + ngrid
      jvecl  = jveck  + ngrid
      jvecij = jvecl  + ngrid
      jveckl = jvecij + ngrid
      jscr   = jveckl + ngrid
      jend   = jscr   + ngrid
      
      call tstcor (jend-1, lenv, 'reconstr')

      irwrec = igetrw (nbasis*nbasis*nbasis*nbasis, 0)
      
C     Read Z
      call fileio(2,-irwz,ngrid*ngrid,v(jz),0)

C     Now the fun starts
      do 100 i = 0, nbasis-1
         call fileio (2, -irwfact, ngrid, v(jveci), i*ngrid)
         do 200 j = 0, nbasis-1
            call fileio (2, -irwfact, ngrid, v(jvecj), j*ngrid)

            do 250 ii = 0, ngrid-1
               v(jvecij + ii) = v(jveci + ii)*v(jvecj + ii)
 250        continue

            do 300 k = 0, nbasis-1
               call fileio (2, -irwfact, ngrid, v(jveck), k*ngrid)
               do 400 l = 0, nbasis-1
                  call fileio (2, -irwfact, ngrid, v(jvecl),
     &                 l*ngrid)

                  do 450 ii = 0, ngrid-1
                     v(jveckl + ii) = v(jveck + ii)*v(jvecl + ii)
 450              continue

                  call dgemv ('n', ngrid, ngrid, 1.0d0, v(jz), ngrid,
     &                 v(jvecij), 1, 0.0d0, v(jscr), 1)
                  elem = ddot (ngrid, v(jscr), 1, v(jveckl), 1)

                  call fileio (1, irwrec, 1, elem, 0)
 400           continue
 300        continue
 200     continue
 100  continue
      
      end
