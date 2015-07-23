      subroutine updmatf (iout, iprint, dimx, x, numx, dimy, y, numy, 
     &             buf, irw, imode, v, lenv)
      implicit none
      
C     Updates a rectangular matrix dimx by dimy stored in the RWF irw by rewriting/adding a block 
C     buf starting at positions x and y
C
C     Arguments:
C                  dimx, dimy - dimenstions of the matrix in the rwf
C                  x, y       - coordinates of the left upper corner of the block
C                  numx, numy - number of elements to update
C                  buf        - rectangular buffer numx by numy containing an update
C                  irw        - RWF to update
C                  imode      - what to do:  0 - write buffer over, 1 - add buffer to what was stored 
C

C     Arguments
      integer      iout, iprint, lenv      
      real*8       v(*)
      
      integer      dimx, dimy, x, numx, y, numy
      integer      imode, irw
      real*8       buf(*)

C     other
      integer      jend, jscr, ii
      integer      ioff

      if ((x + numx .gt. dimx) .or. (y + numy .gt. dimy))
     &     call gauerr ('updmatf: size mismatch')

      if (imode .eq. 0) then 
         do 100 ii = 0, numy - 1
            ioff = x + (y + ii)*dimx
            call fileio (1, -irw, numx, buf(1 + ii*numx), ioff)
 100     continue

      elseif (imode .eq. 1) then 
         jend = 1
         jscr = jend
         jend = jscr + numy
         
         call tstcor (jend-1, lenv, 'updmatf')

         do 200 ii = 0, numy - 1
            ioff = x + (y + ii)*dimx
            call fileio (2, -irw, numx, v(jscr), ioff)
            call daxpy (numx, 1.0d0, buf(1 + ii*numx), 1, v(jscr), 1)
            call fileio (1, -irw, numx, v(jscr), ioff)
 200     continue
      endif
      end
