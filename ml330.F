
      program ML330

      implicit none

      real*8      work(1)

      integer     in, iout, ipunch
      common      /io/ in, iout, ipunch

      common      /gxwork/ work

      integer*8   ioff

      integer     istat, mdv
      integer     gsetjmp

      in = 5
      iout = 6
      ipunch = 7
      istat = gsetjmp (0)

      if ( istat .eq. 0 ) then
        call initscm (1, 0, 0, work, ioff, mdv)
        call dump_all (work(ioff+1), mdv)
      else
        call prtstat (istat, 'ML325')
      endif


      end


