      subroutine updemt (iout, iprint, ngrid, p, nump, q, numq, 
     &             ebuf, irwe, v, lenv)
      implicit none
      
C     Updates the E matrix stored in the RWF irwe adding a block 
C     ebuf

C     Arguments
      integer      iout, iprint, lenv      
      real*8       v(*)
      
      integer      ngrid, p, nump, q, numq, irwe
      real*8       ebuf(nump, *)

C     other
      integer      jend, jscr, iq
      integer      ioff

      jend = 1
      jscr = jend
      jend = jscr + nump

      call tstcor (jend-1, lenv)

      do 100 iq = 0, numq - 1
         ioff = p + (q + iq)*ngrid
         call fileio (2, -irwe, nump, v(jscr), ioff)
         call daxpy (nump, 1.0d0, ebuf(1,iq+1), 1, v(jscr), 1)
         call fileio (1, -irwe, nump, v(jscr), ioff)
 100  continue

      end
