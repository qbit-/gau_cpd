      subroutine aclearf (n, irw, v, lenv)
      implicit none
      
C     Clears a specified RWF explicitly to overcome a BUG in conddf: there may be garbage in the
C     freshly initialized RWFs at ~524 Kword
C     Arguments:
C                  n          - number of elements to clean
C                  irw        - RWF number
C                  

C     Arguments
      integer      n, irw, lenv      
      real*8       v(*)
      
C     other
      
      integer      jscr, ii
      integer      nblk, blksz

      nblk  = n / lenv

      if (nblk .lt. 1) then
         nblk = 1
         blksz = n
      endif

      if (nblk .ne. 1) then
         blksz = n / nblk
      endif
      
      call aclear(blksz, v(1))

      do 100 ii = 1, nblk
         call fileio (-1, irw, blksz, v(1), 0)
 100  continue

      call fileio (3, irw, 0, v(1), 0)
      
      end
