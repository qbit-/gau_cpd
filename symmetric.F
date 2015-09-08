      subroutine symmetric (iout, init, isize, kind, m)
      implicit none
C     Returns symmetric or antisymmetric matrices with random elements
C     
C     Arguments:
C
C     init   - if we initialize the RNG in this run
C     isize  - size
C     kind   - 0 - symmetric, 1 - antisymmetric with negative lower part, 2 - antisymmetric with negative upper part
C     m      - square toeplitz matrix
C     

C     Arguments:

      logical     init 
      integer     iout

      integer     isize, kind
      real*8      m(*)

C     functions
      external    srand
      real*8      rand

C     Other
      integer     ii,jj      
      real*8      c,d

      if (kind .eq. 0) then
         c = 1.0d0
         d = 1.0d0
      elseif (kind .eq. 1) then
         c =  1.0d0
         d = -1.0d0
      elseif (kind .eq. 2) then
         c = -1.0d0
         d = 1.0d0
      endif

      if (init) call srand(0)

      call aclear (isize*isize, m)
      
      do 10 ii = 0, isize-1
         do 20 jj = 0, isize-1
            if ( ii .le. jj ) then
               m(1 + ii + jj*isize) = dble(rand(0))*c
            else
               m(1 + ii + jj*isize) = m(ii*isize + jj + 1)*d
            endif
 20      continue
 10   continue

      if (kind .ne. 0) then
         do 30 ii = 0, isize - 1
            m(1 + ii + ii*isize) = 0.0d0            
 30      continue
      endif
      
      return 
      
      end
