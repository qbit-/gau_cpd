      subroutine toeplitz (iout, isize, kind, m)
      implicit none
C     Returns symmetric (special kinds of Toeplitz) matrices with elements from 1 to isize
C     
C     Arguments:
C     
C     isize  - size
C     m      - square toeplitz matrix
C     

C     Arguments:

      integer     iout

      integer     isize, kind
      real*8      m(*)
      
C     Other
      integer     ii,jj      
      real*8      c,d

      if (kind .eq. 0) then
         c = 1.0d0
         d = 1.0d0
      elseif (kind .eq. 1) then
         c = -1.0d0
         d =  1.0d0
      elseif (kind .eq. 2) then
         c = 1.0d0
         d = -1.0d0
      endif
         
      do 10 ii = 0, isize-1
         do 20 jj = 0, isize-1
            if ( ii .gt. jj ) then
               m(1 + ii + jj*isize) = dble(ii - jj + 1)*c
            else
               m(1 + ii + jj*isize) = dble(jj - ii + 1)*d
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
