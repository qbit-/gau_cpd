      subroutine outcsv (iout, filenm, mode, nelem, m)
      implicit none
C     Writes an array as 1 line in CSV format to file filenm
C     
C     Arguments:
C     
C     filenm - file name to save to
C     mode   - new file/rewrite: 1, append: 2
C     nelem  - length
C     m      - object to save
C     

C     Arguments:

      integer     iout

      integer     mode, nelem
      real*8      m(*)

      character   filenm*(*)
      
C     Other

      integer     lstr1
      integer     ii, fhandle
      character   blank*80
      logical     exist
C     functions
      integer     linend

 1000 format (E20.14,$)
      
      blank = ' '
      lstr1 = Max(LinEnd(filenm),1)

      fhandle = 771
 
      inquire (file=filenm(1:lstr1), exist=exist)
      if (exist) then
         if (mode .eq. 2) then
            open(unit=fhandle,file=filenm(1:lstr1), status='old',
     &           access='append', err = 2000)
         else
            open(unit=fhandle,file=filenm(1:lstr1), status='unknown',
     &           err = 2000)
         endif
      else
         open(unit=fhandle,file=filenm(1:lstr1), status='new',
     &        err = 2000)
      endif
      
      do 100 ii = 1,nelem -1 
         write (fhandle,'(E20.14,$)') m(ii)
         write (fhandle,'(A,$)') ','
 100  continue
      write (fhandle,'E20.14') m(nelem)
      
      close (fhandle)
      return

 2000 write (iout, *) 'outcsv: File access error!'
      
      end
