      subroutine outoctfmt (iout, filenm, name,
     &     ndims, dims, m)
      implicit none
C     Writes a tensor using an Octave simple io text format to file filenm
C     
C     Arguments:
C     
C     filenm - file name to save to 
C     name   - name of the variable
C     ndims  - number of dimensions
C     dims   - dimension sizes
C     m      - object to save
C     

C     Arguments:

      integer     iout
      
      integer     stat, ndims, dims(*)
      real*8      m(*)

      character   filenm*(*), name*(*)
      
C     Other

      integer     lstr1, lstr2
      integer     ii, len, fhandle
      character   blank*80

C     functions
      integer     linend
      
 1000 format ('# Created by Octave ')
 1001 format ('# name: ', A)
 1002 format ('# type: matrix')
 1003 format ('# ndims: ', I3)
 1004 format (' ',I6,$)
 1005 format (' ',E22.16)
      
      blank = ' '
      lstr1 = Max(LinEnd(filenm),1)
      lstr2 = Max(LinEnd(name),1)

      len   = 1
      do ii = 1, ndims
         len = len*dims(ii)
      enddo
      
      fhandle = 771
      open (unit=fhandle,file=filenm(1:lstr1),status='UNKNOWN',
     &     err = 2000)
      
      write (fhandle, 1000)
      write (fhandle, 1001) name(1:lstr2)
      write (fhandle, 1002)
      write (fhandle, 1003) ndims
      write (fhandle, 1004) (dims(ii), ii=1,ndims)
      write (fhandle, *) ''
      write (fhandle, 1005) (m(ii), ii=1,len)
      write (fhandle, *) ''
      write (fhandle, *) ''
      write (fhandle, *) ''

      close (fhandle)
      return

 2000 write (iout, *) 'outoctfmt: File access error!'
      
      end
