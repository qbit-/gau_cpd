      subroutine outoctfmt (iout, filenm, mode, name, 
     &     ndims, dims, m)
      implicit none
C     Writes a tensor using an Octave simple io text format to file filenm
C     
C     Arguments:
C     
C     filenm - file name to save to
C     mode   - new file/rewrite: 1, append: 2
C     name   - name of the variable
C     type   - type of data
C     ndims  - number of dimensions
C     dims   - dimension sizes
C     m      - object to save
C     

C     Arguments:

      integer     iout

      integer     mode
      integer     ndims, dims(*)
      real*8      m(*)

      character   filenm*(*), name*(*)
      
C     Other

      integer     lstr1, lstr2
      integer     ii, len, fhandle
      character   blank*80
      logical     exist
C     functions
      integer     linend
      
 1000 format ('# Created by Octave ')
 1001 format ('# name: ', A)
 1002 format ('# type: matrix')
 1003 format ('# ndims: ', I3)
 1004 format (' ',I6,$)
 1005 format (' ',E24.16E3)
 1010 format ('# type: scalar')
 1020 format ('# rows: ', I3)
 1030 format ('# columns: ', I3)
      
      blank = ' '
      lstr1 = Max(LinEnd(filenm),1)
      lstr2 = Max(LinEnd(name),1)

      len   = 1
      do ii = 1, ndims
         len = len*dims(ii)
      enddo
      
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

      if (mode .eq. 1 .or. .not. exist) write (fhandle, 1000)
      
      write (fhandle, 1001) name(1:lstr2)
      if (ndims .eq. 1) then
         if (dims(1) .eq. 1) then
            write (fhandle, 1010)
         else
            write (fhandle, 1002)
            write (fhandle, 1020) dims(1)
            write (fhandle, 1030) 1
         endif
c$$$         Broken - octave reads by lines in a C-style, and then only reshapes data
c$$$         to specified dimension. Therefore we have to use general format even for matrices
c$$$      elseif (ndims .eq. 2) then
c$$$            write (fhandle, 1002)
c$$$            write (fhandle, 1020) dims(1)
c$$$            write (fhandle, 1030) dims(2)         
      else      
         write (fhandle, 1002)
         write (fhandle, 1003) ndims
         write (fhandle, 1004) (dims(ii), ii=1,ndims)
         write (fhandle, *) ''
      endif
      
      write (fhandle, 1005) (m(ii), ii=1,len)
      write (fhandle, *) ''
      write (fhandle, *) ''
      write (fhandle, *) ''

      close (fhandle)
      return

 2000 write (iout, *) 'outoctfmt: File access error!'
      
      end
