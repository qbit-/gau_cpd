      subroutine readfmt (iout, filenm, fmt, irecln, nelem, v, lenv)
      implicit none
C     Writes an array as 1 line in CSV format to file filenm
C     
C     Arguments:
C     
C     filenm - [in]  file name to read from
C     fmt    - [in]  record format
C     irecln - [in]  number of elements in each record. Must match fmt
C     nelem  - [out] length read
C     v      - [out] buffer
C     lenv   - [in]  length of the buffer. If negative, a warning is printed 
C                    on buffer full

C     Arguments:

      integer     iout

      integer     nelem, lenv, irecln
      real*8      v(*)
      character   filenm*(*)
      character   fmt*(*)
      
C     Other

      integer     lstr1, lstr2
      integer     iostat, fhandle
      integer     ii
      character   blank*80
      logical     exist, warn

C     functions
      integer     linend

      if ( lenv .lt. 0 ) then
         lenv = -lenv
         warn = .true.
      endif
      
      blank = ' '
      lstr1 = Max(LinEnd(filenm),1)
      lstr2 = Max(LinEnd(fmt),1)
      
      fhandle = 771
 
      inquire (file=filenm(1:lstr1), exist=exist)
      if (exist) then
         open(unit=fhandle,file=filenm(1:lstr1), status='old',
     &   err = 100)
      else
         write(iout, '(A,A,A)') 'readcsv: file ',
     &        filenm(1:lstr1), ' not found'
         nelem = -1
         return
      endif

      nelem = 0
      do 10 while (nelem .lt. lenv)
         read(fhandle,fmt(1:lstr2),IOSTAT=iostat,err = 100)
     &    (v(nelem + ii + 1), ii = 0, irecln - 1)
         if (iostat .gt. 0) then
            write (iout, '(A, A, I4, A)') 'readcsv: file ',
     &           filenm(1:lstr1),
     &           ' record ', nelem, ' read error'
         else if (iostat .lt. 0) then
            goto 20
         else
            nelem = nelem + irecln
         endif
 10   continue 
 20   continue
      
      if ((nelem .eq. lenv) .and. warn)
     &     write(iout, '(A,A)') 'readfmt: warning: read stopped',
     &     ' because memory exhausted'
      close (fhandle)
      return

 100  continue   
      write (iout, *) 'readcsv: file', filenm(1:lstr1), ' access error!'
      nelem = -1
      return
      
      end
