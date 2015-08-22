      subroutine frmextquad (iout, iprint, filenm, irwfac, nbas6d, ngrid,
     &     v, lenv)
      implicit none

C     Forms quadrature using weights and grid points supplied externally
C     Arguments:
C     [in]  filenm    - filename with grid points and weigths in format (4E25.16) 
C                 first 3 columns are X,Y,Z coordinates, last column is weights
C     [in]  irwfac      - file unit to output the quadrature
C     [in]  nbas6d    - number of basis functions (Cartesian)
C     [out] ngrid     - number of grid points found in file
C
      
C     Arguments:
      integer    iout, iprint, lenv
      real*8     v(*)
      character  filenm*(*)
      integer    irwfac, nbas6d, ngrid

C     Other:
      integer     jend, nval, ii
      integer     jgrid, jweight, jchi, jtmp
      real*8      xx(1)
      
      jend = 1
      call readfmt(iout, 'grid.dat', '(4E25.16)', 4, nval,
     &     v(jend), lenv)
      ngrid   = nval / 4

C     Allocation
      jgrid   = jend
      jweight = jend + 3*ngrid
      jchi    = jweight + ngrid
      jtmp    = jchi    + ngrid * nbas6d
      jend    = jtmp    + ngrid * nbas6d
      call tstcor (jend, lenv, 'custquad: allocation')

C     Reset file size
      call conddf(irwfac,ngrid*nbas6d)
      
      call trcor1(1,1,ngrid,4,3, v(jgrid), xx, xx, v(jtmp), xx, xx,
     &     .false., xx, xx )
      call amove(3*ngrid, v(jtmp), v(jgrid))
      call amove(ngrid, v(jtmp + 3*ngrid), v(jweight))

      if (iprint .gt. 3) then
         write(iout, *) 'Read ngrid: ', ngrid
         write(iout, *) 'Echo first 5 values of grid and weight:'
         write(iout,'(1X,4E24.16)')
     &        (v(jgrid+ii),v(jgrid+ngrid+ii),v(jgrid+2*ngrid+ii),
     &        v(jweight+ii),ii=0,4)
      
         write(iout, *) 'Echo last 5 values of grid:'
         write(iout,'(1X,3E24.16)')
     &        (v(jgrid+ngrid-5+ii), v(jgrid+2*ngrid-5+ii),
     &        v(jgrid+3*ngrid-5+ii), v(jweight+ngrid-5+ii),ii=0,4)
      endif
      
      call evalbf(iout,iprint,1,ngrid,nbas6d, v(jgrid),v(jchi),xx,xx,xx,
     &     v(jend), v(jend), lenv - jend)
      do 900 ii = 0, nbas6d - 1
         call amul(ngrid, v(jweight),v(jchi+ii*ngrid),v(jchi+ii*ngrid))
 900  continue

      call trcor1 (1,1,nbas6d,ngrid,3, v(jchi), xx, xx, v(jtmp), xx, xx,
     &     .false., xx, xx )
      call fileio (1, -irwfac, nbas6d*ngrid, v(jtmp), 0)

      end
