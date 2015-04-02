      subroutine frmquad (iout, iprint, iradan, iranwt, irangd, iopcl,
     & ngridu, irwfac, irwfact, v ,lenv)

      implicit none
      
C-----------------------------------------------------------------------------+
C     This routine computes the collocation matrix X(mu,G)
C     and saves it to the RWF irwfac
C
C     Currently, Becke scheme is used to get X.
C
C     Arguments:
C     [in]   iradan   - number of radial/angular points XXXYYY
C            requests XXX radial points with 2d angular grid
C            for possible choises see angquad
C            iradan = -XXXYYY requests XXX radial points a spherical product grid
C            with YYY theta points and 2*YYY phi points.
C     [in]  iranwt, irangd - selects cutoffs for the grid (?)
C     [in]  iopcl   - usual WF selection parameter
C     [out] ngridu  - size of the grid
C     [out] irwfac  - RWF for storing X (initialized here)
C     [out] irwfact - RWF for X transpose (initialized here)
C     
C-----------------------------------------------------------------------------+
#include "commonmol.inc"
      
C     arguments
      
      integer     iout, iprint
      integer     iradan, iranwt, irangd
      integer     iopcl, ngridu, irwfac, irwfact
      real*8      v(*)
      integer     lenv

C     quadrature related variables
      
      integer     nbas6d
      integer     ngridmx, ntheta, nphi, mxradu, mxangu, iprune
      integer     iopwrt, leniopwrt, lenscrwrt

      parameter  ( leniopwrt = 5 )
      dimension   iopwrt(leniopwrt)

C     pointers

      integer     jend, junk
      real*8      xx

      integer     npmax

C     functions called in quadrature routine
      
      external    atquadwrt, mofunx, denfnx 
      integer     nproc
      
C     Get 6d basis used with all quadratures

      call getnb6(nbas6d)

C     Get the maximal number of grid points
      
      call DecPrn(iout,iprint,iradan,0,iprune,ntheta,nphi,
     &  mxradu,mxangu,xx) 
      ngridmx = mxradu * mxangu * natoms

      call conddf (irwfac, ngridmx*nbas6d)

C     Prepare the arguments passed/received from the writing routine

      npmax = nproc(0)

C     iopwrt integer array:
C        1 - RWF to store X(mu,G)
C        2 - total number of grid points (incremented by working routine)
      
      iopwrt(1) = irwfac
      iopwrt(2) = 0

C     scratch passed to the writing routine
C     here we just give maximal possible length - may need adjustment in future
      
      lenscrwrt  = ngridmx*nbas6d

      call NIn3Su(iout,iprint,iopcl,.true.,.false.,npmax,.true.,
     $  .true.,.false.,.false.,0,0,0,iradan,iranwt,irangd,
     $  2,0,natoms,nbas6d,IAn,iattyp,atmchg,c,xx,xx,xx,
     $  xx,xx,.false.,xx,atquadwrt,iopwrt,xx,lenscrwrt,MOFunX,junk,xx,
     $  0,DenFnX,junk,xx,0,v(jend),lenv)

C     Read back the output from the writing routine
C     actual total number of grid points:

      ngridu = iopwrt(2)      
      write (iout, *) 'Frmquad: nin3su completed, npoints = ', ngridu

      call conddf (irwfact, ngridu*nbas6d)
C      call transp (nbasis, ngridu, irwfac, irwfact, v, lenv)

      end
