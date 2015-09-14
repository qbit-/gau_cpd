      subroutine frmquad (iout, iprint, iradan, iranwt, irangd, iopcl,
     & dotran, numtr, tform, ngridu, irwfac, v ,lenv)

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
C     [in]  dotran  - use full or transformed basis (MOs, localized AOs etc)
C     [in]  numtr   - number of columns in the transformation matrix 
C     [in]  tran    - transformation matrix
C     [out] ngridu  - size of the grid
C     [out] irwfac  - RWF for storing X (initialized here)
C     
C-----------------------------------------------------------------------------+
#include "commonmol.inc"
      
C     arguments
      
      integer     iout, iprint
      integer     iradan, iranwt, irangd
      integer     iopcl, ngridu, irwfac
      integer     numtr
      real*8      v(*), tform(*)
      integer     lenv
      logical     dotran
      
C     quadrature related variables
      
      integer     nbas6d, iopclx
      integer     ngridmx, ntheta, nphi, mxradu, mxangu, iprune
      integer     iopwrt, leniopwrt, lenscrwrt

      logical    parok, pardig, readb, usebfn, usemo, useocc
      logical    hvrmax
      
      parameter  ( leniopwrt = 5 )
      dimension   iopwrt(leniopwrt)

C     pointers

      integer     jend, junk
      real*8      xx
      parameter   ( junk = 0, xx = 0.0d0 )
      
      integer     npmax

C     functions called in quadrature routine
      
      external    atquadwrt, moquadwrt, denfnx 
      integer     nproc
      
      jend = 1

C     Get 6d basis used with all quadratures

      call getnb6(nbas6d)

C     Get the maximal number of grid points
      
      call DecPrn(iout,iprint,iradan,0,iprune,ntheta,nphi,
     &  mxradu,mxangu,xx) 
      ngridmx = mxradu * mxangu * natoms

      if (dotran) then
         call conddf (irwfac, ngridmx*numtr)
      else
         call conddf (irwfac, ngridmx*nbas6d)
      endif
      
C     Prepare the arguments passed/received from the writing routine

      npmax = nproc(0)

C     iopwrt integer array:
C        1 - RWF to store X(mu,G)
C        2 - total number of grid points (incremented by working routine)
      
      iopwrt(1) = irwfac
      iopwrt(2) = 0

C     scratch passed to the writing routine
C     here we just give maximal possible length - may need adjustment in future

      if (dotran) then
         lenscrwrt  = numtr
      else
         lenscrwrt  = nbas6d
      endif

      iopclx = 0
      parok  = .true.
      pardig = .false.
      readb  = .true.
      usebfn = .not.dotran
      usemo  =  dotran
      useocc = .false.
      hvrmax = .false.
      
      call nin3su (iout,iprint,iopclx,parok,pardig,npmax,readb,
     $  usebfn,usemo,useocc,0,0,0,iradan,iranwt,irangd,
     $  2,numtr,natoms,nbas6d,IAn,iattyp,atmchg,c,xx,tform,xx,
     $  xx,xx,hvrmax,xx,atquadwrt,iopwrt,xx,lenscrwrt,moquadwrt,iopwrt,
     $  xx,lenscrwrt,DenFnX,junk,xx,0,v(jend),lenv - jend)

C     Read back the output from the writing routine
C     actual total number of grid points:

      ngridu = iopwrt(2)      
      write (iout, *) 'Frmquad: nin3su completed, npoints = ', ngridu

      end
