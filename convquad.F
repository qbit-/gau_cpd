      subroutine convquad (iout, iprint, donorm,
     &     dopure, dotran, ngrid, irwfac, v, lenv)

      implicit none
      
C-----------------------------------------------------------------------------+
C     This routine converts raw atomic quadratures to windowed MO quadratures
C     (possibly complex), also making proper adjustments for 5d/6d basis and
C     PRISM scaling of D/F/.. shells.
C     Arguments:
C     [in]   donorm - if we normalize D/F and higher shell to PRISM normalization
C     [in]   dopure - if we convert to pure D/F functions (NYI)
C     [in]   dotran - if we transform the quadrature to MO basis (NYI)
C     [in]   ngrid - number of grid points in the quadrature 
C    
C     [in/out] irwfac  - RWF with X 
C     
C-----------------------------------------------------------------------------+
#include "commonb.inc"
      
C     arguments
      
      integer     iout, iprint
      integer     ngrid, irwfac
      real*8      v(*)
      integer     lenv
      logical     donorm, dopure, dotran
      
C     quadrature related variables
      
      integer     nbasis, nbas6d
      
C     pointers

      integer     jend, jxbuf, jscr, junk
      integer     nxblk
      real*8      xx
      data        junk/0/, xx/0.0d0/
      
      integer     npmax

C     functions called in quadrature routine
      
      integer     nproc
      
      jend = 1

C     Get basis information. Load /B/

      call rdwrb(2,0)
      call getmbn(I5DB1,I7FB1,NShell,ShellT,ShellC,nbasis)      
      call getmbn(1,1,NShell,ShellT,ShellC,nbas6d)
      
C     memory allocation

      jxbuf = jend
      jscr  = jxbuf + nbas6d*ngrid
      jend  = jscr  + nbas6d
      
      call tstcor (jend-1, lenv, 'convquad')
      call fileio (2, -irwfac,nbas6d*ngrid, v(jxbuf),0)

      if (donorm) then
         call adjmo2 (iout, iprint,.false.,1,1,nbas6d,
     &        ngrid,NShell,ShellT,ShellC,v(jxbuf))
      endif
      
      if (donorm .and. dopure) then
            call genaos(I5DB1, I7FB1, NShell, ShellT, ShellC, junk, AOS)
            call redoc1(iout, iprint, .false., 1,1, I5DB1,I7FB1, nbas6d,
     &           nbasis, nbas6d, ngrid, NShell, ShellT, ShellC,
     &           v(jxbuf),v(jscr))
            call pcck(0, junk, v(jxbuf), nbas6d, ngrid,
     &           nbasis, ngrid)

      endif

      call fileio (1, -irwfac,nbas6d*ngrid, v(jxbuf),0)
           
      end
