      subroutine moquadwrt(iout,iprint,ipart,npart,iatom,natoms,nbas6d,
     $  useocc,nmo,ngrid,nsigbf,sigbfn,occmo,intvec,wpvec,ian,iattyp,
     $  atmchg,c,grid,weight,chi,d1chi,d2chi,d3chi,valmo,valdmo,lenscr,
     $  scr)
      implicit none
C
C     Saves quadrature X of transformed atomic functions to the RWF, which is supplied
C     as the first element of intvec
C
      integer  iout, iprint, ipart, npart, iatom, natoms
      integer  nbas6d, nmo, ngrid, nsigbf, intvec, ian, iattyp
      integer  lenscr
      
      real*8  occmo, wpvec, atmchg, c, grid, weight, chi, d1chi, d2chi
      real*8  d3chi, valmo, valdmo, scr
      
      logical useocc
      integer sigbfn(*)
      dimension occmo(*),intvec(*),wpvec(*),ian(*),iattyp(*),atmchg(*),
     $  c(3,*),grid(ngrid,*),weight(*),chi(ngrid,*),d1chi(ngrid,3,*),
     $  d2chi(ngrid,6,*),d3chi(ngrid,10,*),valmo(ngrid,*),
     $  valdmo(ngrid,3,*),scr(*)

C     Additional variables

      integer    irwfac
      integer    n, ig

C     Read/write the input/output array
      
      irwfac    = intvec(1) 
      intvec(2) = intvec(2) + ngrid

      if (lenscr .lt. nmo ) then
         write (iout, *) lenscr, nmo, ngrid
         call gauerr ('moquadwrt: too small buffer passed to moquadwrt')
      endif

      if (iprint .gt. 5) then
         call outmts(iout,'valMOs:',0,0,valmo,ngrid,nmo,ngrid,nmo)
         call outmts(iout,'weights:',0,0,weight,ngrid,1,ngrid,1)
      endif

      call aclear(ngrid*nmo,scr)
      
      do 110 n = 1, nmo
         do 120 ig = 1, ngrid         
            
            scr((ig-1)*nmo + n) =
     &  dsqrt(dsqrt(weight(ig))) * valmo(ig,n)
            
 120     continue
 110  continue
      
      call fileio (1, irwfac, nmo*ngrid , scr, 0) 

      return
      end
