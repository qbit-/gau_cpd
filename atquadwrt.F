      Subroutine atquadwrt(IOut,IPrint,IPart,NPart,IAtom,NAtoms,NBas6D,
     $  NGrid,NSigBf,SigBfn,IntVec,WPVec,IAn,IAtTyp,AtmChg,C,Grid,
     $  Weight,Chi,D1Chi,D2Chi,D3Chi,LenScr,Scr)
      implicit none
C
C     Saves quadrature X of atomic functions to the RWF, which is supplied
C     as the first element of intvec
C
      integer  iout, iprint, ipart, npart, iatom, natoms
      integer  nbas6d, ngrid, nsigbf, intvec, ian, iattyp, lenscr
      real*8   atmchg, c, wpvec, grid, weight, chi, d1chi
      real*8   d2chi, d3chi, scr
      
      Integer SigBfn(*)
      Dimension IntVec(*),WPVec(*),IAn(*),IAtTyp(*),AtmChg(*),C(3,*),
     $  Grid(NGrid,*),Weight(*),Chi(NGrid,*),D1Chi(NGrid,3,*),
     $  D2Chi(NGrid,6,*),D3Chi(NGrid,10,*),scr(*)


C     additional variables

      integer    irwfac
      integer    mu, mu1, ig
      integer    blksz, nblk, ib 

C     Read/write the input/output array
      
      irwfac    = intvec(1) 
      intvec(2) = intvec(2) + ngrid

      if (iprint .gt. 5) then
         call prtchi(iout,0,ngrid,nbas6d,chi,d1chi,d2chi,d3chi)
         call outmts(iout,'weights:',0,0,weight,ngrid,1,ngrid,1)
      endif
      
C     Take sqrt of weights
c$$$      do 10 ig = 1, ngrid
c$$$         weight(ig) = dsqrt(weight(ig))
c$$$ 10   continue
      
C     uncompress and write X to the RWF column by column. We need to uncompress, transpose
C     and weight Chi to form X

      if (lenscr .lt. nbas6d ) then
         write (iout, *) lenscr, nbas6d, ngrid
         call gauerr ('atquadwrt: too small buffer passed to atquadwrt')
      endif

      call aclear(ngrid*nbas6d,scr)

      do 110 mu1 = 1, nsigbf
         do 120 ig = 1, ngrid         
            
            mu = sigbfn(mu1)
            scr((ig-1)*nbas6d + mu) =
     &  dsqrt(dsqrt(weight(ig))) * chi(ig,mu1)
            
 120     continue
 110  continue
         
      call fileio (1, irwfac, nbas6d*ngrid , scr, 0) 
      
      return
      end


