      Subroutine atquadwrt(IOut,IPrint,IPart,NPart,IAtom,NAtoms,NBas6D,
     $  NGrid,NSigBf,SigBfn,IntVec,WPVec,IAn,IAtTyp,AtmChg,C,Grid,
     $  Weight,Chi,D1Chi,D2Chi,D3Chi,LenScr,Scr)
      Implicit Real*8(A-H,O-Z)
C
C     Saves quadrature X of atomic functions to the RWF, which is supplied
C     as the first element of intvec
C
      Integer SigBfn(*)
      Dimension IntVec(*),WPVec(*),IAn(*),IAtTyp(*),AtmChg(*),C(3,*),
     $  Grid(NGrid,*),Weight(*),Chi(NGrid,*),D1Chi(NGrid,3,*),
     $  D2Chi(NGrid,6,*),D3Chi(NGrid,10,*),Scr(*)


C     additional variables

      integer    irwfac
      integer    mu, ig
      integer    blksz, iblk, ib 

C     Read/write the input/output array
      
      irwfac    = intvec(1) 
      intvec(2) = intvec(2) + ngrid

C     uncompress and write X to the RWF column by column. X comes compressed and transposed
C     in Chi

      blksz = lenscr / nbas6d

      do 10 ib = 1, ngrid, blksz
         call aclear(scr(1),nbas6d*blksz)
         
         do 20 mu1 = 1, nsigbf
            do 30 ig = 1, blksz         
               mu = sigbfn(mu1)
               
               scr((ig-1)*nbas6d + mu) = chi((ib-1)*blksz + ig,mu1)
               
 30         continue
 20      continue
         
         call fileio (2, -irwfac, nbas6d*blksz , scr(1), 0) 
 10   continue
      
      return
      end


