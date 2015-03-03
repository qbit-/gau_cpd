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
     $  D2Chi(NGrid,6,*),D3Chi(NGrid,10,*),Scr(NGrid,*)
C
C      write (iout, *) 'atquadwrt called. ngrid = ', ngrid
C      intvec(1) = intvec(1) + ngrid
      Return
      End


