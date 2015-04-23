      subroutine frmria (iout, iprint, next, irwria, ndbf, ndbf6d, keep,
     &     ir12, tolinv, ipflag, allowp, natoms, ian, atmchg, c,
     &     iattyp, v, mdv)
      Implicit Real*8 (A-H,O-Z)
C
C     Form the generalized inverse of the resolution matrix for RI-MP2.
C     Based on the rifrma from l906
C
C     Arguments:
C     
C     [out]  next            - offset to the next free memory
C     [out]  irwria          - RWF with the inverse of DBF overlap (initialized here)
C     [out]  ndbf, ndbf6d    - DBF basis sizes
C     [in]   keep            - if we want to keep the inverse in memory instead of RWF
C     [in]   tolinv          - tolerance for the inverse 
C     [in]   ir12            - metric used:
C                                 1 - 1/r, 2 - overlap
C     [in]   ipflag, allowp  - usual PRISM variables
C     [in]   natoms, ian,    - usual components of /mol/
C            atmchg, c,
C            iattyp
      
#include "commonb2.inc"
      Parameter(IHMeth=0)
      Logical DoAdjA,DoPurA,AllowP(*)
      Integer IAn(*),IAtTyp(*),JJ(1)
      Real*8 AtmChg(*),C(*),V(*),XX(1)
      Save JJ,XX,IJunk,Zero
      Data JJ/0/,XX/0/,IJunk/0/,Zero/0.0d0/
 1000 Format(' RIFrmA: NDBF=',i6,' NDBF6D=',i6,' NRankA=',i6,' DoPurA=',
     $  l1)
C
C     Load /B2/ with the density basis.
      Call RdWrB(-10,0)
C     Figure out NDBF and NDBF6D.
      Call GetMBN(I5DB2,I7FB2,NShelB,ShelTB,ShelCB,NDBF)
      Call GetMBN(1,1,NShelB,ShelTB,ShelCB,NDBF6D)
      Call GenAOS(1,1,NShelB,ShelTB,ShelCB,Junk,AOSB)
C     Compute lenghts and allocate memory.
      NDB6Sq = NDBF6D*NDBF6D
      NDB6TT = (NDBF6D*(NDBF6D+1))/2
      next = 0
      jA   = 1
      jASq = jA   + NDB6Sq
      jEig = jASq + NDB6Sq
      jEnd = jEig + NDBF6D - 1
      Call TstCor(jEnd,MDV,'RIFrmA')

      if ( keep ) then
         next = NDB6Sq
      endif
      
      jV = jEnd + 1
      MDV1 = MDV - jV + 1
C     Form (P|Q).
      DoAdjA = .False.
      DoPurA = .True.
      AccInt = Zero
      Call DFFrmA(IOut,IPrint,DoAdjA,DoPurA,IR12,IHMeth,IPFlag,AllowP,
     $  AccInt,NDBF,NDBF6D,0,1,1,JJ,XX,NAtoms,IAn,AtmChg,C,IAtTyp,I5DB2,
     $  I7FB2,MaxTyB,NShelB,ShelTB,ShelCB,IJunk,IJunk,JJ,V(jA),V(jV),
     $  MDV1)
C     Fix AOS in /B/ as FoFPC in FoFCou converted it back to pure.
      Call GenAOS(1,1,NShell,ShellT,ShellC,Junk,AOS)
C     Form (P|Q)**-1.
      If(DoPurA) then
        NDBFU = NDBF
      else
        NDBFU = NDBF6D
        endIf
      Call GenInv(IOut,IPrint,'A(RI)',0,TolInv,NDBFU,V(jA),V(jASq),
     $  V(jEig),.False.,XX,NRankA,V(jV),MDV1)
      If(IPrint.ge.2) Write(IOut,1000) NDBF,NDBF6D,NRankA,DoPurA
C     Expand overlap to full (Fixme: all code has to be moved to lt matrices)
      call square (v(jA),v(jA),NDBFU,NDBFU,0)
C     Save (P|Q)**-1.      
      if (.not. keep) then
C     NDBUTT = (NDBFU*(NDBFU+1))/2
         NDBUSq = NDBFU*NDBFU
         IRwRIA = IGetRW(NDBUSq,0)
         Call FileIO(1,-IRwRIA,NDBUSq,V(jA),0)
      endif
      Return
      End
