      subroutine frmrib (iout, iprint, next, irwrib,
     &     keep, ipflag, allowp, accdes, ndbf, ndbf6d,
     &     nbas6d, natoms, ian, atmchg, c, iattyp, v, mdv)
      Implicit Real*8 (A-H,O-Z)
C
C     Form the RI integrals file or place into memory
C     Based on rifrmb
C
C     Arguments:
C     
C     [out]  next            - offset to the next free memory
C     [out]  irwrib          - RWF with the DBF integrals (initialized here)
C     [in]   keep            - if we want to keep the integrals in memory instead of RWF
C     [in]   ipflag, allowp  - usual PRISM variables
C     [in]   natoms, ian,    - usual components of /mol/
C            atmchg, c,
C            iattyp 
      
#include "commonb2.inc"
      Parameter (IHMeth=0)
      Logical AllowP(*)
      Integer IAn(*),IAtTyp(*),IPBC(1),JJ(1)
      Real*8 AtmChg(*),C(*),V(*),XX(1)
      Save Zero,JJ,XX,IPBC
      Data Zero/0.0d0/,JJ/0/,XX/0.0d0/,IPBC/0/

      Logical dbg
      Integer Dims(3)

 1000 Format(' RIFrmB:  LenRWB=',i12)
C
C     Compute lengths
      NBSq6D = NBas6D*NBas6D
      NTT6D  = (NBas6D*(NBas6D+1))/2
C
C     Allocate memory.
      next   = 0
      jBAO   = 1
      jXTmp  = jBAO  + NDBF6D*NBSq6D
      jEnd   = jXTmp + NDBF6D*NBSq6D

      Call TstCor(jEnd,MDV,'frmrib')

      if ( keep ) then
         next = NDBF6D*NBSq6D
      endif
      
      jV = jXTmp
      MDV1 = MDV - jV + 1
C     Form (Q|m,n).
      ICntrl = 600
      I1Cent = 5
      IOpCFC = 7
      NBTI = 0
      IRaf = 1
      Call GenAOS(1,1,NShelB,ShelTB,ShelCB,Junk,AOSB)
      Call FoFCou(IOut,IPrint,IHMeth,IOpCFC,ICntrl,IRaf,0,XX,.False.,0,
     $  0,0,IPFlag,AllowP,0,XX,0,0,.True.,.True.,.False.,.False.,
     $  .False.,I1Cent,AccDes,Zero,0,1,NDBF6D,0,1,0,0,0,0,1,1,JJ,JJ,JJ,
     $  XX,JJ,NBTI,JJ,XX,XX,XX,XX,XX,XX,XX,XX,XX,JJ,XX,XX,JJ,XX,XX,
     $  NAtoms,IAn,AtmChg,C,IAtTyp,JJ,XX,XX,IPBC,XX,Energy,Zero,
     $  V(jBAO),0,0,XX,XX,JJ,JJ,XX,V(jV),V(jV),MDV1)
C     Fix AOS in /B/ as FoFCou converted it back to pure.
      Call GenAOS(1,1,NShell,ShellT,ShellC,Junk,AOS)
      Call GenAOS(1,1,NShelB,ShelTB,ShelCB,Junk,AOSB)
      Call FixAOD(.True.,NShell,ShellT,ShellC,NBas6D,NDBF6D,V(jBAO))
      Call GenAOS(I5DB2,I7FB2,NShelB,ShelTB,ShelCB,Junk,AOSB)
      Call ReDoC1(IOut,IPrint,.False.,1,1,I5DB2,I7FB2,NDBF6D,NDBF,
     $  NDBF6D,NBSq6D,NShelB,ShelTB,ShelCB,V(jBAO),V(jXTmp))
      Call Pcck(0,Junk,V(jBAO),NDBF6D,NBSq6D,NDBF,NBSq6D)
      If(IPrint.ge.8) Call OutMtS(IOut,'(P|m,n) integrals:',0,1,V(jBAO),
     $  NDBF,NBSq6D,NDBF,NBSq6D)

C      (P|m,n) -> (m,n|P)
      Call TrCor1(1,1,NBSq6D,NDBF,3,V(jBAO),XX,XX,V(jXTmp),XX,XX,
     $     .False.,XX,XX)
      Call Amove(NDBF*NBSq6D, V(jXTmp), V(jBAO))
      
C     Allocate the B file and save integrals if requested
      if(.not. keep) then
         LenRWB = NDBF*NBSq6D
         IRwRIB = IGetRW(LenRWB,0)
         Call FileIO(1,-irwrib,NDBF*NBSq6D,V(jBAO),0)
      endif

      Dbg = .False.
      if (Dbg) then
         Dims(1) = nbas6d
         Dims(2) = nbas6d
         Dims(3) = ndbf
         call outoctfmt (iout, 'intsri.mat', 1, 'intsri', 3, dims,
     $    v(jBAO))
      endif
      
      Return
      End
