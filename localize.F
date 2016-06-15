      Subroutine Localize(IOut,IPrint,IOpCl,ILocal,MethLc,UnitS,NAtoms,
     $  NBsUse,NBasis,CA,CB,IAtom,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     This is a trimmed adaptation of the Local1 subroutine from l401
C     It does localization according to the Methlc
C     Separate or joint localization for occ/virtual is selected by ilocal
C     
C     If units is set unit overlap is assumed. Otherwise AO overlap is used.
C
C     Ilocal option:   1 - localize occ/virtual separately. Size of occupieds
C                          is set to NAE and NBE
C                      2 - localize occupied/virtuals together
C      
C     Methlc option:   0 - Boys
C                      1 - Population   
C                      2 - Boys + Population
C     IAtom is a scratch array used to form a map from basis functions to atoms
C     Its maximum size is not figured out, but assume natom*nbasis
C      
#include "commonb.inc"
#include "commonmol.inc"
      
      Parameter (IRwS=514,IODipX=518)
      Logical UnitS, CBfn, Spinor
      Dimension CA(NBasis,NBasis),CB(NBasis,NBasis),MOUse(NBasis),
     $  IAtom(*),V(*)
 1000 Format(' Localizing orbitals',10I4,/,
     $      ('                    ',10I4))
      Save One
      Data One/1.0d0/

C     Read /mol/
      call rwmol(2,irwmol)

C     
C     Allocate space.
C
      Call DecoSc(IOpCl,NSpBlk,NRI,NDimBl,CBfn,NSpBlX,Spinor)
      IUHF = NSpBlk - 1
      If(NRI.eq.2.or.NDimBl.eq.2.or.CBfn)
     $  Call GauErr('Missing non-real code in Local1.')
      MethL0 = Mod(MethLc,10)
      If(MethL0.eq.0) then
        NDMAO = 3
        NRIJ = 3
      else if(MethL0.eq.1) then
        NDMAO = 1
        NRIJ = NAtoms
      else if(MethL0.eq.2) then
        NDMAO = 4
        NRIJ = 3 + NAtoms
      else
        Call GauErr('MethLc out of range in Local1.')
        endIf
      NTT = (NBasis*(NBasis+1))/2
      ICL = 1
      IDMAO = ICL + NBasis*NBasis
      IIOrd = IDMAO + NDMAO*NTT
      IIIR = IIOrd + InToWP(NBasis)
      IRIJ = IIIR + InToWP(NBasis)
      IQPIX = IRIJ + NRIJ*NTT
      IEnd = IQPIX + 2*NBasis - 1
      Call TstCor(IEnd,MDV,'Local1')
      If(MethL0.eq.1.or.MethL0.eq.2) then
        If(UnitS) then
          Call AClear(NTT,V(I DMAO))
          Call APlusI(NBasis,One,V(IDMAO))
        else
          Call FileIO(2,-IRwS,NTT,V(IDMAO),0)
          endIf
        IndDip = NTT
      else
        IndDip = 0
        endIf
      If(MethL0.eq.0.or.MethL0.eq.2)
     $  Call FileIO(2,-IODipX,3*NTT,V(IDMAO+IndDip),0)

      If(ILocal.eq.1) then
         NOAUse = NAE
         NOBUse = NBE
      else if(ILocal.eq.2) then
         NOAUse = NBsUse
         NOBUse = NBsUse
      else
        Call GauErr('ILocal out of range in Local1.')
        endIf

      Call RdWrB(3,0)
      Call MakMap(I5DB1,I7FB1,NShell,ShellT,ShellC,JAn,IAtom)
C
C     Localize alpha occupieds.
C
      Call FixPha(1,NBasis,NBasis,0,CA)
      If(NOAUse.gt.0) then
        Do 20 I = 1, NOAUse
   20     MOUse(I) = I
        If(IPrint.gt.1) Write(IOut,1000) (MOUse(I),I=1,NOAUse)
        Call LocMO(IOut,IPrint,MethLc,.True.,.False.,0,NOAUse,NBasis,
     $    NAtoms,NBasis,NTT,MOUse,IAtom,CA,V(ICL),V(IDMAO),V(IIOrd),
     $    V(IIIR),V(IRIJ),V(IQPIX),V,V,V,V,V,V)
        endIf
C
C     Localize alpha virtuals.
C
      NVA = NBsUse - NOAUse
      If(NVA.gt.0) then
        Do 30 I = 1, NVA
   30     MOUse(I) = I + NOAUse
        If(IPrint.gt.1) Write(IOut,1000) (MOUse(I),I=1,NVA)
        Call LocMO(IOut,IPrint,MethLc,.True.,.False.,0,NVA,NBasis,
     $    NAtoms,NBasis,NTT,MOUse,IAtom,CA,V(ICL),V(IDMAO),V(IIOrd),
     $    V(IIIR),V(IRIJ),V(IQPIX),V,V,V,V,V,V)
        endIf
C
C     Localize beta occupieds.
C
      If(IUHF.ne.0) then
        Call FixPha(1,NBasis,NBasis,0,CB)
        If(NOBUse.gt.0) then
          Do 40 I = 1, NOBUse
   40       MOUse(I) = I
          If(IPrint.gt.1) Write(IOut,1000) (MOUse(I),I=1,NOBUse)
          Call LocMO(IOut,IPrint,MethLc,.True.,.False.,0,NOBUse,NBasis,
     $      NAtoms,NBasis,NTT,MOUse,IAtom,CB,V(ICL),V(IDMAO),V(IIOrd),
     $      V(IIIR),V(IRIJ),V(IQPIX),V,V,V,V,V,V)
          endIf
C
C     Localize beta virtuals.
C
        NVB = NBsUse - NOBUse
        If(NVB.gt.0) then
          Do 50 I = 1, NVB
   50       MOUse(I) = I + NOBUse
          If(IPrint.gt.1) Write(IOut,1000) (MOUse(I),I=1,NVB)
          Call LocMO(IOut,IPrint,MethLc,.True.,.False.,0,NVB,NBasis,
     $      NAtoms,NBasis,NTT,MOUse,IAtom,CB,V(ICL),V(IDMAO),V(IIOrd),
     $      V(IIIR),V(IRIJ),V(IQPIX),V,V,V,V,V,V)
          endIf
        endIf
      Return
      End
