*Deck NumIn3
      Subroutine NumIn3(IOut,IPrint,IOpCl,InPar,ParDig,DoSuper,UseBfn,
     $  UseMO,UseOcc,IUseDn,IDerBf,IDerRh,IRadAn,IRanWt,IRanGd,IRTypI,
     $  NMO,NAtoms,NAtomA,NClRep,NBas6D,NDBF,MapAtm,IAn,IAtTyp,AtmChg,C,
     $  NShell,MaxTyp,C1,C2,C3,C4,EXX,JAn,ShellT,ShellC,ShellA,ShlADF,
     $  ShellN,X,Y,Z,RLam,NShelB,MaxTyB,C1B,C2B,C3B,C4B,EXXB,JAnB,
     $  ShelTB,ShelCB,ShelAB,ShlADB,ShelNB,XB,YB,ZB,Thresh,SizInc,CMO,
     $  OccMO,PA,PB,BegShl,BegShB,AtmSiz,ShlSiz,AtmSzD,ShlSzB,DoInAt,
     $  NDisMt,DisMat,NOpAtG,NEqAtm,RotOp,MaxNGI,NAtG,AtEC,DoAtom,
     $  BasNrm,HvRMax,RQMax,BfnFnc,IntBfn,WPBfn,LScrBf,MOFunc,IntMO,
     $  WPMO,LScrMO,DenFnc,IntDen,WPDen,LScrDn,IPart,NPrtUL,NPart,NTotA,
     $  ElSumD,ElSumM,ElSumA,DynPar,LinDyn,NxtVal,I,V,LenV)
      Implicit Real*8(A-H,O-Z)
C
C     Work routine for numerical integration using Becke's scheme.
C     See NIn3Su for details on the arguments.
C
#include "allxc.inc"
      Character*(*) CountN
      Logical UseBfn,UseMO,UseOcc,DoSuper,DoRho,DoWrt,UseP,DoPhi,DoV,
     $  ScrnBf,Packed,NonLoc,DoIntr,DoTau,DoLap,HvRMax,DynPar,LinDyn,
     $  SSWeig,DubWgt,PrtVal,AtBlock,ScrnGd,UseDWt,DoMicB,DoInAt(*),
     $  DoVarX,DoVarG,DoSpar,StoBW2,DoNuFx,DoCPKS,DoGmPe,DoGmPx,PeGmx,
     $  DoHyp,TDKSLg,DoPolD,DoGxPT,DoGmx,PKSDen,DoNDXC,DoVN2,SkipAt,
     $  UseTI,OneRad,DoTIDE,NoAtTx,PrnSph,DoAtom(*),CBfn,Spinor,PBC,
     $  DoHir,DoAtSD,DoZCmp,ParDig,InPar,NeedCR,UseMOs,UseSig
C     ###%%% Removed dorho from parameters to allow changing it
C      Parameter (DoRho=.True.,DoV=.True.,IfTau=0,IfLap=0,DoVarX=.False.,
C     ###%%%
      Parameter (DoV=.True.,IfTau=0,IfLap=0,DoVarX=.False.,
     $  Packed=.True.,MaxPtI=0,MaxMBI=0,MxDnXC=DEFMXDNXC,NDPol=0,
     $  MxTyXC=DEFMXTYXC,MxVrXC=(MxDnXC*MxTyXC),DoVarG=.False.,
     $  DoNuFx=.False.,DoCPKS=.False.,DoGmPe=.False.,DoGmPx=.False.,
     $  PeGmx=.False.,DoHyp=.False.,TDKSLg=.False.,DoPolD=.False.,
     $  DoGxPT=.False.,DoGmx=.False.,PKSDen=.False.,DoNDXC=.False.,
     $  AtBlock=.False.,IDerNu=0,NAtCel=0,DoVN2=.False.,DoSpar=.False.,
     $  StoBW2=.False.,NMtPBC=1,OneRad=.False.,MaxAn=200,PBC=.False.,
     $  DoZCmp=.True.,CountN='NumIn3 -- NxtVal')
      Logical DidAtm(0:MaxAn),LBit
      External BfnFnc,MOFunc,DenFnc
      Dimension IAn(*),IAtTyp(*),AtmChg(*),C(3,*),C1(*),C2(*),C3(*),
     $  C4(*),EXX(*),JAn(*),CMO(*),OccMO(*),PA(*),PB(*),IntBfn(*),
     $  WPBfn(*),IntMO(*),WPMO(*),IntDen(*),WPDen(*),V(*),I(*),JJ(1),
     $  ShlSzB(*),RotOp(*),XX(1),CntCel(3),NAtG(*),DisMat(*),AtEC(*),
     $  X(*),Y(*),Z(*),RQMax(*),AtmSiz(*),ShlSiz(*),NEqAtm(*),MapAtm(*),
     $  C1B(*),C2B(*),C3B(*),C4B(*),EXXB(*),JAnB(*),RLam(*),XB(*),YB(*),
     $  ZB(*),AtmSzD(*),BasNrm(*)
      Integer ShellT(*),ShellC(*),ShellA(*),ShlADF(*),ShellN(*),
     $  BckCel(1),BckDim(3),IntCel(3),IVDen(MxVrXC),IVTyp(MxVrXC),
     $  IVC2T(MxVrXC),IVT2C(MxVrXC),IVLoc(MxVrXC),IUDen(MxVrXC),
     $  IUTyp(MxVrXC),IUC2T(MxVrXC),IUT2C(MxVrXC),IULoc(MxVrXC),GLinCo,
     $  BegShl(*),BegShB(*),ShelTB(*),ShelCB(*),ShelAB(*),ShlADB(*),
     $  ShelNB(*)
      Save Zero,JJ,XX,CntCel,BckCel,BckDim,IntCel,One,Two
      Data Zero/0.0d0/,JJ/0/,One/1.0d0/,XX/0.0d0/,RInc/1.0001d0/,
     $  CntCel/3*0.0d0/,BckCel/1/,BckDim/3*0/,IntCel/3*0/,Two/2.0d0/
C
      Call TStamp(1,'Top of NumIn3.')
      Call DecoSC(IOpCl,NSpBlk,NRI,NDimBl,CBfn,NSpBlX,Spinor)
      Call AClear(NBas6D,BasNrm)
      IOpClX = NSpBlX - 1
      Call LClear(MaxAn+1,DidAtm)
      Call DecPrn(IOut,IPrint,IRadAn,0,IPrune,NTheta,NPhi,MaxRad,MaxAng,
     $  ECTolR)
      NeedCR = InPar.and..not.ParDig
      MinPt = MaxAng
      DoWrt = .False.
C     ###%%% Added a possibility to turn off the formation of Density
      DoRho = IUseDn.ne.0
C     ###%%%
      UseP = Mod(IUseDn,2).eq.1
      IPrinL = Max(IPrint-3,0)
      PrtVal = IPrinL.ge.7
      UseDWt = IDerNu.gt.0
      UseTI = UseDWt
      DoTIDE = IDerNu.eq.1
      NoAtTx = .False.
      NMat = 1
      ISuper = Mod(IUseDn,8)/2
      DoHir = LBit(IUseDn,2)
      DoPhi = Mod(IUseDn,16).ge.8
      MaxNG = IAbs(MaxNGI)
      MaxTyX = MaxTyp
      DoIntr = NDBF.eq.0.and.DoSuper.and.MaxNGI.lt.0
      DoAtSD = DoSuper.and..not.DoIntr.and.NDBF.eq.0
      IRType = IRTypI
      If(IRType.eq.0) IRType = 2
      If(ISuper.ge.2) then
        NRhoS = NAtomA
      else if(ISuper.eq.1) then
        NRhoS = 1
      else
        NRhoS = 0
        endIf
      SSWeig = IRanWt.eq.0.or.IRanWt.eq.-1.or.IRanWt.le.-4
      LScr = Max(LScrBf,LScrMO,LScrDn)
      PrnSph = .True.
      IDRhTy = Min(IDerRh,1)
      NonLoc = IDerRh.gt.0
      DoTau = IfTau.eq.1
      DoLap = IfLap.eq.1
      Call NDerX1(IOpCl,IDerNu,DoNuFx,NonLoc,IfTau,IfLap,MxDnXC,MxTyXC,
     $  IVDen,IVTyp,IVC2T,IVT2C,IVLoc,NVarVT,NVarVC,IUDen,IUTyp,
     $  IUC2T,IUT2C,IULoc,NVarUT,NVarUC,MaxDrW,MaxDrB,MaxDrR,IDRhTy)
      MaxDrB = Max(MaxDrB,IDerBf)
      If(UseMO.and.IDerBf.gt.0.and.Mod(IDRhTy,2).eq.0)
     $  IDRhTy = IDRhTy + 1
      IOptMO = IDerBf + 1
      MMax = ((MaxTyX+1)*(MaxTyX+2)*(MaxTyX+3))/6
      jFacAt = 1
      jRL    = jFacAt + NAtoms
      jScale = jRL + MaxTyX + 1
      iCnvrt = IFRPos(jScale+MMax)
      jStart = ITRPos(iCnvrt+MMax)
      Call TstCor(jStart,LenV,'NumIn3-1')
      Call ASet(NAtoms,One,V(jFacAt))
      Call GenScl(.True.,MaxTyX,V(jRL),V(jScale))
      Call LdCnvr(1,.True.,MaxTyX,I(iCnvrt))
      MaxDrW = IDerNu
      If(UseMO) then
        IUseMO = 3
      else
        IUseMO = 1
        endIf
      NMatX = NMat
      NBsAll = NBas6D
      Call Quad1(IOut,IPrint,IOpCl,IPart,DoVarX,DoVarG,Packed,DoSuper,
     $  DoSpar,StoBW2,DoNuFx,DoCPKS,DoGmPe,DoGmPx,PeGmx,DoHyp,TDKSLg,
     $  DoPolD,DoGxPT,DoGmx,PKSDen,NDPol,LScr,NRhoS,0,DoRho,DoNDXC,
     $  DoVN2,DynPar,LinDyn,CountN,NxtVal,NPrtUL,IRadAn,IRanGd,IRanWt,
     $  MaxDrB,MaxDrR,0,IDRhTy,IDerNu,MaxDrW,NDisMt,MinPt,MaxPtI,MaxMBI,
     $  NVarVT,IVLoc,NVarUT,NShell,ShellT,ShellC,MaxTyp,NBas6D,NBsAll,
     $  NShelB,ShelTB,ShelCB,MaxTyB,NDBF,MaxNGI,NAtomA,NAtAll,NClRep,
     $  NMtPBC,IChHar,MulHar,IAn,AtmChg,Thresh,IUseMO,UseP,NMat,NMatX,
     $  NMO,JJ,PA,XX,XX,PB,XX,XX,jStart,ScrnBf,ScrnGd,DoMicB,SSWeig,
     $  DoIntr,DubWgt,MaxPt,MaxMB,KSigSt,NMatP,NMatDP,NDPMax,ND1E,
     $  NShlDB,NDBFU,MaxFPS,NFPS,NFPSI,MaxFTr,NFPS2,IPrune,IPrun2,
     $  IPrun3,MaxRad,MaxAng,MxDrW1,NColWt,MinBat,NSgDSh,ICBat,ICAtom,
     $  ICMBat,MyTask,NTotA,NTotC,NTotB,NTotMB,MaxSA1,MaxSA2,NSgAtW,
     $  MaxL,RCrit,RMax2,ScalSA,ScalSB,ScalPh,CutAt1,CutAt2,ElSumL,
     $  ThrQad,CM,PMax,ECTol,LenIJV,IndIJV,I,V,LenV)
      jBStrt = jBScr1
      Call AlBWSc(0,NAtomA,jBStrt,jBScr1,jBScr2,jBScr3,jBScr4,jBScr5,
     $  jBScr6,jBScr7,jBScr8,jBScr9,jBScrA,jBScrB,jBScrC,jBScrD,jBScrE,
     $  JTest)
C
C     Do the integration.
C
      NTotA = 0
      ElSumD = Zero
      ElSumM = Zero
      ElSumA = Zero
      NOpAtG = 1
      jV = jChi
      MDV1 = LenV - jV + 1
      Do 200 IAtom = 1, NAtomA
        Call Quad2(IOut,IPrint,IPart,IAtom,AtBlock,DoVarX,UseTI,UseDWt,
     $    DoTIDE,NoAtTx,DoWrt,OneRad,PrnSph,IRadAn,IPrune,IPrun2,IPrun3,
     $    MaxRad,MaxAng,IRType,MaxPt,NMatDP,NOpAtG,ND1E,NAtoms,NAtomA,
     $    IAn,C,CM,PJig,TJig,CutAt1,CutAt2,MapAtm,DoAtom,JJ,JJ,MaxAn,
     $    DidAtm,AtmSiz,V(jFacAt),HvRMax,RQMax,IAtomR,ICAtom,IAtomB,
     $    NAngLa,NAngSy,IAtmTI,IAtmTx,NTheta,NPhi,NBatch,JJ,I(iNRdPt),
     $    I(iNAnPt),V(jRadPt),V(jRadWt),XX,XX,SkipAt)
        If(SkipAt) goto 200
        IRadEn = 0
C       Loop over batches of spheres until done...
        Do 130 IBatch = 1, NBatch
          Call Quad3(IOut,IPrint,SSWeig,Packed,AtBlock,ScrnBf,ScrnGd,
     $      DubWgt,PrtVal,UseDWt,DoMicB,DynPar,NAtomA,NAtomA,NAtCel,
     $      NShell,NShlDB,NBas6D,MaxPt,MxSigA,NClRep,MxDrW1,IDerNu,
     $      MinBat,MaxMB,NColWt,IPart,NPart,MyTask,NTheta,NPhi,IAtom,
     $      IAtomR,ICAtom,IBatch,NBatch,I(iNRdPt+IBatch-1),
     $      I(iNAnPt+IBatch-1),IRadEn,ICBat,NMicBt,NSigAt,NSigSh,NSigBf,
     $      NSgDSh,NSgDBf,MaxSA1,NSgAtW,NTotA,NTotB,NGridT,NAngLa,
     $      NAngSy,NTotMB,SizInc,IRanWt,RMax2,RInc,ThrQad,RCrit,BegShl,
     $      BegShB,JAn,I(iNBsSh),I(ISigAt),I(iSigSh),I(iSgCel),
     $      I(iSigBf),I(iSgBfC),I(iSgBfS),I(iSgBfA),I(iMuSt2),I(iSgSAt),
     $      i(iInvSA),I(iNBsSD),I(iSgDSh),I(iSgDCl),I(iSgDBf),I(iSgDBC),
     $      I(iSgDBS),I(iSgDBA),I(iDMuS2),I(iSgDSA),I(iAtSgW),I(iNSame),
     $      I(iKeep),DoInAt,I(iNIdx),IAn,C,AtmSiz,ShlSiz,ShlSzB,
     $      V(jRadPt),V(jRadWt),V(jAngPt),V(jAngWt),V(jGrid),V(jWeigt),
     $      V(jD1Wgt),V(jD2Wgt),NOpAtG,NEqAtm,RotOp,NDisMt,DisMat,JJ,
     $      I(iNNAtm),V(jBScr1),V(jBScr2),V(jBScr3),V(jBScr4),V(jBScr5),
     $      V(jBScr6),V(jBScr7),V(jBScr8),V(jBScr9),V(jBScrA),V(jBScrB),
     $      V(jBScrC),V(jBScrD),V(jBScrE),V(jRg),V(jChiMx),V(jV),MDV1)
          NSgDS2 = 0
          NDone = 0
          Do 120 MBatch = 1, NMicBt
            ICMBat = ICMBat + 1
            NGrid = I(iNIdx+MBatch-1)
            If(NGrid.gt.MaxMB) Call GauErr('NGrid>MaxMB in NumIn3.')
            jkG1 = jGrid + 3*NDone
            kWeigt = jWeigt + NDone*NColWt
            jD2ChD = jD2Chi + 6*NBsAll*NGrid
            jD3ChD = jD3Chi + 10*NBsAll*NGrid
            jRMO = jChi + NBsAll*NGrid
            jRDMO = jD1Chi + 3*NBsAll*NGrid
C           Set up information for this microbatch.
            Call Quad4(IOut,ScrnBf,Packed,DoAtSD,PBC,IDerNu,IAtom,
     $        IBatch,MBatch,NGrid,NAtAll,NBsAll,NBsTot,MxSigA,NBas6D,
     $        NClRep,NDBFU,NAtomA,C,NAtCel,NSigSh,NSigAt,IAtmTI,IAtSTI,
     $        NShlDB,NSgDSh,NShell,I(iNBsSh),I(iSigSh),ShlSiz,I(iSgCel),
     $        I(iSigBf),I(iSgBfC),I(iSgBfS),I(iMuSt2),I(iSgBfA),
     $        I(iSigAt),JAn,I(iNBsSD),I(iSgDSh),JAnB,ShlSzB,I(iSgDCl),
     $        I(iSgDBF),I(iSgDBC),I(iSgDBS),I(iDMuS2),I(iSgDBA),Thresh,
     $        V(jkG1),V(kWeigt),KSigS2,NSigBf,KSigBf,NSgAt2,NSigS2,
     $        NSgDBf,KSgDBf,NSgDS2,MaxSA2,NBsCur,NDBCur,BatThr,
     $        I(iShReg),I(iSgAt2),I(iSigS2),I(iSgCl2),I(iSgSA2),
     $        V(jInSgA),AtmSzD,I(iSgDS2),I(iSgDC2),I(iSgDA2))
C           Evaluate basis functions.
            Call Quad5(IOut,IPrint,IOpClX,NDimBl,Packed,DoVN2,IDerNu,
     $        MaxDrB,MaxDrR,IDRhTy,IUseMO,MaxL,IAtom,IBatch,MBatch,
     $        NGrid,NAtAll,NBsAll,NBsTot,NDBF,IAtmTx,IAtSTI,UseMOs,
     $        UseSig,CntCel,NShell,EXX,C1,C2,C3,C4,X,Y,Z,ShellA,ShellN,
     $        ShellT,ShellC,ShlADF,NShelB,EXXB,C1B,C2B,C3B,C4B,XB,YB,ZB,
     $        ShelAB,ShelNB,ShelTB,ShelCB,ShlADB,I(iNBsSh),I(iSigBf),
     $        I(iSgBfC),I(iSgBfS),I(iMuSt2),I(iNBsSD),I(iSgDBF),
     $        I(iSgDBC),I(iSgDBS),I(iDMuS2),I(iSgDBA),IVLoc,I(iCnvrt),
     $        V(jScale),I(iMDDrv),I(iLocDr),I(nMDDrv),I(iLocTg),
     $        I(iNOrb),I(iOrb),I(iLook),ScalSA,ScalSB,V(jkG1),V(kWeigt),
     $        NSigBf,KSigBf,NSgAt2,NSigS2,NSgDBf,KSgDBf,NSgDS2,NBsCur,
     $        NDBCur,NMO,I(iSgAt2),I(iSigS2),I(iSgCl2),I(iSgDS2),
     $        I(iSgDC2),V(jR),V(jExpT),V(jZeroM),NColW,V(jW),V(jAlpha),
     $        V(jDCoef),V(jChi),V(jD1Chi),V(jD2Chi),V(jD3Chi),V(jD4Chi),
     $        V(jD2ChD),V(jD3ChD),BasNrm,V(jPhiS),V(jD1PhN),V(jVarVS),
     $        V(jVaVSN),V(jVVSN2))
            NTotC = NTotC + NGrid
            If(NSgDS2.eq.0.and.(KSigS2.eq.0.or.NSigBf.eq.0)) goto 120
C           Evaluate density.
            Call Quad6(IOut,IPrint,IOpClX,NSpBlX,NDimBl,DoRho,DoTau,
     $        DoLap,.False.,.False.,DoSpar,DoZCmp,DoVN2,DoSuper,DoHir,
     $        DoIntr,DoPhi,DoV,MaxDrR,IDRhTy,UseMOs,UseSig,IOptMO,NGrid,
     $        NAtAll,NBsAll,NBas6D,NMtPBC,NDBF,NDPMax,0,NAtomA,IAn,
     $        IAtTyp,AtmChg,C,IAtmTx,IAtSTI,IAtomB,IntCel,BckCel,BckDim,
     $        MaxNG,NAtG,AtEC,NShell,KSigSt,RLam,I(iNBsSh),I(iSigBf),
     $        I(iSgBfC),I(iSgBfS),I(iMuSt2),I(iSgBfA),IVLoc,IULoc,
     $        Thresh,ScalSA,ScalSB,ScalPh,PMax,V(jkG1),V(kWeigt),NSigBf,
     $        KSigBf,NSgAt2,NSigS2,NMO,0,CMO,XX,JJ,PA,PB,BatThr,RJR,
     $        I(iShReg),I(iSgAt2),I(iSigS2),I(iSgSA2),V(jChi),V(jD1Chi),
     $        V(jD2Chi),V(jD3Chi),V(jChiMx),I(iIMap),I(iMap),I(iNuRgn),
     $        I(iNuOff),V(jTRA),V(jTRB),V(jTDA),V(jTDB),MaxFPS,NFPS,
     $        NFPSI,I(iIScrF),V(jScrFP),LTBB,V(jTBB),V(jCPA),V(jCQA),
     $        V(jP0Max),V(jPhiS),V(jD1PhN),V(jVarV),V(jVarVS),V(jVarVN),
     $        V(jVaVSN),V(jVVN2),XX,XX)
            If(UseMO) then
              ElSumM = ElSumM + SumMOE(UseOcc,NGrid,NMO,OccMO,V(kWeigt),
     $          V(jRMO))
              endIf
            If(IUseDn.ne.0) then
              Call IndXCV(IVLoc,(jVarV-1),NGrid,IVRA,IVRB,IVD1A,IVD1B,
     $          IVTA,IVTB,IVLA,IVLB,IVRX,IVRY,IVD1X,IVD1Y,IVTX,IVTY,
     $          IVLX,IVLY,LVarV)
              P1x = SProd(NGrid,V(kWeigt),V(IVRA))
              If(IOpCl.eq.0) then
                ElSumD = ElSumD + P1x + P1x
              else
                ElSumD = ElSumD + P1x + SProd(NGrid,V(kWeigt),V(IVRB))
                endIf
              Call IndXCV(IVLoc,(jVarVS-1),NGrid,IVSRA,IVSRB,IVSD1A,
     $          IVSD1B,IVSTA,IVSTB,IVSLA,IVSLB,IVSRX,IVSRY,IVSD1X,
     $          IVSD1Y,IVSTX,IVSTY,IVSLX,IVSLY,LVarVS)
              If(ISuper.ne.0) then
                If(ISuper.eq.1) then
                  LimJAt = 0
                else
                  LimJAt = NAtomA - 1
                  endIf
                Do 60 JAt = 0, LimJAt
                  ElSumA = ElSumA +
     $              Two*SProd(NGrid,V(kWeigt),V(IVSRA+JAt*NGrid))
   60             Continue
                endIf
              endIf
            If(NeedCR) then
C$OMP Critical (NumIn3_1)
              If(UseBfn) Call BfnFnc(IOut,IPrint,IPart,NPart,IAtom,
     $          NAtomA,NBas6D,NGrid,NSigBf,I(iSigBf),IntBfn,WPBfn,IAn,
     $          IAtTyp,AtmChg,C,V(jkG1),V(kWeigt),V(jChi),V(jD1Chi),
     $          V(jD2Chi),V(jD3Chi),LScrBf,V(jScr))
              If(UseMO) Call MOFunc(IOut,IPrint,IPart,NPart,IAtom,
     $          NAtomA,NBas6D,UseOcc,NMO,NGrid,NSigBf,I(iSigBf),OccMO,
     $          IntMO,WPMO,IAn,IAtTyp,AtmChg,C,V(jkG1),V(kWeigt),
     $          V(jChi),V(jD1Chi),V(jD2Chi),V(jD3Chi),V(jRMO),V(jRDMO),
     $          LScrMO,V(jScr))
              If(IUseDn.ne.0) Call DenFnc(IOut,IPrint,IPart,NPart,IOpCl,
     $          IAtom,NAtomA,NBas6D,NGrid,NSigBf,I(iSigBf),IntDen,WPDen,
     $          IAn,IAtTyp,AtmChg,C,V(jkG1),V(kWeigt),V(jChi),V(jD1Chi),
     $          V(jD2Chi),V(jD3Chi),V(IVRA),V(IVRB),V(IVD1A),V(IVD1B),
     $          V(IVSRA),V(IVSD1A),V(jPhiS),V(IVRX),V(IVRY),V(IVD1X),
     $          V(IVD1Y),LScrDn,V(jScr))
C$OMP End Critical (NumIn3_1)
            else
              If(UseBfn) Call BfnFnc(IOut,IPrint,IPart,NPart,IAtom,
     $          NAtomA,NBas6D,NGrid,NSigBf,I(iSigBf),IntBfn,WPBfn,IAn,
     $          IAtTyp,AtmChg,C,V(jkG1),V(kWeigt),V(jChi),V(jD1Chi),
     $          V(jD2Chi),V(jD3Chi),LScrBf,V(jScr))
              If(UseMO) Call MOFunc(IOut,IPrint,IPart,NPart,IAtom,
     $          NAtomA,NBas6D,UseOcc,NMO,NGrid,NSigBf,I(iSigBf),OccMO,
     $          IntMO,WPMO,IAn,IAtTyp,AtmChg,C,V(jkG1),V(kWeigt),
     $          V(jChi),V(jD1Chi),V(jD2Chi),V(jD3Chi),V(jRMO),V(jRDMO),
     $          LScrMO,V(jScr))
              If(IUseDn.ne.0) Call DenFnc(IOut,IPrint,IPart,NPart,IOpCl,
     $          IAtom,NAtomA,NBas6D,NGrid,NSigBf,I(iSigBf),IntDen,WPDen,
     $          IAn,IAtTyp,AtmChg,C,V(jkG1),V(kWeigt),V(jChi),V(jD1Chi),
     $          V(jD2Chi),V(jD3Chi),V(IVRA),V(IVRB),V(IVD1A),V(IVD1B),
     $          V(IVSRA),V(IVSD1A),V(jPhiS),V(IVRX),V(IVRY),V(IVD1X),
     $          V(IVD1Y),LScrDn,V(jScr))
              endIf
  120       NDone = NDone + NGrid
          If(DynPar.and.NMicBt.ge.0)
     $      MyTask = GLinCo(5,LinDyn,CountN,NxtVal,NPrtUL,1)
  130     Continue
  200   Continue
      Return
      End
