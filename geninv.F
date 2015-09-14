*Deck GenInv
      Subroutine GenInv(IOut,IPrint,Label,NSkip,Tol,N,A,Scr1,Eig,DoAIRt,
     $  AIRt,NRank,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Form the generalized inverse of A.  Eigenvalues smaller than Tol
C     are considered to be zero, and the first NSkip eigenvectors are
C     ignored regardless.  A is lower triangular and is replaced by
C     its generalized inverse.  If DoAIRt, then AIRt is loaded with
C     A**(-1/2).
C
      Logical DoAIRt
      Parameter (MinReg=20,LenScr=(MinReg*10))
      Character*(*) Label
      Dimension A(*), Scr1(N,N), Eig(N), AIRt(*), XX(1), V(MDV),
     $  ScrLoc(LenScr), JJ(1)
      Save XX, JJ
      Data XX/0.0d0/, JJ/0/
 1000 Format(' Eigenvalues and Eigenvectors of ',A,':')
 1010 Format(' Eigenvalues of ',A,':',/,10(1PD10.2))
 1020 Format(' SDV returned info=',I12,' to GenInv.')
C
      If(NSkip.lt.0) Call GauErr('NSkip<0 in GenInv.')
      NTT = (N*(N+1))/2
      Call Square(A,Scr1,N,N,0)
      If(N.ge.MinReg) then
        Call GGESVD('O','N',N,N,Scr1,N,Eig,XX,N,XX,N,A,NTT,Info)
      else
        Call GGESVD('O','N',N,N,Scr1,N,Eig,XX,N,XX,N,ScrLoc,LenScr,
     $    Info)
        endIf
      If(Info.ne.0) then
        Write(IOut,1020) Info
        Call GauErr('SVD failed for GenInv.')
        endIf
      If(IPrint.ge.0) Call TStamp(1,'After SVD in GenInv.')
      Call FixPha(1,N,N,0,Scr1)
      If(IPrint.ge.3) then
        Write(IOut,1000) Label
        Call MatPt1(IOut,Scr1,N,N,N,N,0,0,JJ,JJ,0,Eig,1)
      else if(IPrint.ge.2) then
        Write(IOut,1010) Label, (Eig(I),I=1,N)
        endIf
      Call GenIn1(.True.,.False.,DoAIRt,N,N,NSkip,NRank,Tol,Eig,Scr1,A,
     $  V,AIRt,V,MDV)
      Return
      End
