*Deck GenIn1
      Subroutine GenIn1(DoAI,DoARt,DoAIRt,N,NV,NSkip,NRank,Tol,Eig,EVec,
     $  AI,ARt,AIRt,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Given the eigenvectors from SVD, form the generalized inverse of a
C     symmetric matrix.  NRank is returned.
C
      Logical DoAI, DoARt, DoAIRt, UseGMM
      Dimension Eig(*), EVec(N,*), AI(*), ARt(*), AIRt(*), V(*)
      Save Zero, One
      Data Zero/0.0d0/, One/1.0d0/
C
      NRank = 0
      NTT = (N*(N+1))/2
      NSq = N*N
      UseGMM = MDV.ge.NSq
      If(UseGMM) then
        Do 10 I = 1, (NV-NSkip)
          If(Eig(I).gt.Tol*Eig(1)) then
            NRank = NRank + 1
            If(NRank.ne.I) Call GauErr('GenIn1 is confused-10.')
            endIf
   10     Continue
        If(DoAIRt) then
          Call ScalDg(N,N,NRank,EVec,EVec,Eig,.False.,-1,IErr)
          Call XGEMM(1,'N','T',N,N,NRank,One,EVec,N,EVec,N,Zero,V,N)
          Call Linear(V,AIRt,N,N)
          Call ScalDg(N,N,NRank,EVec,EVec,Eig,.False.,-1,IErr)
        else
          Call ScalDg(N,N,NRank,EVec,EVec,Eig,.False.,-2,IErr)
          endIf
        If(DoAI) then
          Call XGEMM(1,'N','T',N,N,NRank,One,EVec,N,EVec,N,Zero,V,N)
          Call Linear(V,AI,N,N)
          endIf
        If(DoARt) then
          Call ScalDg(N,N,NRank,EVec,EVec,Eig,.False.,3,IErr)
          Call XGEMM(1,'N','T',N,N,NRank,One,EVec,N,EVec,N,Zero,V,N)
          Call Linear(V,ARt,N,N)
          Call ScalDg(N,N,NRank,EVec,EVec,Eig,.False.,-1,IErr)
        else
          Call ScalDg(N,N,NRank,EVec,EVec,Eig,.False.,2,IErr)
          endIf
      else
        Do 20 I = 1, (NV-NSkip)
          If(Abs(Eig(I)).gt.Tol*Eig(1)) then
            NRank = NRank + 1
            If(NRank.ne.I) Call GauErr('GenIn1 is confused-20.')
            Eig(I) = One / Eig(I)
            endIf
   20     Continue
        If(DoAI) then
          Call AClear(NTT,AI)
          Do 30 I = 1, N
            II = (I*(I-1))/2
            Do 30 K = 1, NRank
              T = Eig(K)*EVec(I,K)
              Do 30 J = 1, I
   30           AI(II+J) = AI(II+J) + EVec(J,K)*T
          endIf
        If(DoAIRt) then
          Call AClear(NTT,AIRt)
          Do 40 K = 1, NRank
            S = Sqrt(Eig(K))
            Do 40 I = 1, N
              II = (I*(I-1))/2
              T = S*EVec(I,K)
              Do 40 J = 1, I
   40           AIRt(II+J) = AIRt(II+J) + EVec(J,K)*T
          endIf
        Do 50 I = 1, NRank
   50     Eig(I) = One / Eig(I)
        If(DoARt) then
          Call AClear(NTT,ARt)
          Do 60 K = 1, NRank
            S = Sqrt(Eig(K))
            Do 60 I = 1, N
              II = (I*(I-1))/2
              T = S*EVec(I,K)
              Do 60 J = 1, I
   60           ARt(II+J) = ARt(II+J) + EVec(J,K)*T
          endIf
        endIf
      Return
      End
