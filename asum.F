      function asum (n, dx, incx)
      Implicit Real*8(A-H,O-Z)
C     
C     Sum elements of a real vector (NOT an absolute value sum)
C
      Dimension DX(IncX,*)
      Save Zero
      Data Zero/0.0d0/
C                                                                                                               
      ASum = Zero
      Do 10 I = 1, N
 10         ASum = ASum + DX(1,I)
      Return
      end
