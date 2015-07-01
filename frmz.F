      subroutine frmz(iout, iprint, ngrid, irwspo, irwe,
     &     irwz, v, lenv) 

      implicit none 

C     Performs the last step of the integrals decomposition:
C     Z = sum_{pp,qq} S^(-1)(p,pp) E(pp,qq) S(-1)(qq,q)
C
C     This routine is here only to allow for block BLAS or
C     some matrix multiplication optimizations in future
C
C     Arguments:
C
C     [in]  ngrid  - size of the square matrices S and E
C     [in]  irwspo - RWF with inverse overlap
C     [in]  irwe   - RWF with projected integrals
C     [out] irwz   - RWF with the resulting 2e electron operator Z

C     Arguments:

      integer    iout, iprint, lenv
      real*8     v(*)

      integer    ngrid, irwspo, irwe, irwz

C     Other

      integer    jend, jscr1, jscr2, jscr3 

      jend = 1
      jscr1 = jend 
      jscr2 = jscr1 + ngrid*ngrid
      jscr3 = jscr2 + ngrid*ngrid
      jend  = jscr3 + ngrid*ngrid

      call tstcor (jend, lenv, 'frmz')
      call conddf (irwz, ngrid*ngrid)
      
      call fileio (2, -irwspo, ngrid*ngrid, v(jscr1), 0)
      call fileio (2, -irwe, ngrid*ngrid, v(jscr2), 0)

      call dgemm ('n', 'n', ngrid, ngrid, ngrid, 1.0d0,  
     & v(jscr1), ngrid, v(jscr2), ngrid, 0.0d0, 
     &     v(jscr3), ngrid)

      call dgemm ('n', 'n', ngrid, ngrid, ngrid, 1.0d0,  
     & v(jscr3), ngrid, v(jscr1), ngrid, 0.0d0, 
     &     v(jscr2), ngrid)

c$$$      call dsymm ('l', 'u', ngrid, ngrid, 1.0d0,
c$$$     &     v(jscr1), ngrid, v(jscr2), ngrid, 0.0d0, 
c$$$     &     v(jscr3), ngrid)
c$$$
c$$$      call dsymm ('l', 'u', ngrid, ngrid, 1.0d0,  
c$$$     &     v(jscr1), ngrid, v(jscr3), ngrid, 0.0d0, 
c$$$     &     v(jscr2), ngrid)
c$$$
      call fileio (1, -irwz, ngrid*ngrid, v(jscr2), 0)

      end
      
