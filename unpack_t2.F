      subroutine unpack_t2 (iout, no, nv,
     &     x, v, lenv)

      implicit none


C +------------------------------------------------------------+
C |  unpack_t2   ---  unpacks t2 antisymmetric part            |
C |                   to a rectangular full tensor             |
C |                   T2 is assumed to be stored as            |
C |                   i<j,a<b linear array with i,j < no       |
C |                   a,b < nv and is unpacked to              |
C |                   no*no*nv*nv                              |
C |                   Full array is assumed to be properly     |
C |                   allocated                                |
C +------------------------------------------------------------+
C     Arguments:
C     [in]  no, nv             - dimensions           
C     [in/out]  x              - array to unpack with filled with packed data
C     arguments

      integer    lenv
      integer    iout
      real*8     v(*)

      integer    no, nv 
      real*8     x(*)

C     other variables
      
      integer    jend, ii
      integer    not, nvt, nosq, nvsq
      real*8     xx 
      
      jend = 1
      not  = no*(no-1)/2
      nvt  = nv*(nv-1)/2
      nosq = no*no
      nvsq = nv*nv

   
      call aclear(nosq*nvsq - not*nvt, x)
      call unpcck(0, 0, x, nosq, nvt, not, nvt)

      do 10 ii = 1, nvt
         call square(x(1 + (ii-1)*nosq), x(1 + (ii-1)*nosq),
     &        no, no, -1)
 10   continue

      call tstcor(nosq*nvsq, lenv, 'unpack_t2')
      
      call trcor1(1, 1, nvt, nosq, 3, x, xx, xx, v(jend),
     &     xx, xx, .false., xx, xx)

      call aclear(nosq*nvsq, x(1))
      do 20 ii = 1, nosq
         call square(v(jend + (ii-1)*nvt), x(1 + (ii-1)*nvsq),
     &        nv, nv, -3)
 20   continue
      
      end
