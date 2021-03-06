      subroutine cpdexmem (iout, iprint, initf, nbasis, ngrid,
     & x, z, nmat, nmats, nmatt, p, f, v, lenv)

      implicit none
      
C-----------------------------------------------------------------------------+
C     This routine does the exchange part of the Fock matrices using CPD
C     decomposed integrals. It is intended to be called from cpdfkmem,
C     and works with one spin case per run.
C     The integrals are Mulliken order, stored as:
C
C     (mn|rs) = X_{I,m} X_{I,n} Z_{I,J} X_{J,r} X_{J,s}
C     (repeated indices are summed)
C      
C     Densities are supplied in lower triangular 
C     format in vectors p, fock matrices are returned lower triangular in 
C     f.
C     
C     The formulas for contraction can be found in ctrcou and ctrexc
C      
C     Arguments:
C     [in]  initf   - if we need to initialize or update the Fock matrices of
C           the triplet and antisymmetric parts. Singlet part is always updated
C     [in]  nbasis  - basis size
C     [in]  ngrid   - quadrature size
C     [in]  x       - nbasis*ngrid quadrature
C     [in]  z       - ngrid*ngrid Z matrix
C     [in]  nmat    - total number of lower triangular matrices in P
C     [in]  nmats   - number of singlets stored in 1..nmats first lower triangular matrices in P
C                     (are contracted with both Coulomb and exchange)
C     [in]  nmatt   - number of triplets stored in nmats+1..nmatt positions in P
C                     (are contracted with exchange only)
C     [in]  p       - vector of nmat lower triangular matrices density matrices
C     [out] f       - vector of nmat lower triangular Fock matrices
C     
C-----------------------------------------------------------------------------+
C     arguments

      logical     initf
      integer     iopcl
      integer     iout, iprint
      integer     nbasis, ngrid
      integer     nmat, nmats, nmatt
      real*8      x(*), z(*)
      real*8      p(*), f(*)
      real*8      v(*)
      integer     lenv

C     pointers & sizes
      
      integer     lenx, lenz
      integer     lenscrex
      
C     other variables

      integer     ntt, nbsq, ngsq
      integer     ii
      integer     ioff, nmata
      
C     define sizes

      ntt  = nbasis*(nbasis+1)/2
      nbsq = nbasis*nbasis
      ngsq = ngrid*ngrid

      lenx = nbasis*ngrid
      lenz = ngsq
      
C     Prepare for parallel run (NYI)
      
      lenscrex = lenx + ngsq
      
      call tstcor (lenscrex, lenv, 'cpdexmem')

C     add exchange to symmetric singlets
      do 100 ii = 0, nmats - 1
         call ctrexc (iout, iprint, .false., 0, nbasis, ngrid, x,
     &        z, p(1+ii*ntt), f(1+ii*ntt), v) 
 100  continue
      
C     do exchange for symmetric triplets
      
      ioff = ntt*nmats + 1
      do 200 ii = 0, nmatt - 1
         call ctrexc (iout, iprint, initf, 0, nbasis, ngrid, x,
     &        z, p(ioff+ii*ntt), f(ioff+ii*ntt), v)
 200  continue

C     do exchange for antisymmetric part

      ioff  = ntt*(nmats + nmatt) + 1
      nmata = nmat - nmats - nmatt
      do 300 ii = 0, nmata - 1
         call ctrexc (iout, iprint, initf, 1, nbasis, ngrid, x,
     &        z, p(ioff+ii*ntt), f(ioff+ii*ntt), v)
 300  continue
      
      return
      
      end
