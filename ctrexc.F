      subroutine ctrexc (iout, iprint, initf, isym, nbasis, ngrid, x,
     &     z, p, f, v)

C-----------------------------------------------------------------------------------------+
C     Helper routine for exchange part contraction, intended to run in parallel.
C     Operation is expressed as:
C
C     F_{m,n} = (ms|rn) * P_{s,r}, where (ms|rn) = X_{m,I} X_{s,I} Z_{I,J} X_{r,J} X_{n,J}
C     (repeated indices are summed)
C     ---------------------------
C     L_{s,I} = (P * X)_{s,J}
C     O_{I,J} = X_{s,I} L_{s,J} = (X^{t} L)_{I,J},     where o is an elementwise (Hadamard) product
C     Q_{I,J}   = Z_{I,J} o O_{I,J} (no sum)
C     F_{m,n} = Q_{I,J} X_{m,I} X_{n,J} = X Q X^{t}
C     diag(n) - form a diagonal matrix from a vector n
C      
C     Arguments:
C     [in]  initf   - if we need to update the current contents of F or initialize it
C     [in]  isym    - symmetry of the lower triangulat densities: 0 - symmetric, 1 - antisymmetric
C     [in]  nbasis  - basis size
C     [in]  ngrid   - quadrature size
C     [in]  x       - quadrature
C     [in]  z       - Z matrix
C     [in]  p       - current singlet density to contract, lower triangular 
C     [out] f       - current Fock matrix to output
C
C      
C-----------------------------------------------------------------------------------------+

C     arguments

      logical      initf
      integer      iout, iprint, isym, nbasis, ngrid
      real*8       x(*), z(*)
      real*8       p(*), f(*), v(*)

C     pointers

      integer      jend, jscr, jmvec, jnvec, jl

C     other
      
      integer      ii, nbsq, ngsq, lenx
      
      nbsq    = nbasis*nbasis
      ntt     = nbasis*(nbasis+1)/2
      ngsq    = ngrid*ngrid
      lenx    = nbasis*ngrid
      scrsz   = nbsq + lenx
      
      jscr   = 1
      jl       = jscr + nbsq
      jo       = jl     + lenx
      jend     = jo     + ngsq
      call tstcor (jend - 1, scrsz, 'ctrexc')

      if (isym .eq. 0) then
         call square (p,v(jscr),nbasis,nbasis,3)
      elseif (isym .eq. 1) then
         call square (p,v(jscr),nbasis,nbasis,-3)
      else
         call gauerr ('ctrexc: unknown symmetry switch passed')
      endif

C     Make the exchange negative

      call ascale (nbsq, -1.0d0, v(jscr), v(jscr))
      
C     form L

      if (isym .eq. 0) then
         call dsymm ('l', 'u', nbasis, ngrid, 1.0d0,
     &        v(jscr), nbasis, x, nbasis, 0.0d0, 
     &        v(jl), nbasis)
      else
         call dgemm ('n', 'n', nbasis, ngrid, nbasis, 1.0d0, v(jscr),
     &        nbasis, x, nbasis, 0.0d0, v(jl), nbasis)
      endif
      
C     form O

      call dgemm ('t', 'n', ngrid, ngrid, nbasis, 1.0d0, x,
     &     nbasis, v(jl), nbasis, 0.0d0, v(jo), ngrid)

C     form Q
      
      do 100 ii = 0, ngsq - 1
         v(jo + ii) = v(jo + ii)*z(1+ii)
 100  continue
      
C     form F

      call dgemm ('n', 'n', nbasis, ngrid, ngrid, 1.0d0, x,
     &     nbasis, v(jo), ngrid, 0.0d0, v(jl), nbasis)
      
      call dgemm ('n', 't', nbasis, nbasis, ngrid, 1.0d0, v(jl),
     &     nbasis, x, nbasis, 0.0d0, v(jscr), nbasis)
      
C     return the result

      if (initf) then
         call linear (v(jscr), f, nbasis, nbasis)
      else
         call linear (v(jscr), v(jscr), nbasis, nbasis)
         call aadd (ntt, v(jscr), f, f)
      endif
      
      end
