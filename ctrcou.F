      subroutine ctrcou (iout, iprint, iupdate, isym, nbasis, ngrid, x,
     & z, p, f, v)

C-----------------------------------------------------------------------------------------+
C     Helper routine for Coulomb part contraction, intended to run in parallel.
C     Operation is expressed as:
C
C     F_{m,n} = (mn|rs) * P_{s,r}, where (mn|rs) = X_{m,I} X_{n,I} Z_{I,J} X_{r,J} X_{s,J}
C     (repeated indices are summed)
C     ---------------------------
C     L_{s,J} = (P * X)_{s,J}
C     m_{J}   = trace_{s} (X o L)_{J},     where o is an elementwise (Hadamard) product
C     
C     n_{I}   = Z_{I,J} M_{J}
C     F_{m,n} = n_{I} X_{m,I} X_{n,I} = X diag(n) X^{t}
C     diag(n) - form a diagonal matrix from a vector n
C      
C     Arguments:
C     [in]  iupdate - if we need to update the current contents of F
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

      logical      iupdate
      integer      iout, iprint, isym, nbasis, ngrid
      real*8       x(*), z(*)
      real*8       p(*), f(*), v(*)

C     pointers

      integer      jend, jscr, jmvec, jnvec, jl

C     other
      
      integer      ii, nbsq, ntt, lenx

C     functions
      
      real*8     ddot, dsdot
      
      nbsq    = nbasis*nbasis
      ntt     = nbasis*(nbasis+1)/2
      lenx    = nbasis*ngrid
      scrsz   = nbsq + lenx + ngrid
      
      jscr   = 1
      jl       = jscr + nbsq
      jnvec    = jl     + lenx
      jend     = jnvec  + ngrid
      call tstcor (jend - 1, scrsz, 'ctrexc')
      
      if (isym .eq. 0) then
         call square (p,v(jscr),nbasis,nbasis,3)
      elseif (isym .eq. 1) then
         call ascale (ntt, -1.0d0, p, v(jscr))
         call square (v(jscr),v(jscr),nbasis,nbasis,-1)
      else
         call gauerr ('ctrcou: unknown symmetry switch passed')
      endif
      
C     form L

      if (isym .eq. 0) then
         call dsymm ('l', 'u', nbasis, ngrid, 1.0d0,
     &        v(jscr), nbasis, x, nbasis, 0.0d0, 
     &        v(jl), nbasis)
      else
         call dgemm ('n', 'n', nbasis, ngrid, nbasis, 1.0d0, v(jscr),
     &        nbasis, x, nbasis, 0.0d0, v(jl), nbasis)
      endif

C     form m

      jmvec    = jl  ! Watch out if it can be corrupted by ddot
      do 100 ii = 0, ngrid - 1
         v(jmvec + ii) = ddot (nbasis,
     & v(jl + ii*nbasis), 1, x(1 + ii*nbasis), 1)
 100  continue

C     form n

      call dgemv ('n', ngrid, ngrid, 1.0d0, z, ngrid, v(jmvec),
     &   1, 0.0d0, v(jnvec), 1)
      
C     form F

      call amove  (lenx, x, v(jl))
      call ascaln (nbasis, ngrid, v(jnvec), v(jl))
    
      call dgemm ('n', 't', nbasis, nbasis, ngrid, 1.0d0, v(jl),
     &     nbasis, x, nbasis, 0.0d0, v(jscr), nbasis)

C     return the result

      if (iupdate) then
         call linear (v(jscr), v(jscr), nbasis, nbasis)
         call aadd (ntt, v(jscr), f, f)
      else
         call linear (v(jscr), f, nbasis, nbasis)
      endif
      end
