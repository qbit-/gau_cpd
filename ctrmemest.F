      subroutine ctrmemest(iout, iprint, nbasis, ngrid, incorepos, 
     &     numk, numl, nump, numq, lenwbuf, lenfbuf, lenxbufp, lenxbufq,
     &     lenscr, lenebuf, v, lenv) 

      implicit none 

C     estimates different sizes for the formation of the E matrix 
C
C     Arguments:
C     
C     [in] nbasis - number of basis functions
C     [in] ngrid  - number of grid points in the quadrature
C     [out] incorepos - if fully in-core operation is possible
C     [out] numk - step size over the K index of integrals - slowest dimension
C     [out] numl - step size over the L index of integrals - 2nd slowest dimension
C     [out] nump - step size over the right quadrature grid index
C     [out] numq - step size over the left quadrature grid index
C     [out] lenwbuf - length of the integral buffer 
C     [out] lenfbuf - length of the left projected kl integrals
C     [out] lenxbufp - length of the buffer for right quadratures
C     [out] lenxbufq - length of the buffer for left quadratures
C     [out] lenscr   - scratch space length
C     [out] lenebuf  - space for final E buffer - right and left projected integrals
C     [in]  lenv     - available memory



C     Arguments

      integer    iout, iprint, lenv
      real*8     v(*)

      integer    nbasis, ngrid
      integer    numk, numl, nump, numq

      integer    lenwbuf, lenfbuf, lenxbufp, lenxbufq, lenebuf, lenscr
      logical    incorepos

C     Other

      integer    nbsq, npmax, res
      integer    numkl, pqmax

C     functions 
      
      integer    nproc

      nbsq  = nbasis*nbasis
      npmax = nproc(0) 

C     We make a crude estimation if intergrals can be kept in memory
      if ((lenv - nbsq*nbsq - nbasis*nbsq - 3*nbsq - npmax*nbsq) .gt. 0)
     &  then
         incorepos = .true.
         res = lenv - nbsq*nbsq
      else
         incorepos = .false.
         res = lenv
      endif

C     Just set fixed valus in a dumb way for now

      nump = nbasis
      numq = nbasis
      pqmax = max(nump, numq)

      res = lenv - nump*nbasis*(npmax + 1) - numq*nbasis*(npmax + 1) 
     &     - nump*numq

      if (res .lt. 0)
     &     call gauerr('Insufficient memory in ctrmemest:1')

      numkl = res / (nbsq + nump)
      if (numkl .lt. 1)
     &     call gauerr('Insufficient memory in ctrmemest:2')

      numl  = numkl / nbasis
      if (numl .lt. 1) then 
         numk = mod(numkl, nbasis)
      else
         numk = nbasis
      endif

C     Now assign all sizes

      lenwbuf   = (numk + numl*nbasis)*nbsq
      lenfbuf   = (numk + numl*nbasis)*nump
      lenxbufp  = nump*nbasis
      lenxbufq  = numq*nbasis
      lenscr    = pqmax*nbasis*npmax
      lenebuf   = nump*numq

      end
