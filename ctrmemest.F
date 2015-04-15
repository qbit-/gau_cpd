      subroutine ctrmemest(iout, iprint, nbasis, ngrid, incorepos, 
     &     numkl, nump, numq, lenwbuf, lenfbuf, lenxbufp, lenxbufq,
     &     lenscr, lenebuf, v, lenv) 

      implicit none 

C     estimates different sizes for the formation of the E matrix 
C
C     Arguments:
C     
C     [in] nbasis - number of basis functions
C     [in] ngrid  - number of grid points in the quadrature
C     [out] incorepos - if fully in-core operation is possible
C     [out] numkl - step size over the KL index of integrals. Here we allocate to iterate over whole counts of L only
C     [out] nump  - step size over the right quadrature grid index
C     [out] numq  - step size over the left quadrature grid index
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
      integer    numkl, nump, numq

      integer    lenwbuf, lenfbuf, lenxbufp, lenxbufq, lenebuf, lenscr
      logical    incorepos

C     Other

      integer    numk, numl
      integer    nbsq, npmax, res, pqmax

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

      if (res .le. 0) call gauerr ('Insufficient memory in ctrmemest:1')

C     Just set fixed values in a dumb way for now

      nump = ngrid
      numq = ngrid
      pqmax = max(nump, numq)

      res = lenv - nump*nbasis*(npmax + 1) - numq*nbasis*(npmax + 1) 
     &     - nump*numq

      if (res .lt. 0)
     &     call gauerr('Insufficient memory in ctrmemest:2')

      numkl = res / (nbsq + nump)
      if (numkl .gt. nbsq) numkl = nbsq
      if (numkl .lt. 1)
     &     call gauerr('Insufficient memory in ctrmemest:3')

      numl  = numkl / nbasis
      if (numl .ne. 0) then
         numk  = nbasis
         numkl = numl*nbasis
      else
         numl  = 1
         numkl = numk
      endif
      
C     Now assign all sizes

      lenwbuf = numkl*nbsq
      lenfbuf = numkl*nump
      lenxbufp  = nump*nbasis
      lenxbufq  = numq*nbasis
      lenscr    = pqmax*nbasis*npmax
      lenebuf   = nump*numq

      end
