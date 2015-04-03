      subroutine ctrmemest(iout, iprint, nbasis, ngrid, incorepos, 
     &     numk, numl, nump, numq, lenwbuf, lenxbufp, lenxbufq,
     &     lenscr, v, lenv) 

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
C     [out] nump - step size over the right quadrature projection
C     [out] numq - step size over the left quadrature projection
C     [out] lenwbuf - length of the main working buffer 
C     [out] lenxbufp - length of the buffer for right quadratures
C     [out] lenxbufq - length of the buffer for left quadratures
C     [out] lenscr   - scratch space length
C     [in]  lenv     - available memory



C     Arguments

      integer    iout, iprint, lenv
      real*8     v(*)

      integer    nbasis, ngrid
      integer    numk, numl, nump, numq

      integer    lenwbuf, lenxbufp, lenxbufq, lenscr
      logical    incorepos

C     Other

      integer    nbsq, npmax, res
      integer    numkl, pqmax

C     functions 
      
      integer    nproc

      nbsq  = nbasis*nbasis
      npmax = nproc(0) 

C     We make a crude estimation if intergrals can be kept in memory
      if ((lenv - nbsq*nbsq - nbasis*nbsq - nbsq*npmax - 2*nbsq) .gt. 0)
     &  then
         incorepos = .true.
         res = lenv - nbsq*nbsq
      else
         incorepos = .false.
         res = lenv
      endif

C     try to allocate space for as much kl indices as possible, keeping a minimal space for 
C     scratch 

      numkl = (res - npmax*nbsq - 2*nbsq) / nbsq
      if (numkl .gt. nbsq) numkl = nbsq
      if (numkl .lt. 1) 
     & call gauerr('Insufficient memory in ctrmemest:1')
      
      numk = numkl / nbasis
      if (numk .lt. 1) then 
         numl  = mod (numkl, nbasis) 
         numk  = 1
      else
         numl  = nbasis
      endif

      res = res - nbsq*numk*numl

C     Done with main working buffer. Now try to maximize auxiliarry buffers
C     right quadrature array

      nump = (res - nbsq) / ( nbasis * (npmax + 1) )
      if (nump .gt. numk*numl) nump = numk*numl
      if (nump .gt. nbasis)  nump = nbasis
     
      res = res - nump*nbasis*( npmax + 1 )
      if (res .lt. nbsq)
     &     call gauerr('Insufficient memory in ctrmemest:2')

C     left quadrature array
      
      numq = ( res + nump*npmax*nbasis ) / (nbasis*(npmax + 1))
      if (numq .gt. nbsq/numk) numq = nbsq/numk
      
      pqmax = max(nump,numq)
      res = res + nump*nbasis*npmax - pqmax*nbasis*npmax - numq*nbasis
      
      if (res .lt. 0)
     &     call gauerr('Insufficient memory in ctrmemest:3')

C     Now assign all sizes

      lenwbuf  = numk*numl*nbsq
      lenxbufp = nump*nbasis
      lenxbufq = numq*nbasis
      lenscr   = pqmax*nbasis

      end
