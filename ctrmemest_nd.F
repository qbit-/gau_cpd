      subroutine ctrmemest_nd(iout, iprint, nbasis, ngrid, keepints, 
     &     keepx, numkl, nump, numq, lenwbuf, lenfbuf, lenxbufp,
     &     lenxbufq, lenscr, lenebuf, v, lenv) 

      implicit none 

C     estimates different sizes for the formation of the E matrix 
C
C     Arguments:
C     
C     [in] nbasis - number of basis functions
C     [in] ngrid  - number of grid points in the quadrature
C     [out] keepints - if fully in-core operation is possible
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
      logical    keepints, keepx

C     Other

      integer    nbsq, npmax, res, pqmax

C     functions 
      
      integer    nproc

      nbsq  = nbasis*nbasis
      npmax = nproc(0) 

C     An optimal set of values is a quadratic linear problem. Here we just use heuristics to
C     get some set of values

C     X is used in the loop over P and Q
C     Decide what is bigger - an average O(N^2) integral buffer or the grid size

      if (ngrid/nbasis .lt. 2*nbasis) then
         keepx = .true.
         nump = ngrid
      else
C     grid will be loaded by chunks of O(nbasis) size
         keepx = .false.
         nump = 2*nbasis
      endif

C     The second thing to consider is an integral buffer
C     Total memory needed is:
C     nump*nbasis + nbsq*numkl + nump*numkl + max(numkl*numq,nump*nbasis) + numq*nbasis + nump*numq

C     assume numq = 1 then 
C     approximate max(numkl*numq,nump*nbasis) <= numkl + nump*nbasis

      if (keepx) then
C     In this case we don't need the numq*nbasis term
         if ((lenv - 2*nump*nbasis - nump)/
     &        (nbsq + nump + 1) .gt. nbsq) then
            keepints = .true.
            numkl    = nbsq
         else
            keepints = .false.
            numkl    = (lenv - 2*nump*nbasis - nump)/
     &           (nbsq + nump + 1)
         endif
      else
         if ((lenv - 2*nump*nbasis - nbasis - nump)/
     &        (nbsq + nump + 1) .gt. nbsq) then
            keepints = .true.
            numkl    = nbsq
         else
            keepints = .false.
            numkl    = (lenv - 2*nump*nbasis - nbasis - nump)/
     &           (nbsq + nump + 1)
         endif
      endif

      if (numkl .lt. 1)
     &     call gauerr ('Insufficient memory in ctrmemest:1')
      
C     Now compute the best we can have for numq
      
      if (keepx) then
C     In this case we don't need the numq*nbasis term         
         numq = (lenv - ngrid*nbasis - nbsq*numkl - ngrid*numkl
     &        - ngrid*nbasis)/(numkl + ngrid)
      else
         numq = (lenv - nump*nbasis - nbsq*numkl - nump*numkl -
     &           nump*nbasis)/(numkl + nump + nbasis)
      endif
C     Limit numq if it's larger then the grid size
      if (numq .gt. ngrid) numq = ngrid
      
      if (numq .lt. 1) 
     &     call gauerr('Insufficient memory in ctrmemest:2')
      
C     Now check we are not out of boundaries

      if (keepx) then
         res = lenv - nump*nbasis - nbsq*numkl - nump*numkl -
     &        max(numkl*numq,nump*nbasis) - nump*numq
      else
         res = lenv - nump*nbasis - nbsq*numkl - nump*numkl -
     &        max(numkl*numq,nump*nbasis) - numq*nbasis - nump*numq
      endif
      
      if (res .lt. 0)
     &     call gauerr('Insufficient memory in ctrmemest:3')

C     Now assign all sizes

      if (keepints) then
         lenwbuf = 0
      else
         lenwbuf = numkl*nbsq
      endif
      lenfbuf = numkl*nump
      if (keepx) then
         lenxbufp = 0
         lenxbufq = 0
      else
         lenxbufp  = nump*nbasis
         lenxbufq  = numq*nbasis
      endif
c$$$      lenscr    = pqmax*nbasis*npmax
      lenscr = max(nump*nbasis,numkl*numq)
      lenebuf   = nump*numq

      end
