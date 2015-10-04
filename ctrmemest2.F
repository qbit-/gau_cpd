      subroutine ctrmemest2(iout, iprint, nbasis, ngrid, ndbf,
     &     keepov, keepints, keepx, numkl, nump, numq, lenwbuf, lenfbuf,
     &     lenxbufp, lenxbufq, lenscr, lenebuf, v, lenv) 

      implicit none 

C     estimates different sizes for the formation of the E matrix 
C
C     Arguments:
C     
C     [in] nbasis - number of basis functions
C     [in] ngrid  - number of grid points in the quadrature
C     [in] ndbf   - number of DBF functions if RI is used
C     [out] keepov   - if we keep an overlap for RI in memory 
C     [out] keepints - if fully in-core operation is possible
C     [out] keepx    - if we keep quadrature in memory
C     [out] numkl - step size over the KL index of integrals or DBF index of RI integrals. Here we allocate to iterate over whole counts of L only
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
      logical    keepov, keepints, keepx

C     Other

      integer    ndbf
      integer    nbsq, npmax, res, pqmax

C     functions 
      
      integer    nproc

      nbsq  = nbasis*nbasis
      npmax = nproc(0) 

C     We make a crude estimation if intergrals can be kept in memory

      keepov = .true.
      res    = lenv - ndbf*ndbf
      if (res .le. 0) call gauerr ('Insufficient memory in ctrmemest:0')

      keepx = .true.
      res   = lenv - ngrid*nbasis
      if (res .le. 0) call gauerr ('Insufficient memory in ctrmemest:1')
      
C     Fixme: this is completely wrong
      if ((lenv-ndbf*nbsq-2*nbsq-3*nbsq-npmax*nbsq-nbasis*ngrid).gt.0)
     &        then
         keepints = .true.
         res = res - ndbf*nbsq
      else
         keepints = .false.
      endif
      
      if (res .le. 0) call gauerr ('Insufficient memory in ctrmemest:2')

C     Just set fixed values in a dumb way for now

      nump = ngrid
      numq = ngrid
      pqmax = max(nump, numq)

C     Scratch buffers for loading grid blocks

      res = lenv - nump*nbasis - numq*nbasis 
     &        - nump*numq*npmax - ndbf*pqmax*(npmax + 1)
         
      if (res .lt. 0)
     &     call gauerr('Insufficient memory in ctrmemest:3')

C     Keep W buffer and F buffer. For RI, 2 F buffers are needed

      numkl = res / (nbsq + 2*pqmax)
      
      if (numkl .gt. ndbf) numkl = ndbf
      
C     
      if (numkl .lt. 2)
     &     call gauerr('Insufficient memory in ctrmemest:4')

C     Now assign all sizes

      lenwbuf = numkl*nbsq
      lenfbuf = numkl*nump
      lenxbufp  = nump*max(nbasis,ndbf)
      lenxbufq  = numq*max(nbasis,ndbf)
      lenscr    = pqmax*max(nbasis,ndbf)*npmax
      lenebuf   = nump*numq

      end
