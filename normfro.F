      function normfro (iout, iprint, nbasis, irw, v, lenv)
      implicit none
      
      integer iout, iprint, lenv
      real*8  v(*)
      
      integer nbasis, irw

      integer jend, blksz, jblksz, jblock, nbsq, ii
      real*8  normfro
      real*8  dlange
      real*8  xx
      
      jend  = 1
      nbsq  = nbsq
      
      blksz = (lenv - jend) / (nbsq)
      if (blksz .gt. nbsq) blksz = nbsq
      
      jblock = jend
      jend   = jblock + blksz*nbsq
      call tstcor (jend-1, lenv, 'normfro: memory allocation')
      
      normfro = 0.0d0

      do 100 ii = 0, nbsq, blksz
           if (ii + blksz .gt. nbsq) then
              jblksz = nbsq - ii
           else
              jblksz = blksz
           endif

           call fileio (2, -irw, jblksz*nbsq, v(jblock), ii*blksz) 
           normfro = normfro + dlange ('f', nbsq, jblksz, v(jblock),
     &          nbsq, xx)
 100  continue
         
      end
