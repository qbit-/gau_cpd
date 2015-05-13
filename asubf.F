      subroutine asubf (iout, iprint, nbasis, irw1, irw2, v, lenv)
      implicit none
      
      integer iout, iprint, lenv
      real*8  v(*)
      
      integer nbasis, irw1, irw2

      integer jend, blksz, jblksz, jblock1, jblock2, nbsq, ii
      
      jend  = 1
      nbsq  = nbasis*nbasis
      
      blksz = (lenv - jend) / (2*nbsq)
      if (blksz .gt. nbsq) blksz = nbsq
      
      jblock1 = jend
      jblock2 = jblock1 + blksz*nbsq
      jend    = jblock2 + blksz*nbsq
      call tstcor (jend-1, lenv, 'asubf: memory allocation')

      do 100 ii = 0, nbsq - 1, blksz
           if (ii + blksz .gt. nbsq) then
              jblksz = nbsq - ii
           else
              jblksz = blksz
           endif

           call fileio (2, -irw1, jblksz*nbsq, v(jblock1),ii*blksz*nbsq) 
           call fileio (2, -irw2, jblksz*nbsq, v(jblock2),ii*blksz*nbsq) 
           call daxpy (jblksz*nbsq,-1.0d0,v(jblock1),1,v(jblock2),1)
           call fileio (1, -irw2, jblksz*nbsq, v(jblock2),ii*blksz*nbsq) 
 100  continue
         

      end
