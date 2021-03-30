c       routine to update resonance atom dens

       subroutine radtrap(numst)

       implicit none
       include 'rad.inc'
       integer*4 i,l,cinext,i0,numst,iz_min,j
c       real*8 numi
       real numi


       iz_min = 1
       if (slab_p) iz_min = 0
       
       cinext = 3 - ci

c now initialize histogram of escaping photons:
       do i = 1, ikmax
         nra(0,      i,cinext) = 0.
         nra(izmax+1,i,cinext) = 0.
       enddo

c do numst steps...
       do l=1,numst

c first do spatial motion (ie propagator Q):

         do j = 1, ikmax
           do i = 1, izmax
             numi = nra(i,j,ci) + prra(i,j)   !those there + production...
             nra(i,j,ci) = 0.
             do i0 = iz_min, izmax+1
               nra(i0,j,cinext) = pa(i0,i,j) * numi + nra(i0,j,cinext)
             enddo
           enddo
         enddo
         
c now do the frequency redistribution (ie propagator J):

         do j = 1, ikmax
           do i = 1, izmax
             numi = nra(i,j,cinext)
             nra(i,j,cinext) = 0.
             do i0 = 1, ikmax
               nra(i,i0,ci) = nra(i,i0,ci) + numi * pac(i0,j)
             enddo
           enddo
         enddo
         
       enddo
       
c now update the histogram of escaping photons:
       do i = 1, ikmax
         nesc(0,i) = nesc(0,i) + nra(0,      i,cinext)
         nesc(1,i) = nesc(1,i) + nra(izmax+1,i,cinext)
       enddo
       
       return
       end
