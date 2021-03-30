c  main for complete frequency redistribtion (cylindrac, slab and spherical geometeries)

       implicit none
       include 'rad.inc'
       integer*4 i,j,k
       real*8 num, nump


       write(*,*)
       write(*,"('  PFR v1.2, G.J. Parker (c) 1997')")
       write(*,*)

       call StartUp		!read input files
       call calctransprob	!compute green's functions

       i = 0
       k = ndt
       
       
       if (prod_p) then	!there is a production, run to s.s.
         call writeint(0)
         do while (i.lt.num_ndt)
           i = i+1
           k = ndt
           call radtrap(k)
           if (mod(i,print_ndt).eq.0) call writeint(k*i)

           num = 0.d0
           do k = ikmax, 1, -1
             nump = 0.d0
             do j = 1, izmax
               nump = nump + nra(j,k,ci)
             enddo
             num = num + nump
           enddo
           write(*,"('# ndt ='i5,1p,'    total number =',
     ^             e13.5)") i,num
         enddo
       endif	
       


       if (decay_p) then	!finding fundamental...

         num = 0.d0
         do k = ikmax, 1, -1
           do j = 1, izmax
             nump = dzdp*(j-0.5d0)/zmaxdp
             if (cyl_p) then
c               nra(j,k,ci) = 1d30 * (j*j - (j-1)*(j-1))
       	nump = (nump*2.4048255577/3.0)**2
       	nra(j,k,ci) = ((((((0.00021d0*nump-0.039444d0)*nump
     ^       	  +0.0444479d0)*nump-0.3163866d0)*nump
     ^       	  +1.2656208d0)*nump-2.2499997d0)*nump+1.d0)
     ^                 * (j*j - (j-1)*(j-1)) * 1d30
             elseif(slab_p) then
c               nra(j,k,ci) = 1d30
       	nra(j,k,ci) = 1d30*sin(nump*pi)
             else
c               nra(j,k,ci) = 1d30*(j*j*j-(j-1)**3)
       	nra(j,k,ci) = 1d30*sin(nump*pi)/(nump*pi)
     ^                 	*(j*j*j-(j-1)**3)
             endif
             nra(j,k,ci) = nra(j,k,ci) * int_ls(k)
             num = num + nra(j,k,ci)
           enddo
         enddo
         write(*,"('initial number:',e13.5)") num
         call writeint(0)
         
         do while (num.gt.1d-22 .and. i.lt.num_ndt)
           i = i+1
           k = ndt
           call radtrap(k)
 
            nump = num
           num = 0.d0
           do k = ikmax, 1, -1
             do j = 1, izmax
               num = num + nra(j,k,ci)
             enddo
           enddo
           if (mod(i,print_ndt).eq.0) call writeint(i*ndt)
           write(*,"('# ndt ='i5,1p,'    Tvac/T =',e13.5,
     ^             '    (# =',e13.5,')')") i,
     ^             -log(num/nump)/dble(ndt),num
         enddo
       endif
       
       write(*,*)
       write(*,"('program successfully completed.')")
c      pause
       
       end
