       subroutine CalcTransProb
       
       include 'rad.inc'
       
       
       
       if (.not. read_p) then
         call CalcPFR
       
         if (cyl_p) then
           call calctransprob_cyl
         elseif (slab_p) then
           call calctransprob_slab
         else
           call calctransprob_sph
         endif
       endif
       
       if (save_p) then
         write(*,*)
         write(*,"('  saving J & Q...')")
         call save_new
         write(*,"('  done saving J & Q')")
         write(*,*)
       endif
       
       return
       end








       subroutine CalcPFR
c routine to calculate frequence redistribution propagator (ie J)
c current kernels supported:
c	1) strict Jefferies-White (JW)
c	2) JW but w/ doppler instead of delta function
c	3) angle averaged exact

       implicit none
       include 'rad.inc'
       integer*4 i,j
       real*8 frac0,f0,f1
       real*8 af(kid)
       real voigt, doppler, rint_dopp, erfc
       
       real*8 p0(kid),p1(kid),p2(kid),p3(kid)
       
       write(*,*)
       
       do i = 1, kid
         p0(i) = 0.d0
         p1(i) = 0.d0
         p2(i) = 0.d0
         p3(i) = 0.d0
       enddo

       if (cfr_dopp_p .or. cfr_lrtz_p .or. cfr_vgt_p) then
         if ( cfr_dopp_p ) then
           write(*,"('  computing frequency redistribution propagator'
     ^	  	,' J: (CFR w/ Doppler)')")
         elseif ( cfr_lrtz_p ) then
           write(*,"('  computing frequency redistribution propagator'
     ^	  	,' J: (CFR w/ Lorentzian)')")
         else
           write(*,"('  computing frequency redistribution propagator'
     ^	  	,' J: (CFR w/ Voigt)')")
         endif

c doing complete frequency redistribution (CFR)
c [do you always use a sledge hammer to swat a fly?]

c by definition, final frequency independent of initial frequency.
c probability of going to any specific frequency cell is simply 
c equal to the integral of the line shape over that cell. did that 
c integral already, so it makes this trivial... (it's even normalized!)
         
         do i = 1, ikmax	!initial frequency
           do j = 1, ikmax  !final frequency
             pac(j,i) = int_ls(j)
           enddo
           write(*,"('    frequency cell ',i3,' done')") i	    
         enddo
       
c pure Doppler line shape... redistribution function is erfc( x>)...	
       elseif (pure_dopp_p) then
         write(*,"('  computing frequency redistribution propagator'
     ^	  	,' J: (pure Doppler)')")

         do i = 1, ikmax
           f0    = min( df(1,i), df(2,i) )
           f1    = max( df(1,i), df(2,i) )
           
           p0(i) = f1 * erfc(f1) - f0 * erfc(f0) +
     ^	    	   (exp(-f0*f0) - exp(-f1*f1)) / sqrt(pi)
            p1(i) = (fv(i) - f0) * erfc(fv(i)) +   
     ^              f1 * erfc(f1) - fv(i) * erfc(fv(i)) +
     ^		   (exp(-fv(i)*fv(i)) - exp(-f1*f1)) / sqrt(pi)
           p2(i) = erfc(fv(i))
c            p1(i) = .5d0 * ( erfc(fv(i)) * (fv(i) - f0) +
c     ^		   f1 * erfc(f1) - fv(i) * erfc(fv(i)) +
c     ^		   (exp(-fv(i)*fv(i)) - exp(-f1*f1)) / sqrt(pi) )
c            p2(i) = .5d0 * ( erfc(fv(i)) * (fv(i) - f0) +
c     ^		   f1 * erfc(f1) - fv(i) * erfc(fv(i)) +
c     ^		   (exp(-fv(i)*fv(i)) - exp(-f1*f1)) / sqrt(pi) )
         enddo
         
         do i = 1, ikmax   !initial frequency cell
           
           frac0 = 0.d0
           do j = 1, ikmax		!final frequency cell
             if (j .lt. i) then		! final frequency < initial
               p3(j) = p2(i) * abs( df(2,j) - df(1,j) )
             elseif (j .gt. i) then	! final frequency > initial
               p3(j) = p0(j)
             else			  !same frequency cell...
               p3(j) = p1(j) 
             endif
             frac0 = frac0 + p3(j)
           enddo
           
           frac0 = 1.d0 / frac0
           do j = 1, ikmax
             pac(j,i) = p3(j) * frac0
           enddo

           write(*,"('    frequency cell ',i3,' done')") i	    
         enddo
         

c Jeffries-White: either w/ delta function or shifted doppler...
       elseif (jw_p) then
         if ( jw_dopp_p ) then
           write(*,"('  computing frequency redistribution propagator'
     ^	  	,' J: (Jeffries-White w/ shifted Doppler)')")
         else
           write(*,"('  computing frequency redistribution propagator'
     ^	  	,' J: (Jefferies-White)')")
         endif
       
c   integral over Voigt- already done...
         do i = 1, ikmax
           p0(i) = int_ls(i) * Pcoll
         enddo

c   integral over Doppler
         p1(ikmax) = rint_dopp(df(1,ikmax),1.d32)
         frac0 = p1(ikmax)
         do i = ikmax-1, 1, -1
           p1(i) = rint_dopp(df(1,i),df(2,i))
           frac0 = frac0 + p1(i)
         enddo
         frac0 = 1.d0 / frac0
         do i = 1, ikmax
           p1(i) = p1(i) * frac0 * (1.d0 - pcoll)
         enddo
         
c    critical function
         do i = 1, ikmax
           af(i) = max(0.d0, min(1.d0,
     ^	    	1.d0 - doppler(fv(i)) / voigt(real(fv(i)),real(av))))
         enddo
         
c ok, construct J: for each initial cell:
         do i = 1, ikmax
         
           frac0 = 0.d0
           do j = ikmax, 1, -1	!for each final cell

             p2(j) = (1.d0 - af(i)) * p1(j) + p0(j)
             
             if (jw_dopp_p) then		!shifted doppler
               f0 = fv(i) - df(2,j)
               f1 = fv(i) - df(1,j)
               if (f0*f1.ge.0.d0) then
                 p2(j) = p2(j) + (1.d0 - Pcoll) * af(i) *
     ^	        	rint_dopp(abs(f0),abs(f1))
               else
                 p2(j) = p2(j) + (1.d0 - Pcoll) * af(i) *
     ^	          	(rint_dopp(0.d0,abs(f0)) + 
     ^			 rint_dopp(0.d0,abs(f1)))
               endif
               f0 = -fv(i) - df(2,j)
               f1 = -fv(i) - df(1,j)
               p2(j) = p2(j) + (1.d0 - Pcoll) * af(i) *
     ^	        	rint_dopp(abs(f0),abs(f1))
             elseif (i.eq.j) then		!delta function
               p2(j) = p2(j) + (1.d0 - Pcoll) * af(i)
             endif
             frac0 = frac0 + p2(j)
           enddo
           
           frac0 = 1.d0 / frac0
           do j = ikmax, 1, -1
             pac(j,i) = p2(j) * frac0
           enddo
           
           write(*,"('    frequency cell ',i3,' done')") i
         enddo
         
       elseif (exact_incoh_p) then
c exact J w/ incoherent scattering in frame of the atom:
         write(*,"('  computing frequency redistribution propagator'
     ^	  	,' J: (exact w/ incoherent scattering)')")
         call compute_incoh
       else
c exact J w/ coherent scattering in frame of the atom:
         write(*,"('  computing frequency redistribution propagator'
     ^	  	,' J: (exact w/ coherent scattering)')")
         call compute_coh
       endif

       write(*,"('  done computing J')")
       write(*,*)
       
       return
       end







       subroutine compute_coh
c routine to calculate frequence redistribution propagator (ie J)
c for coherent scattering in atomic rest frame. see AF Molisch, et al, 
c JQRST, v53, n3, pp 269-275, (1995), eqn 10. for a better discussion see
c JT Jefferies, Spectral Line Formation, Blaisdell, Waltman (1968), 
c chapter 5.

       implicit none
       include 'rad.inc'
       integer*4 i,j
       real*8 frac0,x1,x2,xi,norml,pi3_2
       
       real*8 p0(kid),p1(kid)
       real*8 comp_ex
       real voigt
       
       real*8 x(200), w(200)
       integer n
       
       n = 150
       call find_gaulag( x, w, n)
       
       pi3_2 = sqrt( (acos(-1.d0))**3 )
       
c   integral over Voigt- already done...
       do i = 1, ikmax
         p0(i) = int_ls(i) * Pcoll
       enddo
       
       do i = 1, ikmax		!initial frequency
         frac0 = 0.d0
         xi = fv(i)
         norml = 1.d0 / ( pi3_2 * voigt(real(fv(i)),real(av)) )
         do j = 1, ikmax		!final frequency
           x1 = df(1,j)
           x2 = df(2,j)
           if ( i .ne. j ) then
             p1(j) = comp_ex(xi, x1, x2,av,x,w,n,norml) +
     ^	              comp_ex(xi,-x2,-x1,av,x,w,n,norml)
           else
             p1(j) = comp_ex(xi, min(x1,x2), xi,av,x,w,n,norml) +
     ^	              comp_ex(xi, xi, max(x1,x2),av,x,w,n,norml) +
     ^		      comp_ex(xi,-x2,-x1,        av,x,w,n,norml)
           endif
           frac0 = frac0 + p1(j)
         enddo
         
         frac0 = 1.d0 / frac0
         do j = 1, ikmax
           pac(j,i) = p1(j) * frac0 * (1.d0 - pcoll) + p0(j)
         enddo
         
         write(*,"('    frequency cell ',i3,' done')") i
       enddo
       
       return
       end
       


       
       
       function comp_ex(xi,x1,x2,a,x,w,n,nrml)
c function return the 
       
       implicit none
       real*8 comp_ex, xi, x1, x2, a, nrml
       real*8 xl, xu, frac0, frac1, frac2, frac3
       real*8 u1, u2, du, u, xx, du0
       integer num_u, i, nn, num_u0, nrf
       
       integer n
       real*8 x(n), w(n)
       
       nrf   = 0

       num_u = 25
       
       frac2 = 0.d0
       frac3 = 0.d0
       
       du    = ( max( x2, x1) - min( x2, x1) ) / num_u
       
       frac1 = 0.d0

       if ( 0.5 * (x2 + x1) .lt. xi ) then
         u = min( x2, x1 )
       else
         u = max( x2, x1 )
         du = - du
       endif
       du0    = abs(du)
       num_u0 = num_u

       do while ( nrf .lt. 20 )
         do i = 1, num_u
         
           u   = u + du
           if ( du .lt. 0. ) then
             u = max( u, xi )
           else
             u = min( u, xi )
           endif
         
           xl = min( u, xi )
           xu = max( u, xi )
           u1 = 0.5 * ( xu - xl )
           u2 = 0.5 * ( xu + xl )
  
           frac0 = 0.d0
           do nn = n, 1, -1
             xx    = sqrt( x(nn) )
             frac0 = frac0 + w(nn) * exp( - 2.d0 * xx * u1 ) *
     ^	      	 ( atan( (u2 + xx) / a ) - atan( (u2 - xx) / a ) ) / xx
           enddo
           frac1 = frac1 + frac0 * exp( - u1 * u1 )
         
         enddo
       
         frac0 = abs( (frac1+frac2) * 0.5d0 * du0 ) * nrml
         
c	  write(*,"(i2,i8,1p,4e14.7)") nrf,num_u0,
c     ^	  	frac3,frac0,abs(1.d0-frac3/frac0),abs(frac3-frac0)
     
         if ( abs(frac3 - frac0) .gt. 1.d-4 .or. nrf .eq. 0 ) then
           frac2  = frac1 + frac2
           frac1  = 0.d0
  	    frac3  = frac0
         
           num_u0 = num_u0 * 2
           du0    = ( max( x2, x1) - min( x2, x1) ) / num_u0
           du     = 2.d0 * du0

           if ( 0.5 * (x2 + x1) .lt. xi ) then
             u    = min( x2, x1 ) - du0
           else
             u    = max( x2, x1 ) + du0
             du   = - du
           endif
         
           if (nrf .ne. 0) num_u = num_u * 2
           nrf    = nrf + 1
         else
           nrf    = 2000
         endif
       enddo
       
       comp_ex = frac0

       return
       end





       subroutine compute_incoh
c routine to calculate frequence redistribution propagator (ie J)
c for incoherent scattering in atomic rest frame. see AF Molisch, et al, 
c JQRST, v53, n3, pp 269-275, (1995), eqn 14. for a better discussion see
c JT Jefferies, Spectral Line Formation, Blaisdell, Waltman (1968), 
c chapter 5.

       implicit none
       include 'rad.inc'
       integer*4 i,j
       real*8 frac0,incomp_ex
       real*8 xi,x1,x2,norml,pi5_2
       
       real*8 p0(kid),p1(kid)
       real voigt

       real*8 x(200), w(200)
       integer n
       
       n     = 150
       call find_gaulag( x, w, n)

       pi5_2 = sqrt( (acos(-1.d0))**5 )
       
c   integral over Voigt- already done...
       do i = 1, ikmax
         p0(i) = int_ls(i) * Pcoll
       enddo

       do i = 1, ikmax		!initial frequency
         frac0 = 0.d0
         xi = fv(i)
         norml = 1.d0 / ( pi5_2 * voigt(real(fv(i)),real(av)) )
         do j = 1, ikmax		!final frequency
           x1 = df(1,j)
           x2 = df(2,j)
           if ( i .ne. j ) then
             p1(j) = incomp_ex(xi, x1, x2,av,x,w,n,norml) +
     ^	              incomp_ex(xi,-x2,-x1,av,x,w,n,norml)
           else
             p1(j) = incomp_ex(xi, min(x1,x2), xi,av,x,w,n,norml) +
     ^	              incomp_ex(xi, xi, max(x1,x2),av,x,w,n,norml) +
     ^		      incomp_ex(xi,-x2,-x1,        av,x,w,n,norml)
           endif
           frac0 = frac0 + p1(j)
         enddo
         
         frac0 = 1.d0 / frac0
         do j = 1, ikmax
           pac(j,i) = p1(j) * frac0 * (1.d0 - pcoll) + p0(j)
         enddo
         
         write(*,"('    frequency cell ',i3,' done')") i
       enddo

       return
       end

       
       
       function incomp_ex(xi,x1,x2,a,x,w,n,nrml)
c function return the 
       
       implicit none
       real*8 incomp_ex, xi, x1, x2, a, nrml
       real*8 frac0, u1, u2
       integer nn
       
       real*8 y1, y2, v1, v2, w1, w2
       
       integer n
       real*8 x(n), w(n)
       
       y1    = min( x1, x2 )
       y2    = max( x1, x2 )
       
       frac0 = 0.d0

       do nn = n, 1, -1
         u2  = x(nn)
         u1  = sqrt( u2 )
         v2  = (y2 + u1) / a
         v1  = (y1 + u1) / a
         w2  = (y2 - u1) / a
         w1  = (y1 - u1) / a
         frac0 = frac0 + a * w(nn) *
     ^		( v2 * atan( v2 ) - v1 * atan( v1 )
     ^		- w2 * atan( w2 ) + w1 * atan( w1 ) + .5d0 * log( 
     ^            (1.d0 + v1 * v1) * (1.d0 + w2 * w2) /
     ^	        ( (1.d0 + v2 * v2) * (1.d0 + w1 * w1) )) ) *
     ^	        ( atan( (xi + u1) / a ) - atan( (xi - u1) / a ) ) / u1
       enddo
       
       incomp_ex = abs( frac0 * 0.5d0 * nrml )

       return
       end






       subroutine CalcTransProb_slab

       implicit none
       include 'rad.inc'
       integer*4 i,j,nn,k
       real*8 frac0
       real*8 zi,zi1,zi2,zz
       real*8 kval
       real r4
       
       real*8 mb(0:zid),ppa(0:zid+1),weight,dA,expint2

       write(*,*)
       write(*,"('  computing spatial propagator Q: (slab)')")

c from each initial cell (both frequency and spatial)
       do k = 1, ikmax
         kval = kv(k)
         do i=1,izmax
         
           do j = 0, izmax+1
             pa(j,i,k) = 0.
           enddo
           
           dA = dzdp / dble(nave)
           weight = 1.d0 / dble(nave)
           zi2 = dzdp * (i - 1)
           
           do nn = 1, nave
           
             zi1 = zi2
             zi  = dA*0.5d0 + zi1
             zi2 = dA + zi1
          
c to each final boundary
             do j=0,izmax
               zz = abs(kval*(dzdp * j - zi))
c	        mb(j) = sign(1, j-i)*expint2(zz)*0.5d0
               mb(j) = 	        expint2(zz)*0.5d0
             enddo
           
             ppa(0) = abs(mb(0))
             ppa(izmax+1) = mb(izmax)
             frac0 = ppa(0) + ppa(izmax+1)
             do j = 1, izmax

c	        ppa(j) = mb(j-1) - mb(j)
c	        if (j.eq.i) ppa(j) = 1.d0 + ppa(j)

               ppa(j) = abs(mb(j-1) - mb(j))
               if (j.eq.i) ppa(j) = 1.d0 - mb(j-1) - mb(j)
               
               ppa(j) = max( 0.d0, ppa(j))
               frac0 = frac0 + ppa(j)
             enddo
           
             frac0 = 1.d0 / frac0
             do j = 0, izmax+1
               pa(j,i,k) = pa(j,i,k) + weight*ppa(j)*frac0
             enddo
           
           enddo
           
           frac0 = 0.d0
           do j = 0, izmax+1
             frac0 = frac0 + pa(j,i,k)
           enddo
           
           frac0 = 1.d0 / frac0
           r4 = 0.
           do j = 0, izmax+1	!normalize...
             pa(j,i,k)=pa(j,i,k) * frac0
             r4 = max(r4, pa(j,i,k))
           enddo
           if (r4 .eq. 1.) then
             do j = 0, izmax+1
               if (pa(j,i,k) .ne. r4) then
                 pa(j,i,k) = 0.
               else
                 pa(j,i,k) = 1.
               endif
             enddo
           endif
           
         enddo
         write(*,"('    frequency cell ',i3,' done')") k
       enddo

       write(*,"('  done computing spatial propagator Q.')")
       write(*,*)

       return
       end
       





       real*8 function expint2(x)
c return E2(x): exponetial integral of order 2

       real*8 x

       if (x.lt.1.d0) then
         expint2=exp(-x)-x*(-log(x)+((((
     ^		0.00107857d0*x-0.00976004)*x+0.05519968d0)*x-
     ^		0.24991055d0)*x+0.99999193d0)*x-0.57721566d0)
       else
         expint2=exp(-x)*(1.d0 -((((
     ^		x+8.5733287401d0)*x+18.0590169730d0)*x+8.6347608925d0)*
     ^		x+0.2677737343d0)/((((
     ^		x+9.5733223454d0)*x+25.6329561486d0)*x+
     ^		  21.0996530827d0)*x+3.9584969228d0))
       endif

       return
       end





       
       

       subroutine CalcTransProb_sph

       implicit none
       include 'rad.inc'
       integer*4 i,j,nn,k
       real*8 frac0
       real*8 zi,zf,zi1,zi2,kval
       real r4
       
       real*8 mb(0:zid),ppa(0:zid+1),weight,dA,expint2

       write(*,*)
       write(*,"('  computing spatial propagator Q: (spherical)')")

       mb(0) = 0.d0
c from each initial cell (both frequency and spatial)
       do k = 1, ikmax
         kval = kv(k)
         do i=1,izmax

           do j = 0, izmax+1
             pa(j,i,k) = 0.
           enddo

           dA = dzdp*dzdp*dzdp*(i*i*i-(i-1)**3) / dble(nave)
           weight = 1.d0 / dble(nave)

c	    dA = dzdp / dble(nave)
c	    weight = 1.d0 / dble(nave)

  	    zi2 = dzdp * (i - 1)
           
           do nn = 1, nave
           
             zi1 = zi2

             zi  = (dA*0.5d0 + zi1*zi1*zi1)**(1.0/3.0)
             zi2 = min((dA   + zi1*zi1*zi1)**(1.0/3.0),i*dzdp)

c	      zi  = dA*0.5d0 + zi1
c	      zi2 = min(dA   + zi1, i*dzdp)
             
             if (nave.eq.2) then
               if (nn.eq.1) then
                 weight = 0.8d0
                 zi = dzdp*(i-0.5d0)
               else
                 weight = 0.2d0
                 zi = dzdp*( 0.5d0*(i*i*i+(i-1.d0)**3) )**(1.0/3.0)
               endif
             elseif (nave.eq.1) then
               zi1 = dzdp * (i-1)
               zi2 = dzdp * i
               zi  = 0.75d0 * (zi2**4 - zi1**4) / (zi2**3 - zi1**3)
             endif
          
c to each final boundary
             do j = 1, izmax
               zf = dzdp * j
               mb(j) =   abs(
     ^	          (exp(-kval*abs(zf-zi))-exp(-kval*(zf+zi)))/
     ^	          (4.d0*kval*zi) + (zf*zf-zi*zi)/(4.d0*zi) * (
     ^	          (expint2(kval*abs(zf-zi)))/abs(zf-zi) - 
     ^	          (expint2(kval*(zf+zi)))/(zf+zi))   )
             enddo
           
             ppa(0)       = 0.d0
             ppa(izmax+1) = mb(izmax)
             frac0 = ppa(izmax+1)
             do j = 1, izmax
               ppa(j) = abs( mb(j-1) - mb(j) )
               if (j.eq.i) ppa(j) = 1.d0 - mb(j) - mb(j-1)

c	        ppa(j) = mb(j-1) - mb(j)
c	        if (j.eq.i) ppa(j) = 1.d0 + ppa(j)

       	 ppa(j) = max(0.d0, ppa(j))
               frac0 = frac0 + ppa(j)
             enddo
           
             frac0 = 1.d0 / frac0
             do j = 0, izmax+1
               pa(j,i,k) = pa(j,i,k) + weight*ppa(j)*frac0
             enddo
           
           enddo
           
           frac0 = 0.d0
           do j = 0, izmax+1
             frac0 = frac0 + pa(j,i,k)
           enddo
           frac0 = 1.d0 / frac0
           r4 = 0.
           do j = 0, izmax+1	!normalize...
             pa(j,i,k)=pa(j,i,k) * frac0
             r4 = max( r4, pa(j,i,k))
           enddo
           if (r4 .eq. 1.) then
             do j = 0, izmax+1
               if (pa(j,i,k) .ne. r4) then
                 pa(j,i,k) = 0.
               else
                 pa(j,i,k) = 1.
               endif
             enddo
           endif
           
         enddo
         write(*,"('    frequency cell ',i3,' done')") k
       enddo

       write(*,"('  done computing spatial propagator Q.')")
       write(*,*)

       return
       end
       
       
       
       
       

       subroutine CalcTransProb_cyl

       implicit none
       include 'rad.inc'
       integer*4 i,j,nn,j0,i0,k
       real*8 frac0
       real*8 zi,zf,zi1,zi2,kval
       real*8 int_srf, int_vol
       
       real*8 mb(0:zid+1),weight,dA
       real*8 ppa_s(0:zid+1)  !, ppa_v(0:zid+1)
       real*8 pa_s(0:zid+1), pa_v(0:zid+1)
       real r4

       

       call loadgauss()
       
       write(*,*)
       write(*,"('  computing spatial propagator Q:(cylindrical)')")

       mb(0) = 0.d0
c from each initial cell (both frequency and spatial)
       do k = 1, ikmax
         kval = kv(k)
         do i=1,izmax
           i0 = i
           
           do j = 0, izmax+1
             pa_s(j) = 0.d0
             pa_v(j) = 0.d0
           enddo

           dA = dzdp * (i*i - (i-1)*(i-1)) * dzdp / dble(nave)
           weight = 1.d0 / dble(nave)
           zi2 = max(0.d0,dzdp * (i-1))

c	    dA = dzdp / dble(nave)
c	    weight = 1.d0 / dble(nave)
c  	    zi2 = dzdp * (i - 1)
           
           do nn = 1, nave
           
             zi1 = zi2

             zi = sqrt(dA * 0.5d0 + zi1*zi1)
             zi2 = min(i*dzdp, sqrt(dA+zi1*zi1) )

c	      zi  = dA*0.5d0 + zi1
c	      zi2 = min(dA   + zi1, i*dzdp)

             if (nave.eq.2) then
               if (nn.eq.1) then
                 weight = 0.8d0
                 zi = dzdp*(i-0.5d0)
               else
                 weight = 0.2d0
                 zi = dzdp* sqrt( 0.5d0*(i*i+(i-1.d0)**2) )
               endif
             elseif (nave.eq.1) then
               zi  = 2.d0 * (i*i*i-(i-1)*(i-1)*(i-1)) /
     ^		   ( 3.d0 * (i*i - (i-1)*(i-1) )) * dzdp
             endif
          
c to each final boundary
cv	      frac0 = 0.d0
             do j=1,izmax
               j0 = j
               zf = dzdp * j
               mb(j) = int_srf(zi,zf,kval)
csa	        mb(j) = abs(int_srf(zi,zf,kval))
cv	        ppa_v(j) = int_vol(zi,zf,i0,j0,kval,dzdp)
cv	        frac0 = frac0 + ppa_v(j)
             enddo
             
cv	      mb(izmax) = int_srf(zi,dzdp*j,kval)
             
cv	      ppa_v(0) = 0.d0
cv	      ppa_v(izmax+1) = abs(mb(izmax))
cv	      frac0 = 1.d0 / (frac0 + ppa_v(izmax+1))
cv	      do j = 1, izmax+1
cv	        pa_v(j) = pa_v(j) + weight * ppa_v(j) * frac0
cv	      enddo
           
             ppa_s(0) = 0.d0
             ppa_s(izmax+1) = abs(mb(izmax))
             frac0 = ppa_s(izmax+1)
             do j = 1, izmax
               ppa_s(j) = mb(j-1) - mb(j)
               if (j.eq.i) ppa_s(j) = 1.d0 + ppa_s(j)
csa	        ppa_s(j) = abs(mb(j-1) - mb(j))
csa	        if (j.eq.i) ppa_s(j) = 1.d0 - mb(j-1) - mb(j)
               
cv	        if (ppa_s(j) .le. 0.d0) ppa_s(j) = ppa_v(j)
       
       	 if (ppa_s(j) .le. 1e-5) then
       	   j0 = j
       	   zf = dzdp * j
       	   ppa_s(j) = int_vol(zi,zf,i0,j0,kval,dzdp)
       	 endif
       	 
       	 ppa_s(j) = max(0.d0, ppa_s(j))
               frac0 = frac0 + ppa_s(j)	        
             enddo
           
             frac0 = 1.d0 / frac0
             do j = 0, izmax+1
               pa_s(j) = pa_s(j) + weight * ppa_s(j) * frac0
             enddo
           
           enddo
           
           frac0 = 0.d0
           do j = 0, izmax+1
cv	      if (pa_s(j).eq.0.d0) pa_s(j) = pa_v(j)
cvs	      pa(j,i,k) = 0.5d0 * (pa_s(j) + pa_v(j))

cv	      pa(j,i,k) = pa_v(j)	      
             pa(j,i,k) = pa_s(j)
             
             frac0 = frac0 + pa(j,i,k)
           enddo
           
c	    write(*,*) i,k,frac0
           frac0 = 1.d0 / frac0
           r4 = 0.
           do j = 0, izmax+1	!normalize...
             pa(j,i,k) = pa(j,i,k) * frac0
             r4 = max( r4, pa(j,i,k))
           enddo
           if (r4 .eq. 1.) then
             do j = 0, izmax+1
               if (pa(j,i,k) .ne. r4) then
                 pa(j,i,k) = 0.
               else
                 pa(j,i,k) = 1.
               endif
             enddo
           endif
           
         enddo
         write(*,"('    frequency cell ',i3,' done')") k
       enddo

       write(*,"('  done computing spatial propagator Q.')")
       write(*,*)

       return
       end






       function int_srf(zi,zf,kval)
c return the following surface integral: sign(zf-zi) *
c kval * zf * (zf - zi * cos(phi)) [ K_1(kval*q) - Ki_1(kval(q) ] /
c (q * pi), where q*q = zf*zf + zi*zi - 2*zi*zf*cos(phi) and the 
c integral is over phi from 0 to pi. 
c uses gaussian/legendre quadrature...
       
       real*8 int_srf, zi, zf, kval
       real*8 hb2, h2b2, mu, q, sum1, k1,ki1
       integer i
       
       real*8 xphi_s(300),wphi_s(300)
       real*8 xr_v(20,200),wr_v(20,200)
       real*8 xphi_v(300,200),wphi_v(300,200)
       integer phimax,rmax
       
       common /gausscom/ xphi_s, wphi_s, xr_v, wr_v,
     ^		xphi_v, wphi_v, phimax, rmax
       	
       	
       hb2 = -2.d0*zi*zf
       h2b2 = zf*zf + zi*zi
       
       sum1 = 0.d0
       do i = phimax, 1, -1
         mu = xphi_s(i)
         q = sqrt(max(0.d0, h2b2 + hb2*mu))
         sum1 = sum1 + 
     ^	  	(zf - zi*mu)*(k1(kval*q)-ki1(kval*q))*wphi_s(i)/q
       enddo
       
       int_srf = sign(1.d0, zf-zi) * abs(
     ^		zf * kval * sum1 / 3.1415926535897932385d0)

       return
       end





       function int_vol(zi, zf, izi, izf, kval, dz)
c return the following volue integral: 
c kval * r * Ki_1(kval(q) ] / (2 * pi q)
c where q*q = r*r + zi*zi - 2*zi*r*cos(phi) and the 
c integral is over phi from 0 to pi and over r in the radial cell. 
c if initial and final cells are the same, then phi integration is reduced
c and we add on an anayltic solution.
c uses gaussian/legendre quadrature... do phi integral inside r, since 
c phi integral will have many more terms, in general
       
       real*8 int_vol, zi, zf, kval, dz
       integer izf,izi
       real*8 mu, q, sum1, zff, sum0, ki1, k1
       integer i,j

       real*8 xphi_s(300),wphi_s(300)
       real*8 xr_v(20,200),wr_v(20,200)
       real*8 xphi_v(300,200),wphi_v(300,200)
       integer phimax,rmax
       
       common /gausscom/ xphi_s, wphi_s, xr_v, wr_v,
     ^		xphi_v, wphi_v, phimax, rmax
       	
       
       sum1 = 0.d0	
       if (izf .ne. izi) then
         do j = 1, rmax
           sum0 = 0.d0
           zff  = xr_v(j,izf)
           do i = phimax, 1, -1
             mu = xphi_s(i)
             q = sqrt(max(0.d0, zi*zi + zff*zff - 2.d0*zi*zff*mu))
             sum0 = sum0 + zff * ki1(kval*q) * wphi_s(i)/ q
           enddo
           sum1 = sum1 + sum0 * wr_v(j,izf)
         enddo
         int_vol = kval * sum1 / 3.1415926535897932385d0
       else
         do j = 1, rmax
           sum0 = 0.d0
           zff  = xr_v(j,izf)
           do i = phimax, 1, -1
             mu = xphi_v(i,izf)
             q = sqrt(max(0.d0, zi*zi + zff*zff - 2.d0*zi*zff*mu))
             sum0 = sum0 + zff * ki1(kval*q) * wphi_v(i,izf)/ q
           enddo
           sum1 = sum1 + sum0 * wr_v(j,izf)
         enddo
         sum0 = kval * dz
         int_vol = max(0.d0, kval * sum1 / 3.1415926535897932385d0 +
     ^	  	1.d0 + sum0 * ( ki1(sum0) - k1(sum0) ) )
       endif

       return
       end





       subroutine loadgauss()

       implicit none
       include 'rad.inc'
       real*8 pi8,xl,xm,tmp(20),tmp1(20)
       integer i,j
       
       real*8 xphi_s(300),wphi_s(300)
       real*8 xr_v(20,200),wr_v(20,200)
       real*8 xphi_v(300,200),wphi_v(300,200)
       integer phimax,rmax
       
       common /gausscom/ xphi_s, wphi_s, xr_v, wr_v,
     ^		xphi_v, wphi_v, phimax, rmax


       write(*,*)
       write(*,"('  computing Gaussian quadrature terms:')")

       pi8 = acos(-1.d0)

c find phi & weights for the surface integral. we want finest 
c resolution about phi = 0, so range runs from 0 to 2 pi. the 
c integeral is symmetric, so only get half of the phi values...
       
       phimax = 300
       phimax = n_phi_gauss
       call find_gauleg(0.d0, 2.d0*pi8, 
     ^		xphi_s, wphi_s, phimax, .true.)
       
c ok, take the cos of phi since that's what
c is used in the surface integral...
       do i = 1, phimax
         xphi_v(i,1) = xphi_s(i) / pi8 - 1.d0  !maps 0 -> -1., pi -> 0
         wphi_v(i,1) = wphi_s(i) / pi8	     !divide the range dependence out of weights
         xphi_s(i)   = cos( xphi_s(i) )
       enddo
       		
       	
c find r for the volume integral. we don't need many points.
       rmax = 1
       rmax = n_r_gauss
       if (rmax .ne. 1) then
         call find_gauleg(-1.d0, 1.d0, tmp, tmp1, rmax, .false.)
       else
         tmp(1)  = 0.d0
         tmp1(1) = 2.d0
       endif
       
c for each final cell:
       do i = izmax, 1, -1
       
c			find the r values and correct weigths, store
         xm = dzdp * (i - 0.5d0)
         xl = dzdp * 0.5d0
         do j = 1, rmax
           xr_v(j,i) = tmp(j) * xl + xm
           wr_v(j,i) = tmp1(j) * xl
         enddo

c			find the range of angle, angle values and correct weights, store	  
         xl = 2.d0 * atan(0.25d0/(i - 0.5d0))
         xm = pi8
         xl = pi8 - xl
         do j = 1, phimax
           xphi_v(j,i) = cos( xm + xl*xphi_v(j,1) )
           wphi_v(j,i) = wphi_v(j,1) * xl
         enddo
         
       enddo

       write(*,"('  done computing terms')")
       write(*,*)

       return
       end
       
       





       subroutine find_gauleg(x1,x2,x,w,n1,sym_p)
       
c the following program computes the abscissas and weight factors for 
c gaussian integration.  this copied from `numerical recipes (fortran)' 
c by w.h.press, b.p.flannery, s.a.teukolsky and w.t.vetterling,
c Cambridge press, 1989,1990, page 125 (called GAULEG, ie weight function =1)

       implicit real*8 (a-h,o-z)
       real*8 x1,x2,x(n1),w(n1)
       logical sym_p
       parameter (eps=3.d-14)

c n=number of terms... x1(x2) is the lower(upper) range of integration

       n = n1
       if (sym_p) n = 2 * n1

       m  = (n+1) / 2
       xm = 0.5d0 * (x2+x1)
       xl = 0.5d0 * (x2-x1)
       do i = 1, m
         z=cos(3.1415926535897932388462643d0*(i-0.25d0)/(n+0.5d0))
 1	  continue
           p1=1.d0
           p2=0.d0
           do j = 1, n
             p3 = p2
             p2 = p1
             p1 = ((2.0d0*j-1.0d0)*z*p2-(j-1.d0)*p3)/j
           enddo
           pp = n * (z*p1-p2) / (z*z-1.d0)
           z1 = z
           z  = z1 - p1 / pp
         if (abs(z-z1).gt.eps) go to 1
         x(i)     = xm - xl * z
         w(i)     = 2.d0 * xl / ((1.d0-z*z)*pp*pp)
         if (.not. sym_p) then
           x(n+1-i) = xm + xl * z
           w(n+1-i) = w(i)
         endif
c	  write(*,*) i,m
       enddo

       return
       end







       subroutine find_gaulag(x,w,n)
       
       integer n, maxit
       real*8 w(n), x(n)
       real*8 eps
       parameter(eps=3.d-14,maxit=30)
       
       integer i, its, j
       real ai  !, gammln
       real*8 p1, p2, p3, pp, z, z1
       
       do i = 1, n
         if ( i .eq. 1 ) then
           z = 3./(1.+2.4*n)
         elseif ( i .eq. 2 ) then
           z = z+15./(1.+2.5*n)
         else
           ai = i - 2
           z = z+(1.+2.55*ai)/(1.9*ai)*(z-x(i-2))
         endif
         do its = 1, maxit
           p1 = 1.d0
           p2 = 0.d0
           do j = 1, n
             p3 = p2
             p2 = p1
             p1 = ((2*j-1-z)*p2-(j-1)*p3)/j
           enddo
           pp = (n*p1-n*p2)/z
           z1 = z
           z = z1-p1/pp
           if (abs(z-z1) .le. eps) go to 1
         enddo
         pause 'too many iterations in gaulag'
 1	  x(i) = z
c 	  w(i) = -exp(gammln(alf+n)-gammln(float(n)))/(pp*n*p2)
 	  w(i) = -1.d0 / ( pp*n*p2 )
       enddo
       
       return
       end
       
       
       

       
       
       
       
       
       function k1(x)
c return the value of the modified bessel function of the second kind 
c of order one: K_1(x). these expresions are from AJM Hitchcock, 
c Rev. Math. Tables Aids Comput. v11 p86 (1950) and YL Luke,
c Integrals of Bessel Functions (NY: McGraw-Hill) (1962).
       implicit real*8 (a-h,k,o-z)

       if (x.gt.1.0d0) then
         t=1.0d0/x
         k1=exp(-x)*sqrt(t)*((((((((((((
     ^	    -0.0108241775d0*t+0.0788000118d0)*t-0.2581303765d0)*t
     ^	    +0.5050238576d0)*t-0.6632295430d0)*t+0.6283380681d0)*t
     ^	    -0.4594342117d0)*t+0.2847618149d0)*t-0.1736431637d0)*t
     ^	    +0.1280426636d0)*t-0.1468582957d0)*t+0.4699927013d0)*t
     ^	    +1.2533141373d0)
       else
         t=x*0.5d0
         t=t*t
         tt=x/3.75d0
         tt=tt*tt
         k1=
     ^	    ((((((-0.00004686d0*t-0.00110404d0)*t-0.01919402d0)*t
     ^	    -0.18156897d0)*t-0.67278579d0)*t+0.15443144d0)*t
     ^	    +1.d0)/x
     ^	    +((((((
     ^	    +0.00032411d0*tt+0.00301532d0)*tt+0.02658733d0)*tt
     ^	    +0.15084934d0)*tt+0.51498869d0)*tt+0.87890594d0)*tt
     ^	    +0.5d0)*x*log(x*0.5d0)
       endif
       return
       end

       


       function ki1(x)
c return the value of the integral from x to infinity of the 
c modified bessel function of the second kind of order zero. 
c these expresions are from AJM Hitchcock, 
c Rev. Math. Tables Aids Comput. v11 p86 (1950) and YL Luke,
c Integrals of Bessel Functions (NY: McGraw-Hill) (1962).

       implicit real*8 (a-h,k,o-z)
c	parameter (piby2=acos(0.0d0))

       if (x.gt.1.0d0) then
         t=1.d0/x
         ki1=exp(-x)*sqrt(t)*(((((((((((((-2.1974824449d0*t
     ^	    +16.7887658787d0)*t-57.8890096515d0)*t+119.2265927008d0)*t
     ^	    -163.7482102377d0)*t+158.8327470627d0)*t-112.8072652384d0)*t
     ^	    +60.4628882856d0)*t-25.4391219592d0)*t+9.0256045356d0)*t
     ^	    -3.0905443850d0)*t+1.2573312033d0)*t-0.7832344963d0)*t
     ^	    +1.2533139163d0)
       elseif (x.gt.0.0d0) then
         t=x*0.5d0
         tt=t*t
         ki1= dacos(0.0d0)
     ^	    -(((((((
     ^	     0.00000116d0*tt+0.00002069d0)*tt+0.00062664d0)*tt
     ^	    +0.01110118d0)*tt+0.11227902d0)*tt+0.50407836d0)*tt
     ^	    +0.84556868d0)*t)
     ^	    +(((((((
     ^	     0.000000319d0*tt+0.000012590d0)*tt+0.000385833d0)*tt
     ^	    +0.007936494d0)*tt+0.100000003d0)*tt+0.666666667d0)*tt
     ^	    +2.0d0)*t)*log(t)
       else
         ki1=dacos(0.0d0)
       endif
       return
       end



