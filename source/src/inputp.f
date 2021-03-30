       subroutine StartUp

       implicit none
       include 'rad.inc'
       integer*4 i,j
       logical log_p
       character*80 str1
       real*8 vt_v, vt_v0, x_mx, x_mn, frac0, k000 !, frac1
       logical k_lc_p
       real k00
       logical dk_p
       integer freq_grid
       
       real ave_v, voigt, findfreq_v, rint_v, rint_dopp, rint_lrtz
       real dk(2,kid), freq_min, freq_max, k_min, k_max, freq_mid
       real k_mid, x1, x2, x3
       
       namelist /rad/ zmaxdp,k0,tvac,
     ^	  izmax,nave,num_ndt,print_ndt,ndt,
     ^	  decay_p,prod_p,slab_p,cyl_p,sph_p,
     ^	  prodfile,filedat,pathname,freq_grid,
     ^	  ikmax,ikbmax,av,pcoll,esc_clear_p,
     ^	  jw_p, jw_dopp_p, read_p, save_p, oldfile, newfile,
     ^    pure_dopp_p, cfr_dopp_p, cfr_lrtz_p, cfr_vgt_p, k_lc_p,
     ^    exact_coh_p, exact_incoh_p, n_phi_gauss, n_r_gauss,
     ^	  freq_min, freq_max, k_min, k_max, freq_mid, k_mid
     
     
     	zmaxdp	= 1.d0
     	k0	= 1.d5
     	tvac	= 1.d0
     	izmax	= zid
     	nave	= 2

     	num_ndt	= 10000
     	print_ndt = 10
     	ndt	= 100

     	decay_p = .true.
     	prod_p	= .false.

     	slab_p	= .false.
     	cyl_p	= .true.
     	sph_p	= .false.

       n_phi_gauss = 300
       n_r_gauss   = -1

     	pathname = ''
     	prodfile ='none'
     	filedat	= 'rad.out'
     	
     	ikmax	  = kid
     	ikbmax	  = -1    !0.3 * ikmax
     	freq_grid = 1
     	freq_min  = -1.
     	freq_mid  = -1.
     	freq_max  = -1.
     	k_min     = -1.
     	k_mid     = -1.
     	k_max     = -1.
     	
     	av	= 0.1
     	Pcoll	= 0.5
       esc_clear_p = .true.

       jw_p	= .false.
       jw_dopp_p = .false.

       pure_dopp_p = .false.

       cfr_dopp_p  = .false.
       cfr_lrtz_p  = .false.
       cfr_vgt_p   = .false.

       exact_coh_p = .false.
       exact_incoh_p = .false.

       
       k_lc_p  = .false.
       
       read_p	= .false.
       save_p	= .true.
       oldfile	= 'oldprop.dat'
       newfile = 'newprop.dat'
       
       nout = 0
       do j = 1, kid
         do i = 1, zid
           prra(i,j) = 0.
         enddo
       enddo
       
       inquire(file='rad.in',exist=log_p)
       if (log_p) then
         open(14,file='rad.in')
          read(14,rad)
         close(14)
       else
         write(*,"(' ERROR- no input file: rad.in')")
         write(*,"('  run aborting...')")
c        pause
         stop
       endif

       n_phi_gauss = min( 400, max( 50, n_phi_gauss ))
       if (n_r_gauss.le.0) n_r_gauss = int( nave / 2 )
       n_r_gauss = max( 1, n_r_gauss )
       
       if (read_p) call read_old
       
       izmax = min(izmax, zid)

       if (freq_grid.eq.0) then
         if (ikbmax.lt.2.or.ikbmax.gt.ikmax-2) ikbmax=0.6*ikmax
         if (ikmax .gt. kid) then
           ikbmax = real(ikbmax)*kid/ikmax
           ikmax = min(ikmax, kid)
           ikbmax = max(2,min(ikmax-2, ikbmax))
         endif
         dk_p = .false.
       else
         dk_p = .true.
         if (freq_grid .lt. 0. or. freq_grid .gt. 3) freq_grid = 1
       endif
       
       filedat = PathName(1:index(pathname,' ')-1)//
     ^		filedat(1:index(filedat,' ')-1)
       oldfile = PathName(1:index(pathname,' ')-1)//oldfile
       newfile = PathName(1:index(pathname,' ')-1)//newfile
     
       if (prod_p) then
         prodfile = PathName(1:index(pathname,' ')-1)//
     ^		prodfile(1:index(prodfile,' ')-1)
         inquire(file=prodfile,exist=prod_p)
         if (prod_p) then
           open(14,file=prodfile,status='old')
 57	      read(14,"(a)",end=99,err=99) str1
 	      if (str1(1:4) .ne. '>>>>') go to 57
 	      read(14,*,end=99,err=99) (prra(i,1),i=1,izmax)
           close(14)
         endif
       endif
       
       decay_p = (.not. prod_p)
       
       if (cyl_p) then
         slab_p = .false.
         sph_p  = .false.
       elseif (slab_p) then
         cyl_p  = .false.
         sph_p  = .false.
       elseif (sph_p) then
         cyl_p  = .false.
         slab_p = .false.
       else
         cyl_p  = .true.
         slab_p = .false.
         sph_p  = .false.
       endif
       
       if (jw_p) then
         pure_dopp_p = .false.
         cfr_dopp_p  = .false.
         cfr_lrtz_p  = .false.
         cfr_vgt_p   = .false.
         exact_coh_p = .false.
         exact_incoh_p = .false.
       elseif (exact_coh_p) then
         jw_p        = .false.
         pure_dopp_p = .false.
         cfr_dopp_p  = .false.
         cfr_lrtz_p  = .false.
         cfr_vgt_p   = .false.
         exact_incoh_p = .false.
       elseif (exact_incoh_p) then
         jw_p        = .false.
         pure_dopp_p = .false.
         cfr_dopp_p  = .false.
         cfr_lrtz_p  = .false.
         cfr_vgt_p   = .false.
         exact_coh_p = .false.
       elseif (pure_dopp_p) then
         jw_p        = .false.
         cfr_dopp_p  = .false.
         cfr_lrtz_p  = .false.
         cfr_vgt_p   = .false.
         exact_coh_p = .false.
         exact_incoh_p = .false.
       elseif (cfr_dopp_p) then
         jw_p        = .false.
         pure_dopp_p = .false.
         cfr_lrtz_p  = .false.
         cfr_vgt_p   = .false.
         exact_coh_p = .false.
         exact_incoh_p = .false.
       elseif (cfr_lrtz_p) then
         jw_p        = .false.
         pure_dopp_p = .false.
         cfr_dopp_p  = .false.
         cfr_vgt_p   = .false.
         exact_coh_p = .false.
         exact_incoh_p = .false.
       elseif (cfr_vgt_p) then
         jw_p        = .false.
         pure_dopp_p = .false.
         cfr_dopp_p  = .false.
         cfr_lrtz_p  = .false.
         exact_coh_p = .false.
         exact_incoh_p = .false.
       else
         jw_p        = .true.
         pure_dopp_p = .false.
         cfr_dopp_p  = .false.
         cfr_lrtz_p  = .false.
         cfr_vgt_p   = .false.
         exact_coh_p = .false.
         exact_incoh_p = .false.
       endif
       
       if (.not.read_p .and. k_lc_p) then
         if (cfr_lrtz_p) then
           k0 = k0 * av * pi
         elseif (pure_dopp_p .or. cfr_dopp_p) then
           k0 = k0 * sqrt(pi)
         else
           k0 = k0 / voigt(0., real(av))
         endif
       endif

       if (cfr_lrtz_p) then
         k00 = k0 / ( av * pi )
       elseif (pure_dopp_p .or. cfr_dopp_p) then
         k00 = k0 / sqrt(pi)
       else
         k00 = k0 * voigt(0., real(av))
       endif

       zmax = zmaxdp
       dzdp = zmaxdp / izmax
       dz   = dzdp
       
       if (cyl_p) then
         write(*,"('  cylindrical geometry, radius =',
     ^	  	1p,e13.5,' (cm)')") zmaxdp
         write(*,"('  [',i3,' azimuthal and',i3,' radial Gaussian',
     ^		' quadrature point(s)]')") n_phi_gauss,n_r_gauss
         do i = 1, izmax
           prra(i,1) = prra(i,1) * pi * dzdp * dzdp * (2.d0*i-1.d0)*tvac
         enddo
       elseif (slab_p) then
         write(*,"('  slab geometry, gap =',1p,e13.5,' (cm)')") zmaxdp
         do i = 1, izmax
           prra(i,1) = prra(i,1) * dzdp * tvac
         enddo
       else
         write(*,"('  spherical geometry, radius =',
     ^	  	1p,e13.5,' (cm)')") zmaxdp
         do i = 1, izmax
           prra(i,1) = prra(i,1) * tvac *
     ^	    	4.d0/3.d0*pi*dzdp*dzdp*dzdp*(i*i*i - (i-1)**3)
         enddo
       endif
       
       do j = 1, ikmax
         do i=1,izmax+1
           nra(i,j,1)=0.0
           nra(i,j,2)=0.0
         enddo
         nesc(0,j) = 0.d0
         nesc(1,j) = 0.d0
       enddo
       ci=1

     	if (decay_p) then
     	  write(*,"('  decay simulation- no external prodution')")
       else
     	  write(*,"('  steady state simulation- ',
     ^     	  'external prodution via data file: ',a)") 
     ^	          prodfile(1:max(1,index(prodfile,' ')-1))
       endif
       write(*,"('  # of spatial cells = ',
     ^		i3)") izmax
       write(*,"('  # of frequency cells = ',
     ^		i3)") ikmax
       write(*,"('  integrated absorption coefficient, k-bar =',
     ^		1p,e13.5,' (1/cm)')") k0
       write(*,"('  line center absorption coefficient, k_0 =',
     ^		1p,e13.5,' (1/cm)')") k00
       write(*,"('  Voigt paramenter =',1p,e13.5)") av
       write(*,"('  probability of dephasing collision =',
     ^		1p,e13.5)") Pcoll
        if (jw_p) then
          if (jw_dopp_p) then
            write(*,"('  Jefferies-White frequency redistribution ',
     ^           'w/ shifted Doppler')")
          else
            write(*,"('  Jefferies-White frequency redistribution')")
          endif
        elseif (pure_dopp_p) then
          write(*,"('  pure Doppler redistribution [erfc(x>)]')")
        elseif (cfr_dopp_p) then
          write(*,"('  CFR with Doppler lineshape')")
        elseif (cfr_lrtz_p) then
          write(*,"('  CFR with Lorentz lineshape')")
        elseif (cfr_vgt_p) then
          write(*,"('  CFR with Voigt lineshape')")
        elseif (exact_incoh_p) then
          write(*,"('  angle averaged exact frequency redistribution ',
     ^          	'(incoherent scattering)')")
        else
          write(*,"('  angle averaged exact frequency redistribution ',
     ^          	'(coherent scattering)')")
        endif
       write(*,"('  vacuum lifetime, Tvac =',
     ^		1p,e13.5,' (s)')") tvac
       write(*,"('  spatial distribution printed every ',
     ^		i6,' timesteps')") ndt * print_ndt
       write(*,"('  maximum simulation time is ',
     ^		i8,' timesteps')") ndt * num_ndt
       write(*,"('  initial cell is averaged ',
     ^		i3,' times')") nave
        if (esc_clear_p) then
          write(*,"('  exit spectrum cleared after each print out')")
        else
          write(*,"('  exit spectrum is running total')")
        endif
     
        vt_v0 = av * tan( asin(1.d0) * (1.d0 - 1.d-5) )
        if (dk_p) vt_v0 = k0 *
     ^    av / acos(-1.d0) / ( vt_v0*vt_v0 + av*av)
     
       x1 = min( freq_min, freq_mid, freq_max )
       x3 = max( freq_min, freq_mid, freq_max )
       x2 = min(           freq_mid, freq_max )
       if ( x1 .eq. x2 ) 
     ^		x2 = min( max( freq_mid, freq_max), freq_min )
       if ( x1 .ge. 0. ) then
         freq_min = x1
         freq_mid = x2
         freq_max = x3
       elseif ( x2 .ge. 0.) then
         if ( freq_min .ge. 0. ) freq_min = x2
         if ( freq_mid .ge. 0. ) then
           freq_mid = x2
           if ( freq_min .eq. x2 ) freq_mid = x3
         endif
         if ( freq_max .ge. 0. ) freq_max = x3
       endif
       
       x1 = min( k_min, k_mid, k_max )
       x3 = max( k_min, k_mid, k_max )
       x2 = min(        k_mid, k_max )
       if ( x1 .eq. x2 ) x2 = min( max( k_mid, k_max), k_min )
       if ( x1 .ge. 0. ) then
         k_min = x1
         k_mid = x2
         k_max = x3
       elseif ( x2 .ge. 0.) then
         if ( k_min .ge. 0. ) k_min = x2
         if ( k_mid .ge. 0. ) then
           k_mid = x2
           if ( k_min .eq. x2 ) k_mid = x3
         endif
         if ( k_max .ge. 0. ) k_max = x3
       endif
       
            
        if (cfr_lrtz_p) then		!lorentzian line shape

c frequency where 99% travel 2 L...
         vt_v = -log(0.99d0) / (2.d0*zmaxdp) / k0
         x_mx = min( vt_v0, 
     ^	  	sqrt( abs( av/(pi*vt_v) - av*av )) )
         if (dk_p) x_mx = max( vt_v0, vt_v * k0 )
         if (freq_grid .eq. 2 .and. k_min .ge. 0. ) x_mx = k_min
         if (freq_grid .eq. 3 .and. freq_max .ge. 0. )
     ^      x_mx = k0 * av / (pi * (freq_max*freq_max + av*av))
       
c frequency where k0 Lv(x) L = 1...
         vt_v = 1.d0 / (k0 * zmaxdp)
         x_mn = sqrt( abs( av/(pi*vt_v) - av*av ))
         if (dk_p) x_mn = vt_v * k0
         if (freq_grid .eq. 2 .and. k_mid .ge. 0. ) x_mn = k_mid
         if (freq_grid .eq. 3 .and. freq_mid .ge. 0. )
     ^      x_mn = k0 * av / (pi * (freq_mid*freq_mid + av*av))

         if (freq_min .ge. 0.) 
     ^	    freq_min = k0 * av / (pi * (freq_min*freq_min + av*av))
       
       elseif (pure_dopp_p .or. cfr_dopp_p) then  !doppler lineshape (gaussian)

c frequency where 99% travel 2 L...
         vt_v = -log(0.99d0) / (2.d0*zmaxdp) / k0
         x_mx = min( vt_v0, 
     ^	  	sqrt( abs( log( vt_v * sqrt(pi) ) )) )
         if (dk_p) x_mx = max( vt_v0, vt_v * k0 )
         if (freq_grid .eq. 2 .and. k_min .ge. 0. ) x_mx = k_min
         if (freq_grid .eq. 3 .and. freq_max .ge. 0. )
     ^      x_mx = k0 * exp(-freq_max*freq_max) / sqrt(pi)
       
c frequency where k0 Lv(x) L = 1...
         vt_v = 1.d0 / (k0 * zmaxdp)
         x_mn = sqrt( abs( log( vt_v * sqrt(pi) ) ))
         if (dk_p) x_mn = vt_v  * k0
         if (freq_grid .eq. 2 .and. k_mid .ge. 0. ) x_mn = k_mid
         if (freq_grid .eq. 3 .and. freq_mid .ge. 0. )
     ^      x_mn = k0 * exp(-freq_mid*freq_mid) / sqrt(pi)

         if (freq_min .ge. 0.) 
     ^	    freq_min = k0 * exp(-freq_min*freq_min) / sqrt(pi)
       
       else 
     
c frequency where 99% travel 2 L...
         vt_v = -log(0.99d0) / (2.d0*zmaxdp) / k0
         x_mx = min( real(vt_v0), findfreq_v(vt_v,av) )
         if (dk_p) x_mx = max( vt_v0, vt_v  * k0 )
         if (freq_grid .eq. 2 .and. k_min .ge. 0. ) x_mx = k_min
         if (freq_grid .eq. 3 .and. freq_max .ge. 0. )
     ^      x_mx = k0 * voigt(real(freq_max),real(av))
       
c frequency where k0 Lv(x) L = 1...
         x_mn = findfreq_v(1.d0 / (k0 * zmaxdp), av)
         if (dk_p) x_mn = 1.d0 / zmaxdp
         if (freq_grid .eq. 2 .and. k_mid .ge. 0. ) x_mn = k_mid
         if (freq_grid .eq. 3 .and. freq_mid .ge. 0. )
     ^      x_mn = k0 * voigt(real(freq_mid),real(av))

         if (freq_min .ge. 0.) 
     ^	    freq_min = k0 * voigt(real(freq_min),real(av))
       
       endif ! voigt lineshape
       
c find ikbmax such that the fractional change of k(x) is  
c about equal around x', where k(x')R ~ 1 (or user specified or half way between extremes)
       if (dk_p) then
         k000   = min( dble(k00), -2.d0*log(1.d-6)/dzdp)
         if (freq_grid .eq. 2 .and. k_max .ge. 0. ) k000 = k_max
         if (freq_grid .eq. 3 .and. freq_min .ge. 0.) k000 = freq_min
         if ( freq_grid .eq. 1 .or. 
     ^	      (freq_grid .eq. 2 .and. k_mid .lt. 0.) .or.
     ^	      (freq_grid .eq. 3 .and. freq_mid .lt. 0. ) )
     ^	         x_mn   = exp( 0.5d0*(log(k000) + log(x_mx)) )
c	  write(*,*) k000,x_mn,x_mx
         x1 = min( x_mx, x_mn, k000 )
         x3 = max( x_mx, x_mn, k000 )
         x2 = min(       x_mn, k000 )
         if ( x1 .eq. x2 ) x2 = min( max( x_mn, k000), x_mx )
         x_mx = x1
         x_mn = x2
         k000 = x3	  
         vt_v   = log(k000 / x_mn)
         frac0  = log(x_mn / x_mx)
         if ( ikbmax .lt. 2 ) ikbmax =  nint(
     ^	    (ikmax * vt_v + frac0) / (vt_v + frac0) )
       endif
       ikbmax = max( 2, min( ikmax-2, ikbmax ))
     
       if (dk_p) then
c geometrically increasing away from x_mn to line center
         frac0 = (k000/x_mn)**(1.0/(ikbmax-1))
         dk(1,ikbmax) = x_mn
         do i = ikbmax-1, 1, -1
           dk(1,i) = dk(1,i+1) * frac0
         enddo
         dk(1,1) = k00

c geometrically increasing away from x_mn to wings
         frac0 = (x_mx/x_mn)**(1.0/(ikmax-ikbmax))
         do i = ikbmax+1, ikmax
           dk(1,i) = dk(1,i-1) * frac0
         enddo
         dk(2,ikmax) = dk(1,ikmax) * frac0

         df(1,1) = 0.d0

c	  if (freq_grid.eq.3) then
c	  if (.false.) then

c k(i+1) = k(i)^^x * f, from k(x*) to k(0)
c	    frac1 = max( freq_pow, 1.0001 )
c	    frac0 = (frac1 - 1.d0) / (frac1**(ikbmax-1)-1.d0)
c	    frac0 = (k000/x_mn**(frac1**(ikbmax-1)))**frac0
c	    dk(1,ikbmax) = x_mn
c	    do i = ikbmax-1, 1, -1
c	      dk(1,i) = frac0 * dk(1,i+1) ** frac1
c	    enddo

c	    frac0 = (frac1 - 1.d0) / (frac1**(ikmax-ikbmax)-1.d0)
c	    frac0 = (x_mx/x_mn**(frac1**(ikmax-ikbmax)))**frac0
c	    do i = ikbmax+1, ikmax
c	      dk(1,i) = frac0 * dk(1,i-1) ** frac1
c	    enddo
c	    dk(2,ikmax) = frac0 * dk(1,ikmax) ** frac1

c m*i+b = exp(-k(i)) ?
c	    frac1 =  exp(-k00)
c	    frac0 = (exp(-x_mx) - frac1) / real(ikmax)
c	    do i = 1, ikmax
c	      dk(1,i) = -log(frac0*(i-1.d0)+frac1)
c	    enddo
c	    dk(2,ikmax) = -log(frac0*ikmax+frac1)
c	    
c	    frac0 = 0.5d0/exp(-dk(1,ikmax-ikbmax))
c	    do i = 1, ikmax
c	      dk(1,i) = dk(1,i) * frac0
c	      dk(2,i) = dk(2,i) * frac0
c	    enddo
           	    
c	  endif

         do i = 1, ikmax-1
           dk(2,i)  = dk(1,i+1)
         enddo

       else
         
c geometrically increasing away from x_mn to x_mx
         frac0 = (x_mx/x_mn)**(1.0/(ikmax-ikbmax+1))
         df(1,ikbmax) = x_mn
         do i = ikbmax+1, ikmax
           df(1,i) = df(1,i-1) * frac0
         enddo
         df(2,ikmax) = df(1,ikmax) * frac0

c geometrically increasing away from x_mn to 0
         frac0 = (2.d0) ** (1.0/ikbmax)
         do i = ikbmax-1, 2, -1
           df(1,i) = df(1,i+1) * frac0
         enddo
         do i = ikbmax-1, 2, -1
           df(1,i) = 2.d0*df(1,ikbmax) - df(1,i)
         enddo

         df(1,1) = 0.d0

         do i = 1, ikmax-1
           df(2,i)  = df(1,i+1)
         enddo

       endif
       
       if (cfr_lrtz_p) then

         if (dk_p) then
           do i = 2, ikmax
             x_mx = dk(1,i) / k0
             df(1,i) = sqrt( abs( av/(pi*x_mx) - av*av ))
             df(2,i-1) = df(1,i)
           enddo
           x_mx = dk(2,ikmax) / k0
           df(2,ikmax) = sqrt( abs( av/(pi*x_mx) - av*av ))
         else
           do i = 1, ikmax
             dk(1,i) = k0 * av / (pi * (df(1,i)*df(1,i) + av*av))
             dk(2,i) = k0 * av / (pi * (df(2,i)*df(2,i) + av*av))
           enddo
         endif

         df(2,ikmax) = max(2.d0*df(2,ikmax), 1.d32)
         frac0 = 0.d0
         do i = 1, ikmax
           int_ls(i) = rint_lrtz( df(1,i), df(2,i), av)
           frac0     = frac0 + int_ls(i)
           fv(i)     = av / (2.d0 * pi) *
     ^	       abs( log ( (av*av + df(2,i)*df(2,i)) /
     ^                    (av*av + df(1,i)*df(1,i)) ))
     ^	    	   / int_ls(i) 
           kv(i)     = k0 * av / (pi * (fv(i)*fv(i) + av*av))
         enddo
       
       elseif (pure_dopp_p .or. cfr_dopp_p) then

         if (dk_p) then
           do i = 2, ikmax
             x_mx = dk(1,i) / k0
             df(1,i) = sqrt( abs( log( x_mx * sqrt(pi) ) ))
             df(2,i-1) = df(1,i)
           enddo
           x_mx = dk(2,ikmax) / k0
           df(2,ikmax) = sqrt( abs( log( x_mx * sqrt(pi) ) ))
         else
           do i = 1, ikmax
             dk(1,i) = k0 * exp(-df(1,i)*df(1,i)) / sqrt(pi)
             dk(2,i) = k0 * exp(-df(2,i)*df(2,i)) / sqrt(pi)
           enddo
         endif

         df(2,ikmax) = max(2.d0*df(2,ikmax), 1.d32)
         frac0 = 0.d0
         do i = 1, ikmax
           int_ls(i) = rint_dopp( df(1,i), df(2,i))
           frac0     = frac0 + int_ls(i)
           fv(i)     = 0.5d0 * 
     ^	       abs(exp(-df(1,i)*df(1,i)) - exp(-df(2,i)*df(2,i)))
     ^	    	   / (sqrt(pi) * int_ls(i) )
           kv(i)     = k0 * exp(-fv(i)*fv(i)) / sqrt(pi)
         enddo
         
       else		!Voigt lineshape

         if (dk_p) then
           do i = 2, ikmax
             x_mx = dk(1,i) / k0
             df(1,i) = findfreq_v(x_mx,av)
             df(2,i-1) = df(1,i)
           enddo
           x_mx = dk(2,ikmax) / k0
           df(2,ikmax) = findfreq_v(x_mx,av)
         else
           do i = 1, ikmax
             dk(1,i) = k0 * voigt(real(df(2,i)),real(av))
             dk(2,i) = k0 * voigt(real(df(2,i)),real(av))
           enddo
         endif

c now find the average frequency (ie averaged over Voigt) in each cell...
         do i = 1, ikmax
           fv(i) = ave_v(df(1,i),df(2,i),av)
           kv(i) = k0 * voigt(real(fv(i)),real(av))
         enddo
       
         int_ls(ikmax) = max(0.d0, 0.5d0 - 
     ^		atan(df(1,ikmax)/av)/acos(-1.d0))
         frac0 = int_ls(ikmax)
         do i = 1, ikmax-1
           int_ls(i) = rint_v(df(1,i),df(2,i),400,av)
           frac0 = frac0 + int_ls(i)
         enddo

       endif ! use voigt lineshape...
       
       frac0 = 1.d0 / frac0
       do i = 1, ikmax
         int_ls(i) = int_ls(i) * frac0
       enddo

       open(22,file=filedat(1:index(filedat,' ')-1)//'.freq')
       do i = 1, ikmax
         str1 = ' '
         if (i .eq. ikbmax) str1 ='*'
c	  write(*,"(a,i3,1p,7e13.6)") str1(1:1),i,kv(i),fv(i),
c     ^	    int_ls(i),
c     ^	    df(1,i),df(2,i),(dk(1,i)-dk(2,i))/kv(i),
c     ^	    (df(2,i)-df(1,i))/fv(i)
         write(22,"(a,i3,1p,7e13.6)") str1(1:1),i,kv(i),fv(i),
     ^	    int_ls(i),
     ^	    df(1,i),df(2,i),(dk(1,i)-dk(2,i))/kv(i),
     ^	    (df(2,i)-df(1,i))/fv(i)
       enddo
       write(22,"(2x,'#',4x,'<k>',10x,'<f>',10x,'int L',8x,'f<',
     ^	    11x,'f>',11x,'(f>-f<)/<f>',2x,'(k>-k<)/<k>')")
       close(22)
c      pause
       
       do i = ikmax, 1, -1
         do j = 1, izmax
           prra(j,i) = prra(j,1) * int_ls(i) 
         enddo
       enddo

       return

c done!, go compute
 99	continue
       write(*,"(' ERROR- fatal error while reading ',a)") prodfile
       write(*,"('  run aborting...')")
c      pause
       stop
       
       end





       subroutine read_old
       
       implicit none
       include 'rad.inc'
       integer i, j, k

       
       inquire(file=oldfile,exist=read_p)
       if (read_p) then
          open(14,file=newfile,status='old')
           read(14,*,err=99,end=99) zmaxdp,k0,tvac,izmax,nave,
     ^	      slab_p,cyl_p,sph_p,ikmax,ikbmax,av,pcoll,jw_p,jw_dopp_p,
     ^	      pure_dopp_p, cfr_dopp_p, cfr_lrtz_p, cfr_vgt_p,
     ^        exact_coh_p, exact_incoh_p, n_phi_gauss, n_r_gauss
           do i = 1, ikmax
             read(14,*,err=99,end=99) (pac(k,i),k=1,ikmax)
           enddo
           do k = 1, ikmax
             do i = 1, izmax
               read(14,*,err=99,end=99) (pa(j,i,k),j=0,izmax+1)
             enddo
           enddo
          close(14)
          write(*,"('    [old distribution read]')")
       endif

       return
       
 99	continue
       close(14)
       read_p = .false.
       
       return
       end
       






       subroutine save_new
       implicit none
       include 'rad.inc'
       integer i, j, k

       
       open(14,file=newfile,status='unknown')
        write(14,*) zmaxdp,k0,tvac,izmax,nave,
     ^	      slab_p,cyl_p,sph_p,ikmax,ikbmax,av,pcoll,jw_p,jw_dopp_p,
     ^	      pure_dopp_p, cfr_dopp_p, cfr_lrtz_p, cfr_vgt_p,
     ^        exact_coh_p, exact_incoh_p, n_phi_gauss, n_r_gauss
        do i = 1, ikmax
          write(14,"(1p,300e14.7)") (pac(k,i),k=1,ikmax)
        enddo
        do k = 1, ikmax
          do i = 1, izmax
            write(14,"(1p,300e14.7)") (pa(j,i,k),j=0,izmax+1)
          enddo
        enddo
       close(14)
       return
       
       end
       
       
       
