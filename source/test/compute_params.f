c this program produces a simple table for Pcool, av and k_bar at different
c pressures for Argon. NOT needed for main pfr program.


	implicit real*8 (a-h,o-z)
	real*4 a

	bk = 1.380662d-23
	gu = 3.
	gl = 1.
	vth = sqrt(2.d0 * bk * 300.d0 / 
     ^		(39.948 * 1.660565d-27))
	tv = 8.6d-9
	wl = 106.66d-9
	pi = acos(-1.d0)
	
	write(*,"(8x,'Torr',19x,'Pcoll',20x,'av',19x,'k_bar')")
	do i = -2, 1
	  do j = 0, 10
	    a = 0.10 * j + i
	    p = 10.0**a
	    den = p *133.3333d0 / (bk * 300.d0)	!density (1/m^3)
	    
	    bar_k = gu * wl * wl * wl * den /
     ^	    		(8.d0 * pi * gl * tv * vth)
	    av    = (1.d0 + gu * wl * wl * wl * den / (gl*65.71d0))
     ^	    		* wl / (4.d0 * pi * tv * vth)
	    pc    = 1.d0 / (1.d0 + gl*65.71d0 / (gu*wl*wl*wl*den))
	    
	    write(*,"(1p,4e23.15)") p,pc,av,bar_k/100.d0
	  enddo
	enddo

	pause

	end
