       subroutine WriteInt(i0)

       implicit none
       include 'rad.inc'
       integer*4 j,i,i0
       character*90 namefn,FileName

       real*4 z, vol, den(3,zid+kid), mxd
       real*8 xx
       real  mxf0, mxf1, mxf2
       real*8 xf0,  xf1

       filename=namefn(filedat,nout)
       open(9,file=filename,status='unknown')
       write(9,"('output after ',i9,' time steps')") i0
       write(9,"('density of excited resonant atoms:')")
       write(9,"(6x,'x',10x,'density',6x,'number',8x,'ratio')")
       
       mxd = 0.
       xf0 = 0.d0
       do j=1,izmax
         if (cyl_p) then
           vol=pi*(j*j-(j-1)*(j-1))*dzdp*dzdp
         elseif (slab_p) then
           vol = dzdp
         else
           vol = 4.d0/3.d0*pi*dzdp*dzdp*dzdp*(j*j*j - (j-1)**3)
         endif
         z = max(0.d0,min(zmaxdp,(j-0.5d0)*dzdp))
         
         xx = 0.d0
         do i = ikmax, 1, -1
           xx = xx + nra(j,i,ci)
         enddo
         xf0      = xf0 + xx
         den(1,j) = z
         den(2,j) = xx/vol
         den(3,j) = xx
         mxd = max( mxd, den(2,j) )
       enddo
       
       if (mxd .ne. 0.0) mxd = 1.0 / mxd
       do j = 1, izmax
         write(9,"(1p,4e13.5)") (den(i,j),i=1,3),den(2,j)*mxd
       enddo
       write(9,"(18x,'total #:',1p,e13.5)") xf0

       xf0 = 0.d0
       xf1 = 0.d0
       mxf0 = 0.
       mxf1 = 0.
       mxf2 = 0.
       do i = 1, ikmax
         xf0 = xf0 + nesc(0,i)
         xf1 = xf1 + nesc(1,i)
         den(1,i) = nesc(0,i) / (df(2,i)-df(1,i))
         den(2,i) = nesc(1,i) / (df(2,i)-df(1,i))
         den(3,i) = (nesc(0,i)+nesc(1,i)) / (df(2,i)-df(1,i))
         mxf0 = max( mxf0, den(1,i) )
         mxf1 = max( mxf1, den(2,i) )
         mxf2 = max( mxf2, den(3,i) )
       enddo
       if (mxf0.ne.0.0) mxf0 = 1.0 / mxf0
       if (mxf1.ne.0.0) mxf1 = 1.0 / mxf1
       if (mxf2.ne.0.0) mxf2 = 1.0 / mxf2

       write(9,*)
       write(9,"(' spectrum escaping:')")
       if (slab_p) then
         write(9,"(5x,'freq',10x,'density',7x,'number',9x,'ratio',8x,
     ^      'density',7x,'number',9x,'ratio')")
         write(9,"(18x,'  (z=0)',8x,'(z=0) ',9x,'(z=0)',7x,
     ^      '  (z=L)',8x,'(z=L) ',9x,'(z=L)')")
         do i = 1, ikmax
                                                                        
           write(9,"(1p,7e14.5)") fv(i),den(1,i),nesc(0,i),
     ^      	den(1,i)*mxf0,den(2,i),nesc(1,i),den(2,i)*mxf1
         enddo
         write(9,"(20x,'total #:',1p,e14.5,28x,e14.5)") xf0,xf1
       else
         write(9,"(4x,'freq',9x,'density',6x,'number',9x,'ratio')")
         do i = 1, ikmax
           write(9,"(1p,4e13.5)") fv(i),den(3,i),
     ^  	 nesc(0,i)+nesc(1,i),den(3,i) * mxf2
         enddo
         write(9,"(18x,'total #:',1p,e13.5)") xf0+xf1
       endif
       
       if (esc_clear_p) then
         do i = 1, ikmax
           nesc(0,i) = 0.d0
           nesc(1,i) = 0.d0
         enddo
       endif
       close(9)
       nout = nout + 1
       return
       end










       function ave_v(y18,y28,a8)
c function returns the average of freq over the voigt profile from y18 to y28,
c integral of voigt lineshape is approximated via an extended Simpson's Rule
c (see Numerical Recipes for example, W.H.Press, B.P. Flannery...)
c 401 is the number of mesh points...

       real*8 y18,y28,a8
       
       
       sqrtpi = 1.d0 / sqrt(3.1415926535897932385d0)

       a = a8
       x1=min(abs(y18),abs(y28))
       x2=max(abs(y18),abs(y28))
c	n=n1+1-mod(n1,2)		!want n to be odd...
       n=401
       dx=(x2-x1)/n
       s=voigt(x2,a)/3.
       s1=s*x2
       do i=2,n-1
         x2=x2-dx
         if (mod(i,2).eq.0) then
           x=4./3.
         else
           x=2./3.
         endif
         x=x*voigt(x2,a)
         s=s+x
         s1=s1+x2*x
       enddo
       x=voigt(x1,a)/3.
       s=s+x
       s1=s1+x1*x
       ave_v=s1/s
       return
       end




       function doppler(x)
c function returns value of doppler lineshape at frequency x.
c
c	parameter (sqrtpi=1.0/sqrt(3.1415926535897932385))
       real*8 x

       doppler = exp(-x*x)/sqrt(3.1415926535897932385d0)
       return
       end





       function rint_dopp(y1,y2)
c function returns the integral of the doppler profile from y1 to y2,
c integral of doppler lineshape from 0 to x is 0.5*erf(x)...  
c use approx of erf given in Abramowitz and Stegun (handbook of mathematicla 
c functions, chapter 7, eq 7.1.26, pg 299, error < 1.5e-7

       real*8 y1,y2

       z1 = min(abs(y1), abs(y2))
       z2 = max(abs(y1), abs(y2))
       x1=1./(1.+0.3275911*z1)
       x2=1./(1.+0.3275911*z2)
       s1=((((1.061405429*x1-1.453152027)*x1+1.421413741)*x1
     ^    -0.284496736)*x1+0.254829592)*x1
       s2=((((1.061405429*x2-1.453152027)*x2+1.421413741)*x2
     ^    -0.284496736)*x2+0.254829592)*x2
       rint_dopp=0.5*abs(s1*exp(-z1*z1)-s2*exp(-z2*z2))
       return
       end





       function erfc(y1)
c function returns the value of the complimentary error function, erfc(y), 
c where erfc(y) = integral from y to inifinity of 2 exp(x^2) / sqrt(pi)
c erfc(y) = 1 - erf(y), where erf(y) is the error function. 	
c use approx of erf given in Abramowitz and Stegun (handbook of mathematicla 
c functions, chapter 7, eq 7.1.26, pg 299, error < 1.5e-7
       
       real*8 y1

       z1 = real(y1)
       x1 = 1./(1.+0.3275911*z1)
       s1 = ((((1.061405429*x1-1.453152027)*x1+1.421413741)*x1
     ^    -0.284496736)*x1+0.254829592)*x1
       erfc = s1*exp(-z1*z1)
       
       return
       end



       function rint_lrtz(y1,y2,y)
c function returns the integral of the lorentz profle from y1 to y2,
c integral of lorentz profile from 0 to x is arctan(x/a)/pi, w/ parameter a.

       real*8 y1,y2,y
       parameter (vpi=1./3.1415926535897932385)

       x1=min(abs(y1),abs(y2)) / y
       x2=max(abs(y1),abs(y2)) / y
       rint_lrtz=(atan(x2)-atan(x1))*vpi
       return
       end











       function voigt(x,y)
c computes the complex probability function w(z)=exp(-z*z)*erfc(-iz)
c in the upper half plane.  maximum relative error of both real and 
c imaginary parts is < 10^-4.  from J. Humlicek, JQSRT, vol 27, no. 4, p 437.
c voigt profile = real(w(z))/sqrt(pi), z= freq + i av, av= voigt parameter...
c x8=freq, y8= voigt parameter
       
       complex z,t,u
       real*4 x,y, voigt
       real*8 sqrtpi
c	parameter (sqrtpi=1.0/sqrt(3.1415926535897932385))

       sqrtpi = 1.d0 / sqrt(3.1415926535897932385d0)
       z = cmplx(x,y)
       t=cmplx(y,-x)
       s=abs(x)+y
 1	if (s.lt.5.5)goto2

       u=t*t
       voigt=real(t*(1.410474+u*.5641896)/(0.75+u*(3.+u)))*sqrtpi
       return
 2	if (y.lt..195*abs(x)-.176)goto 3

       voigt=real(
     ^  (16.4955+t*(20.20933+t*(11.96482+t*(3.778987+t*.5642236))))/
     ^  (16.4955+t*(38.82363+t*(39.27121+t*(21.69274+t*
     ^  (6.699398+t))))))*sqrtpi
       return
 
 3	u=t*t
       voigt=real(
     ^  cexp(u)-t*(36183.31-u*(3321.9905-u*(1540.787-u*(219.0313-u*
     ^  (35.76683-u*(1.320522-u*.56419))))))/(32066.6-u*(24322.84-u*
     ^  (9022.228-u*(2186.181-u*(364.2191-u*(61.57037-u*
     ^  (1.84139-u))))))))*sqrtpi

       return
       end








       

       
       function findfreq_v(vt8,a8)
c give the value of Lv (=vt8) and the Voigt parameter (a8), 
c find the frequeny x which gives this value!
       
       real*8 vt8,a8
       real x1, x2, x0, voigt
       parameter(eps=1.e-5)
       
       y = vt8
       a = a8
       
c	write(*,*) 'findfreq_v ',y,a
       
       x1 = 0.
       x2 = 5.0
       
       xv = voigt(x1,a)
       if (y.ge.xv) then
         findfreq_v = 0.
         return
       endif
       
c	write(*,*) x1,x2,xv
       
c	xv = voigt(x2,a)
       do while (voigt(x2,a) .gt. y)
         x1 = x2
         x2 = x2 * 2.0
c	  xv = voigt(x2,a)
c	  write(*,*) x1,x2,xv
       enddo

c	write(*,*) x1,x2,xv
       
       x0=0.5*(x1+x2)
       xv=voigt(x0,a)
       do while (abs(xv-y)/y .gt.eps)
         if (xv.gt.y) then
           x1=x0
         else
           x2=x0
         endif
         x0=0.5*(x1+x2)
         xv=voigt(x0,a)
         
c	  write(*,*) x0, xv
       enddo
       
c	write(*,*) x0
       findfreq_v=x0
       return
       	
       end
       
       
       

       function rint_v(y18,y28,n1,a8)
c function returns the integral of the voigt profile from y1 to y2,
c integral of voigt lineshape is approximated via an extended Simpson's Rule
c (see Numerical Recipes for example, W.H.Press, B.P. Flannery...)
c n1 is the number of mesh points...

       real*8 y18,y28,a8
c	parameter (sqrtpi=1.0/sqrt(3.1415926535897932385))

       sqrtpi = 1.d0 / sqrt(3.1415926535897932385d0)
       a = a8
       x1=min(abs(y18),abs(y28))
       x2=max(abs(y18),abs(y28))
       n=n1+1-mod(n1,2)		!want n to be odd...
       dx=(x2-x1)/n
       s=voigt(x2,a)/3.
       do i=2,n-1
         x2=x2-dx
         if (mod(i,2).eq.0) then
           x=4./3.
         else
           x=2./3.
         endif
         s=s+x*voigt(x2,a)
       enddo
       s=s+voigt(x1,a)/3.
       rint_v=s*dx !*sqrtpi
       return
       end










       character*90 function NameFn(prefixn,suffn)
c
c			construct file name of form:  prefixn'.'suffn
c					num < 10**9
       integer*4 suffn,num,npos,j
       character*80 prefixn
       character*9 suff

       num = suffn
       suff = '         '
       npos = 1
       do j=8,0,-1
       	if (num .ge. 10**j) then
       		suff(npos:npos) = char(48+int(dble(num/10**j)))
       		npos=npos+1
       		num = int(num - 10**j*int(dble(num/10**j)))
       	else
       	  if (npos .ne. 1) then
       		suff(npos:npos) = '0'
       		npos=npos+1
       	  endif
       	endif
       enddo
       if (npos .eq. 1) suff='0'
c	NameFn=prefixn//'.'//suff
       namefn=prefixn(1:index(prefixn,' ')-1)//'.'//suff
       return
       end

