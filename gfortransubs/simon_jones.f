	subroutine sjeqd(n,x,lambda,hsjd)
C
C       This subroutine calculates the Sheather-Jones
C       "solve-the-equation" bandwidth for a Gaussian kernel estimate.
C       (Journal of the Royal Statistical Society, Series B,
C       1991, Vol. 53, pp 683-690). 
C       Double precision is used throughout.
C       IMPORTANT INPUT VARIABLES
C       n = sample size
C       x = vector containing the data
C       lambda = interquartile range (iqr)
C       IMPORTANT OUTPUT VARIABLES
C       hsjd = value of Sheather-Jones bandwidth
C       OTHER IMPORTANT VARIABLES
C       hnml = bandwidth based on iqr assuming data is normal
C              (used as a starting value in this subroutine)
C       chat = estimate of constant c in bandwidth a
C       (i.e. a = c. h**(5/7))
C 
        implicit real*8(a-h,o-z)
	real*8 x(1600),a,anew,lambda,l2,c,c7,tol,h,hsjd,y,y2,e,s,t
        real*8 k1,k2,k3,g,gprime,threen
	integer n,loop
	parameter(tola=1.0d-5,tolg=1.0d-4,loop=15)
	xn = dble(n) 
        k3 = (xn - 1.d0)/dsqrt(2.d0)
        threen = 3.d0*xn
        pi = 3.141592654d0
        rt2pi = dsqrt(2.d0*pi)
        hnml = 0.79d0*lambda/(xn**0.2)
	l2=lambda*lambda
C       *** ESTIMATE VALUE OF c ***
        ac = (0.920d0*lambda)/(xn**(1.d0/7.d0))
        bc = (0.912d0*lambda)/(xn**(1.d0/9.d0))
        sa = 0.d0
        sb = 0.d0
        do 2 i = 1, n-1
           do 3 j = i+1, n
              da = ((x(i) - x(j))/ac)**2
              db = ((x(i) - x(j))/bc)**2
              ea = dexp(-0.5d0*da)
              eb = dexp(-0.5d0*db)
              sa = sa + (da**2 - 6.d0*da + 3.d0)*ea
              sb = sb + (db**3 - 15.d0*db**2 + 45.d0*db - 15.d0)*eb
3          continue
2       continue
        rhat2 = (2.d0*sa)/(xn*(xn-1.d0)*(ac**5)*rt2pi)
        rhat2 = rhat2 + 3.d0/(rt2pi*(xn-1.d0)*(ac**5))
        rhat3 = (-2.d0*sb)/(xn*(xn-1.d0)*(bc**7)*rt2pi)
        rhat3 = rhat3 + 15.d0/(rt2pi*(xn-1.d0)*(bc**7))
        chat = 1.357d0*((dabs(rhat2/rhat3))**(1.d0/7.d0))
        chat7 = chat**7.d0
C	*** USE NEWTON-RAPHSON METHOD TO CALCULATE h ***
C	INITIAL VALUES FOR h AND a
        firsth = hnml 
	ainit = chat*(firsth**(5.d0/7.d0))
99      a = 1.5d0*ainit
	do 4 m=1,loop
		s=0.d0
		t=0.d0
		do 5 i=1,n-1
			do 5 j=i+1,n
				y = (x(i)-x(j))/a
				y2 = y*y
				e = dexp(-0.5d0*y2)
				s = s + (y2*y2 - 6.d0*y2 + 3.d0)*e
				t = t + y2*(y2*y2-10.d0*y2+15.d0)*e
5		continue
		t = 2.0d0 * t
		s = 2.0d0 * s + threen
		k1 = a*a/chat7
		k1 = dsign(dabs(k1)**0.2d0,k1)
		k2 = k3/s
		k2 = dsign(dabs(k2)**0.2d0,k2)
		gprime = -0.2d0*(k2*(t/s-5.0d0) + 7.0d0*k1)
		g = a*(k2-k1)
		anew = a - g/gprime
C		CALCULATE h FROM a
		k1 = anew*anew/chat7
		k1 = dsign(dabs(k1)**0.2d0,k1)
        	h=k1*anew
                if(anew.le.0.0d0)goto 98
		if(dabs(anew-a).lt.tola)goto100
                if(dabs(g).lt.tolg)goto100
		a=anew
4	continue
98      write(6,910)
	write(7,910)
        ainit = 1.5d0*ainit
        goto 99
	stop

100	continue
        hsjd = h
910	format(' NEWTON`S METHOD UNABLE TO FIND SOLUTION TO sjeqd')
        return
	end

