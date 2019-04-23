      SUBROUTINE RT(x, d)
C
C     generate a random number x with t distribution of degree of freedom  
C     d
      integer d,i
      REAL x, tmp(d),y
      external rnorm
C
      y=0
      do 10 i=1,d
         call rnorm(x,tmp(i))
         y=y+tmp(i)**2
 10   continue
      x=x/sqrt(y/d)
      return
      end
