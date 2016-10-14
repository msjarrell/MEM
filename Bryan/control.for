      subroutine control

      
c     This code is the step-size control procedure for the MaxEnt
c     update procedure.  It is so fast that it does not need
c     to be optimized.  The priority is on simplicity.
      
      include 'defs.h'
      
      double precision mumin, mumax, ls, los
      integer i, muflag
      integer eql

      mu = 0.0d0
      mumin=0.0d0
      mumax=0.0d0
      muflag=0
      
c     Generate los, the maximum step size
      los=0.0
      do 10 i=1,nf
         los=los+f(i)
 10   continue
      los=los*step
C$$$
C$$$
c      write(6,*) 'Los, step: ', los,step
      
c     Calculate du(i) and ls
 20   ls=0.0
      do 30 i=1,ns
         du(i)=YIRHS(i)/(mu+lambda(i)+alpha)
         ls=ls+du(i)**2
 30   continue
C$$$
C$$$
c      write(6,*) 'Ls, mu: ', ls,mu

c     Now test if MU-chop is finished
      if (ls.le.los.and.muflag.eq.0.or.eql(ls,los,precision).eq.1) then
c        exit with a successful du
         go to 40
      endif
      
      if (ls.gt.los) then
c        mu must be increased (mu=0.0, first iterate)
         if (muflag.eq.0) then
            mu=10.
            muflag=1
         else
            call uchop (mu,mumin,mumax)
         endif
         go to 20
      else 	
c        mu must be decreased
         call dchop (mu,mumin,mumax)
         go to 20
      endif
 40   continue
      
c      write (6,*) 'Mu=',mu
c      
      return
      end


      subroutine dchop (x,xmin,xmax)
c     This function chops x downward towards xmin

      double precision x, xmin, xmax
      xmax=x
      x=0.5*(xmax+xmin)

      return
      end
      
      
      subroutine uchop (x,xmin,xmax)
c     This function chops x upward towards xmax; however,
c     if xmax is zero, then it just sets x=2x.
      
      double precision x, xmin, xmax
      xmin=x
      if (xmax.lt.1.0e-12) then
         x=2.0*xmin
      else
         x=0.5*(xmin+xmax)
      endif
      
      return
      end
      
      
      integer function eql (x,y,prec)
c     This integer function determines if x and y are equal
c     to some precision set by prec.
      
      double precision x, y, prec
      eql=0
      if (abs((x-y)/(x+y)).lt.prec) eql=1
      
      return
      end

