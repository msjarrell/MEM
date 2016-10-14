      subroutine update

c     For a fixed alpha, solves the MaxEnt equation by the
c     quasi-Newton method with a Marquardt cut-off on the step size.

      include 'defs.h'

      double precision G(nfmax), XK(nfmax,nfmax), A(nfmax,nfmax)
      double precision YI(nfmax,nfmax), YIT(nfmax,nfmax), 
     &                 MYIT(nfmax,nfmax)
      double precision ds(nfmax), dl(nfmax)
      double precision temp1(nfmax,nfmax), temp2(nfmax,nfmax)      
      double precision fv1(nfmax),fv2(nfmax)

      external funct

      double precision norm1, norms, norml,dnrm2
      integer l,j,i,ier
      double precision ddot
      double precision ax,bx,tol,cx,xmax,golden,probmax,test1

c     Beginning of loop
 1    continue

c     Form the reconstructed data Tf=T*f
      call MXVA (T,1,ndmax,f,1,Tf,1,ndata,nf)

c     Form the gradient of L
      do 10 l=1,ndata
         dl(l) = (Tf(l)-data(l))/error(l)**2
 10   continue

c     Form G = sigma*v[T]*dl
      call MXVA(v,ndmax,1,dl,1,G,1,ns,ndata)
      do 30 j = 1,ns
         G(j) = s(j)*G(j)
 30   continue

c     Form RHS = -alpha*u - G
      do 40 i = 1, ns
         RHS(i) = - alpha*u(i) - G(i)
 40   continue

c     Construct the Jacobian in its diagonal space
c     ref: Golub & vanLoan, Matrix Computations, pg.469; Bryan loc.cit.
c     First form K =U[T]*f*U
      do 70 i=1, nf
         do 60 j = 1, ns
            temp1(i,j) = f(i)*Umat(i,j)
 60      continue
 70   continue

c     XK(ns,ns) = Umat(ns,nf)*temp1(nf,ns)
c      call MXMA(Umat,nfmax,1,temp1,1,nfmax,XK,1,nfmax,ns,nf,ns)
      call dgemm('T','N',ns,ns,nf,1.0d0,Umat,nfmax,temp1,nfmax,
     1           0.0d0,XK,nfmax)
      do i=1,nfmax
      do j=1,nfmax
        temp1(i,j)=XK(i,j)
      end do
      end do
      call ssvdc(temp1,nfmax,ns,ns,Z,fv1,P,nfmax,temp2,
     1           nfmax,fv2,10,ier)

c     Form A = sqrt(Z)*P[T]*M*P*sqrt(Z)
      do 90 i = 1, ns
         do 80 j = 1, ns
	    if(Z(j).lt.0.0) write(6,*) 'negative z(j), j=',j
            temp2(i,j) = P(i,j)*sqrt(Z(j))
 80      continue
 90   continue
c      call MXMA(XM,1,nfmax,temp2,1,nfmax,temp1,1,nfmax,ns,ns,ns)
      call dgemm('N','N',ns,ns,ns,1.0d0,XM,nfmax,temp2,nfmax,0.0d0,
     1      temp1,nfmax)
c      call MXMA(temp2,nfmax,1,temp1,1,nfmax,A,1,nfmax,ns,ns,ns)
      call dgemm('T','N',ns,ns,ns,1.0d0,temp2,nfmax,temp1,nfmax,0.0d0,
     1           A,nfmax)
      do i=1,nfmax
      do j=1,nfmax
        temp1(i,j)=A(i,j)
      end do
      end do
	call ssvdc(temp1,nfmax,ns,ns,LAMBDA,fv1,R,nfmax,temp2,
     1             nfmax,fv2,10,ier)

c     Form Y[-1] = R[T]*sqrt(Z)*P[T] and Y[-T] = P*sqrt(Z)*R
      do 110 i = 1, ns
         do 100 j = 1, ns
            temp1(i,j) = sqrt(Z(i))*P(j,i)
            temp2(i,j) = sqrt(Z(i))*R(i,j)
 100     continue
 110  continue
c      call MXMA(R,nfmax,1,temp1,1,nfmax,YI,1,nfmax,ns,ns,ns)
      call dgemm('T','N',ns,ns,ns,1.0d0,R,nfmax,temp1,nfmax,0.0d0,
     1      YI,nfmax)
c      call MXMA(P,1,nfmax,temp2,1,nfmax,YIT,1,nfmax,ns,ns,ns)
      call dgemm('N','N',ns,ns,ns,1.0d0,P,nfmax,temp2,nfmax,0.0d0,
     1           YIT,nfmax)

c     Form MYIT = M*Y[-T]
c      call MXMA(XM,1,nfmax,YIT,1,nfmax,MYIT,1,nfmax,ns,ns,ns)
      call dgemm('N','N',ns,ns,ns,1.0d0,XM,nfmax,YIT,nfmax,0.0d0,
     1           MYIT,nfmax)

c     Form YIRHS = -alpha*Y[-1]*u - Y[-1]*g
      call MXVA(YI,1,nfmax,RHS,1,YIRHS,1,ns,ns)


c     Adjust step-size so it is consistent with the linerization
c          returns Y[-1}*du = Y[-1]*(-alpha*u-G)/(alpha*I+LAMBDA+mu)
      call control


c     Compute new image
c        u = u + du
c        du = (-alpha*u - g - M*Y[-T]*(Y[-1]*du))/(alpha+mu)
      call MXVA(MYIT,1,nfmax,du,1,fv1,1,ns,ns)
      do 120 i = 1, ns
         u(i) = u(i) + (RHS(i) - fv1(i))/(alpha+mu)
 120  continue
      call MXVA(Umat,1,nfmax,u,1,fv1,1,nf,ns)
      do 130 i = 1, nf
	 if(fv1(i).lt.-60.0) then
c	    write(6,*) 'fv1 woops'
	    iuflow=1
	    fv1(i)=-60.0
	 end if
         f(i) = model(i)*exp(fv1(i))
	 f_real(i) = model(i)*exp(fv1(i))/dw(i)
 130  continue

c     Form Lo
      lo=0.0
c        Form the reconstructed image Tf=T*f
      call MXVA (T,1,ndmax,f,1,Tf,1,ndata,nf)
      do 140 l=1,ndata
         lo=lo+0.5*((Tf(l)-data(l))/error(l))**2
 140  continue
      aim=2.0*lo/float(ndata)

c     Form So
      so=0.0
      do 150 i=1,nf
         so=so-f(i)*log(f(i)/model(i))+f(i)-model(i)
 150  continue

c     Measure tests (remember the metric! K=Y[-T]*Y[-1])
      call MXVA(YI,1,nfmax,u,1,ds,1,ns,ns)
      call MXVA(YI,1,nfmax,G,1,dl,1,ns,ns)
      norms = dnrm2(ns,ds,1)
      norml = dnrm2(ns,dl,1)

      call MXVA(YI,1,nfmax,RHS,1,fv1,1,ns,ns)
      test1= 2.*ddot(ns,fv1,1,fv1,1)/(alpha*norms+norml)**2

      if (mu.gt.0.001.or.test1.gt.testf) go to 1

c     Convergence has been obtained
      if(cflag.eq.0)then
        ax=0.001
        bx=alpha
        cx=1.0e6
        tol = 1.e-04
        probmax = golden(ax,bx,cx,funct,tol,xmax)
        alpha0 = xmax
c       The next value of alpha is given by
        alpha=(1.0-ftr)*alpha1 +ftr*alpha0
      else
	alpha0=alpha
      end if

c     Now determine the new value of ngood.
      norm1=0.0
      do 160 i= 1,ns
         norm1=norm1+(LAMBDA(i))/(LAMBDA(i)+alpha)
 160  continue
       ngood=norm1
c      write (6,*) ' ngood = ',ngood


c     For auto error scaling
      sig=sqrt(2.0*lo/(ndata-norm1))

c     Now compute to posterior probability for the present value of alpha
c         and accumulate average values
      post = 1.
      do 180 i = 1, ns
         post = post*alpha/(alpha+LAMBDA(i))
 180  continue
      if(aflag.eq.0.or.aflag.eq.1) then
c     Implement the Jeffrey prior.
        post = sqrt(post)*exp(alpha*so-lo)/alpha
      else
        post = sqrt(post)*exp(alpha*so-lo)
      end if

c      write(6,200) alpha, post

      weight = weight + post*abs(alpha-alpha1)
      if(aflag.eq.1.or.aflag.eq.3) then
         padm(iter) = post 
         alpt(iter) = alpha1
         alphabar = alphabar + post*abs(alpha-alpha1)*alpha
         sigbar = sigbar + post*abs(alpha-alpha1)*sig
         ngbar = ngbar + post*abs(alpha-alpha1)*norm1
         do 190 i = 1, nf
            fbar(i) = fbar(i) + post*abs(alpha-alpha1)*f(i)
 190     continue
      endif
c
      return
 200  format(' P(alpha=',f10.2,'|D,m): ',e10.3)
      end


      double precision function funct(x)
c     assumes the Jeffery prior
      include 'defs.h'
      double precision x
      integer i
      if(x.le.0.0) then
	funct=-1.0d100
	write(6,*) 'woops, too small'
	return
      end if
      if(x.ge.1.0d200) then
	write(6,*) 'woops,too large'
	stop
      end if
      funct = 1.
      do 10 i = 1, ns
         funct = funct/(x+lambda(i))
 10   continue
c     the mimus sign makes the search look for a maximum
c         Assumes the Jeffery prior.  Its effect is divided out so 
c         alpha can equal zero
      if(aflag.eq.0.or.aflag.eq.1) then
        funct = -sqrt(funct*(x**(ns-2)))*exp(x*so-lo)
      else
        funct = -sqrt(funct*(x**(ns)))*exp(x*so-lo)
      end if

c
      return
      end
      

      double precision FUNCTION GOLDEN(AX,BX,CX,F,TOL,XMIN)
      implicit double precision (a-h,o-z)
      PARAMETER (R=.61803399,C=.38196602)
      X0=AX
      X3=CX
      IF(ABS(CX-BX).GT.ABS(BX-AX))THEN
        X1=BX
        X2=BX+C*(CX-BX)
      ELSE
        X2=BX
        X1=BX-C*(BX-AX)
      ENDIF
      F1=F(X1)
      F2=F(X2)
1     IF(ABS(X3-X0).GT.TOL*(ABS(X1)+ABS(X2)))THEN
        IF(F2.LT.F1)THEN
          X0=X1
          X1=X2
          X2=R*X1+C*X3
          F0=F1
          F1=F2
          F2=F(X2)
        ELSE
          X3=X2
          X2=X1
          X1=R*X2+C*X0
          F3=F2
          F2=F1
          F1=F(X1)
        ENDIF
      GOTO 1
      ENDIF
      IF(F1.LT.F2)THEN
        GOLDEN=F1
        XMIN=X1
      ELSE
        GOLDEN=F2
        XMIN=X2
      ENDIF
      RETURN
      END


