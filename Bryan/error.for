      subroutine errcalc
c     This code calculates the error propagation of the integrated image
c     over the image of range [low,high].  In addition the integrand may 
c     be multiplied by a function g defined at the bottom of this file.
c
      include 'defs.h'
      
      double precision temp(nfmax,nfmax), ddQinv(nfmax,nfmax), 
     &                 Y(nfmax,nfmax)
      integer ilow, ihigh, i, j, k, ivalue
      double precision g
c
c
c	In order to propagate the error, we need the covariance of the
c	image.  In the gaussian approximation, this would be
c
c             T         -1                      -1
c	<df df > ~= -ddQ   = (- alpha ddS - ddL)
c
c                                                                   T T  
c	         = {f}/alpha -{f} U Y {lambda/alpha(alpha+lambda)} Y U {f}
c
c                                            T  T
c                = {f}U Y {1/(alpha+lambda)}Y  U  {f}
c
c	first form Y(ns,ns)   Y= P*{1/sqrt(z)}*R

	do i=1,ns
	do j=1,ns
	  Y(i,j)=0.0
	do k=1,ns
	   Y(i,j)=Y(i,j) + P(i,k)*R(k,j)/sqrt(Z(k))
	end do
	end do
	end do

c	Now form the matrix temp(ns,ns) = Y*{1/(alpha+lambda)}*Y[T]
	
	do i=1,ns
	do j=1,ns
	   temp(i,j)=0.0
	do k=1,ns
	   temp(i,j)=temp(i,j) +
     &               Y(i,k)*Y(j,k)/(alpha+lambda(k))
	end do
	end do
	end do

c	Now form the product U*temp*U[T]
c
c	first form U*temp(nf,ns) in Y; Y(nf,ns) =  U(nf,ns) temp(ns,ns)
        call dgemm('N','N',nf,ns,ns,1.0d0,Umat,nfmax,temp,nfmax,
     .             0.0d0,Y,nfmax)

c	second, form U*temp*U[T] in ddQinv
        call dgemm('N','T',nf,nf,ns,1.0d0,Y,nfmax,Umat,nfmax,
     .             0.0d0,ddQinv,nfmax)

c	Form {f}*U*Y*{1/(alpha+lambda)}*Y[T]*U[T]*{f} (nf,nf),
c       or Y*U[T] in ddQinv.  This should be the covariance of the image.
	do i=1,nf
	  do j=1,nf
	    ddQinv(j,i)=f(j)*ddQinv(j,i)*f(i)
	  end do
	end do

c	Now <delta f_i delta f_j> = ddQinv(i,j).  The error involved in 
c	the integral
c
c		/
c	G   =	| df f(w) g(w)  
c		/
c
c	May be approximated by  < (delta G)^2 >
c
c                         //
c	< (delta G)^2 > = || dx dy g(x) g(y) <delta f(x) delta f(y)>
c                         //
c
c	Now split the region from w(1) to w(nf) into nregions parts

      do 5 ivalue=1,nregions
        low(ivalue)=w(1)+(ivalue-1)*(w(nf)-w(1))/float(nregions)
        high(ivalue)=w(1)+(ivalue)*(w(nf)-w(1))/float(nregions)

        do i=1,nf
         if (w(i).gt.low(ivalue)) go to 10
        end do
   10   ilow=i

        do i=1,nf
         if (w(i).gt.high(ivalue)) go to 20
        end do
   20   ihigh=i-1

        Gvalue(ivalue)=0.0
        do i=ilow,ihigh
          Gvalue(ivalue)=Gvalue(ivalue)+f(i)*g(w(i))
        end do

        dGvalue(ivalue)=0.0
        do i=ilow,ihigh
        do j=ilow,ihigh
	  dGvalue(ivalue)=dGvalue(ivalue)+g(w(i))*g(w(j))*ddQinv(i,j)
        end do
        end do
	dGvalue(ivalue)=sqrt(dGvalue(ivalue))

        width=(high(1)-low(1))*0.5

 5    continue
      open(unit=77,file='error.dat',status='unknown')
      write(77,*) '#		/'
      write(77,*) '#	G   =	| df f(w) g(w)'
      write(77,*) '#		/'
c23456789012345678901234567890123456789012345678901234567890123456789012
      write(77,*) '#  w_ave     G      width      dG        w_low   w_hi 
     &gh'
      do ivalue=1,nregions
	write(77,110) 0.5*(low(ivalue)+high(ivalue)),Gvalue(ivalue),
     &             width,dGvalue(ivalue),
     &             low(ivalue),high(ivalue)
      end do
 110   format(1x,f7.3,f8.4,2x,f7.3,2x,e10.4,2x,f7.3,2x,f7.3,2x)
c
      return
      end


      double precision function g(x)
c     Define the function to be masked.  For example, if we
c     wanted to calculate the first moment of f, then g=x
c     or x**2 for the second moment. etc.
      double precision x

      g=1.0d0

      return
      end



	subroutine Av_spec

        include 'defs.h'
	double precision fprop(nfmax),fnow(nfmax),unow(nfmax),
     &                   uprop(nfmax),fv1(nfmax)
	double precision Qnow,Qprop,Rweight,acrat(nfmax),loprop,
     &       soprop,lonow,sonow,acratalpha,r1,fprodnow,fprodprop,
     &       acratswp,alphanow,alphaprop,ran1
	integer isweep,i,iu,l,irun

	acratswp=0.0
        r1=ran1(-idummy)
	open(unit=88,file='Av_spec.dat',status='unknown')
	write(88,*) '#  From  subroutine Av_spec'
	write(6,*) ' run   delu   acrat  delalpha acratalpha alpha'

	fprodnow=1.0
	do i=1,nf
	  fave(i)=0.0
	  faves(i)=0.0
	  fnow(i)=f(i)
	  fprodnow=fprodnow*f(i)
	end do
	fprodnow=1.0/sqrt(fprodnow)

	do i=1,ns
	  unow(i)=u(i)
	  uprop(i)=unow(i)
	  acrat(i)=0.0
	end do

	Qnow=alpha*so-lo
	alphanow=alpha

	do 1000 irun=1,nruns+nwarms
	do 999 isweep=1,nsweeps
	do 998 iu=1,ns
	  uprop(iu)=unow(iu) + abs(u(iu))*(ran1(idummy)-0.5)*delu
	  call MXVA(Umat,1,nfmax,uprop,1,fv1,1,nf,ns)
c	   fv1(:)=fv1(:)+Umat(:,iu)*(uprop(iu)-unow(iu))
          do 130 i = 1, nf
	    if(fv1(i).lt.-500.0) then
c	       write(6,*) 'fv1 woops'
	       fv1(i)=-500.0
	    end if
            fprop(i) = model(i)*exp(fv1(i))
 130      continue
c         Form Lo
          loprop=0.0
c         Form the reconstructed image Tf=T*f
          call MXVA (T,1,ndmax,fprop,1,Tf,1,ndata,nf)
          do 140 l=1,ndata
            loprop=loprop+0.5*((Tf(l)-data(l))/error(l))**2
 140      continue

c     Form So
          soprop=0.0
	  fprodprop=1.0
          do 150 i=1,nf
            soprop=soprop-fprop(i)*log(fprop(i)/model(i))+
     &            fprop(i)-model(i)
	  fprodprop=fprodprop*fprop(i)
 150      continue
	  fprodprop=1.0/sqrt(fprodprop)

	  Qprop=alphanow*soprop-loprop
	  Rweight=exp(Qprop-Qnow)*(fprodprop/fprodnow)
	  if(Rweight/(1.0+Rweight).gt.ran1(idummy)) then
c	  make changes in the fields
	    acrat(iu)=acrat(iu)+1.0
	    acratswp=acratswp+1.0
	    lonow=loprop
	    sonow=soprop
	    Qnow=Qprop
	    fprodnow=fprodprop
	    unow(iu)=uprop(iu)
	    do i=1,nf
	      fnow(i)=fprop(i)
	    end do
	  else
c	    try again
	  end if
 998	continue  ! end update of field u(iu)

c	propose a change in alpha
	alphaprop=alphanow + alpha*(ran1(idummy)-0.5)*delalpha
	Qprop=alphaprop*sonow-lonow
	Rweight=exp(Qprop-Qnow)*(alphaprop/alphanow)**(nf/2)
	if(Rweight/(1.0+Rweight).gt.ran1(idummy)) then
c	  make changes in the fields
	  alphanow=alphaprop
	  Qnow=Qprop
	  acratalpha=acratalpha+1.0
	else
c	    try again
	end if

	

	if(irun.gt.nwarms) then
	  do i=1,nf
	    fave(i)=fave(i)+fnow(i)/dw(i)
	    faves(i)=faves(i)+(fnow(i)/dw(i))**2
	    alphaave=alphaave+alphanow
	    alphaaves=alphaaves+alphanow**2
	  end do
	end if

 999	continue  ! end MC sweep

	acratswp=acratswp/float(nsweeps*ns)
	acratalpha=acratalpha/float(nsweeps)
	write(6,556) irun,delu,acratswp,delalpha,acratalpha,alphanow
 556	format(i3,2x,e7.2,2x,f5.3,2x,e7.2,2x,f5.3,2x,f9.3)

c	adjust delu and delalpha
	delu=0.5*delu*(acratswp+1.5)
	acratswp=0.0
	delalpha=0.5*delalpha*(acratalpha+1.5)
	acratalpha=0.0


 1000	continue  ! end MC run


	write(6,*) '  '
	write(6,*) '   i        u(i)               acrat(i) '
	do i=1,ns
	   acrat(i)=acrat(i)/float((nruns+nwarms)*nsweeps)
	   write(6,*) i,u(i),acrat(i)
	end do
	write(6,*) '  '

	write(88,*) '#     w     f_asve(w)   sdev_f_ave     f(w)'
	do i=1,nf
	   fave(i)=fave(i)/float(nruns*nsweeps)
	   faves(i)=faves(i)/float(nruns*nsweeps)
	   fsdev(i)=sqrt(abs(faves(i)-(fave(i))**2))
	   write(88,99) w(i),fave(i),fsdev(i),f(i)/dw(i)
	end do
 99	format(f9.4,2x,e10.4,2x,e10.4,2x,e10.4)
	alphaave=alphaave/float(nruns*nsweeps)
	alphaaves=alphaaves/float(nruns*nsweeps)
	alphasdev=sqrt(alphaaves-alphaave**2)

	return
	end


      double precision FUNCTION RAN1(IDUM) 
      DIMENSION R(97) 
      PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6) 
      PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6) 
      PARAMETER (M3=243000,IA3=4561,IC3=51349) 
      DATA IFF /0/ 
      IF (IDUM.LT.0.OR.IFF.EQ.0) THEN 
        IFF=1 
        IX1=MOD(IC1-IDUM,M1) 
        IX1=MOD(IA1*IX1+IC1,M1) 
        IX2=MOD(IX1,M2) 
        IX1=MOD(IA1*IX1+IC1,M1) 
        IX3=MOD(IX1,M3) 
        DO 11 J=1,97 
          IX1=MOD(IA1*IX1+IC1,M1) 
          IX2=MOD(IA2*IX2+IC2,M2) 
          R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1 
11      CONTINUE 
        IDUM=1 
      ENDIF 
      IX1=MOD(IA1*IX1+IC1,M1) 
      IX2=MOD(IA2*IX2+IC2,M2) 
      IX3=MOD(IA3*IX3+IC3,M3) 
      J=1+(97*IX3)/M3 
      IF(J.GT.97.OR.J.LT.1)PAUSE 
      RAN1=dble(R(J))
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1 
      RETURN 
      END 
