	program rdbins
c	This code is used to process the binned data into data which can be 
c	used by readmc, the preprocessor of the analytic continuation code.
c	This this code takes the binned data, and produces the averages, 
c	and covariances.  The binn size is variable within the code, i.e. 
c	bins can be averaged, and multiple data files can be opened and 
c	averaged.  In addition, rdbins is capable of qualifying the data 
c	by performing a variety of tests.   These include tests to measure 
c	the kurtosis and skew of the distribution, and a calculation of 
c	the probability that the data is gaussian (normally distributed). 
c	This last probability should be at least 0.5.
c
c	f77 -r8 -O -o rdbins rdbins.for  

	implicit none
	integer nl,i,j,k,l,skip,nuse,ind,naver,ihisto,isymm
 	integer nlmax,nbinsmax,nbins,nbinsuse,ibin
	integer nfiles,nfile,nl_bin
	parameter(nlmax=200,nbinsmax=4000)
	real beta,tau,dtau
	real data(nlmax,nbinsmax),data2(nbinsmax),aver(nlmax)
	real ave,adev,sdev,var,skew,curt,chsq,prob,error
	real store(nlmax),ds(nlmax,nbinsmax)
        character*72 linemc

      write(6,11) 
 11   format('               enter output data filename:  ',$)
      read(5,201) linemc
      if(index(linemc(1:72),'bins').ne.0 )then
	  write(6,*) 'This is an input file'
	  stop
      end if
      open(unit=9,file=linemc,status='unknown')

      write(6,12) 
 12   format('          enter the number of input files:  ',$)
      read(5,*) nfiles

      if(nfiles.eq.1)then
          write(6,13) 
 13       format('      enter the number of bins to be used:  ',$)
	  read(5,*) nbinsuse
      else
	  nbinsuse=nbinsmax
      end if

c	skip is used when the data is too dependent.  If
c	skip is greater than one, some of the time slices in
c	each bin will be skipped.  Typically (now) skip=1
      write(6,14) 
 14   format('                               enter skip:  ',$)
      read(5,*) skip	

      write(6,15) 
 15   format('enter course grain size for bin averaging:  ',$)
      read(5,*) naver	

      ihisto=0
      write(6,16) 
 16   format('                  Print histograms? (y/n):  ',$)
      read(5,201) linemc 
      if(index(linemc(1:72),'y').ne.0) ihisto=1

      isymm=0
      write(6,17) 
 17   format('    Should the data be symmetrised? (y/n):  ',$)
      read(5,201) linemc 
      if(index(linemc(1:72),'y').ne.0) then
	isymm=1
	write(9,*) 'symmetric storage mode'
      end if

      nbins=0
      ibin=0
      do 999 nfile=1,nfiles

      	write(6,18) 
 18   	format('            enter the input data filename:  ',$)
	read(5,201) linemc
      	open(unit=92,file=linemc,status='old')

 400    read(92,201) linemc
        if(index(linemc(1:72),'Results').ne.0 ) then
          read(92,*) nl,nl_bin,beta
	  nuse=(nl_bin-1)/skip+1
	  dtau=beta/float(nl)
	  if (isymm.eq.1)nuse=nuse/2+1
        else
         goto 400
        endif
	
	do i=1,nbinsuse
c	  find a bin header
 407      read(92,201,end=410) linemc
          if(index(linemc(1:72),'bin').ne.0) then
	    ibin=ibin+1
c	    readin a bin
	    read(92,*) (store(j),j=1,nl_bin)
	    if(isymm.eq.0) then
	      do j=1,nuse
		ds(j,ibin)=store(1+skip*(j-1))
		data(j,ibin)=0.0
              end do
	    else
	      ds(1,ibin)=store(1)
	      do j=2,nuse
		ds(j,ibin)=0.5*(store(1+skip*(j-1))+ ! Ilja Makkonen
     &			store(1+nl_bin-skip*(j-1))  )  ! j --> nl+2-j
		data(j,ibin)=0.0
	      end do
	    end if
          else
            goto 407
          endif
        end do
 410	continue

 999  continue

      nbins=ibin
      nbins=nbins/naver
      write(6,19) nbins
 19   format('                                    nbins:  ',i3)
      write(6,20) nl
 20   format('                                       nl:  ',i3)

      write(9,*) 'Results for nl,nbins,beta'
      write(9,*) (nl_bin-1)/skip+1,nbins,beta
      write(9,*) ' '


c     Now do bin averaging.  The data in bins ibin > naver*int(nbins/naver)
c     is not used.
      ind=0
      do i=1,nbins
	  do j=1,naver
	    ind=ind+1
	     do k=1,nuse
	        data(k,i)=data(k,i)+ds(k,ind)/float(naver)
             end do
          end do
      end do
	

c     Now analyze the data (the analysis tools are described below)
c23456789012345678901234567890123456789012345678901234567890123456789012
      write(9,*) '   tau      Gtau(n)           +/-       curt  skew  pr  
     &ob'
      do l=1,nuse
	  do i=1,nbins
	    data2(i)=data(l,i)
          end do
	  call moment(data2,nbins,ave,adev,sdev,var,skew,curt)
	  tau=dtau*float(skip*(l-1)) ! thanks Ilja Makkonen
          aver(l)=ave
	  if(adev.gt.1.0e-20) then
	   call histo(data2,nbinsmax,nbins,ave,var,sdev,chsq,prob,ihisto)
	  end if
	  error=sdev/sqrt(float(nbins))
	  write(9,202) tau,ave,error,curt,skew,prob
      end do
      write(9,*) ' '
      write(9,*) ' Covariance data '
      write(9,*) '   i     j     Gij-GiGj        +/-'
      call covar(data,nlmax,nbinsmax,nuse,nbins,aver,9)
 202  format(x,f8.5,x,e14.7,x,e14.7,x,f5.2,x,f5.2,x,f5.3)
 201  format(a72)

      stop
      end


c	This series of subroutine calls is used to process
c	binned data stored in an array meas(index,binnumber).
c	
c	Subroutine moment calculated the moments of the data
c	for a particular index.  The data for the particular
c	index must be transferred to a 1-d array: 
c	
c	     data(:)=meas(index,:)
c
c	Then execute
c
c	call moment(data,nbins,ave,adev,sdev,var,skew,curt)

c	where nbins is the number of bins of data.  On output
c	ave is the bin average, adev the average deviation,
c	sdev the standard deviation, var the variance, skew
c	the skew and curt the kurtosis of the distribution.
c	(for a more detailed discussion of curt and skew, see
c	Numerical Recipes)

c
c
c	One can also generate a histogram of this same data, and
c	compare it to a gaussian distribution.  On input, nbinsmax
c	is the declared dimension of data, nbins the number of bins,
c	ave, var, and sdev that calculated by moments.  On output
c	prob is the probability that the binned data is normal
c	(fora discussion of chsq see Numerical Recipes).  prob
c	should be about 0.5 or bigger.
	
c	call histo(data,nbinsmax,nbins,ave,var,sdev,chsq,prob)

c	If you generate the binned data in two ways, say shuffled and
c	sequential binning, then you can compare the two data sets
c	with an FTEST.  If the binning procedure removed all the 
c	corellations between adjacent bins, then the two data set
c	are identical, and the calculated probability would be 1.

c	Use the ftest on the two data sets of G(tau)
c	write(9,*) ' '
c	write(9,*) ' FTEST DATA '
c	write(9,*) ' '
c	write(9,*) '   l      f -----> prob  '
c      	do 35 i=1,nlg
c	   data(:)=gmeas1(i,:)
c	   data2(:)=gmeas2(i,:)
c	   call ftest(data,run,data2,run,f,prob)
c	   write(9,87) i,f,prob
c 35	continue
c 87      format(1x,i3,3x,f7.4,3x,f8.5)

c	Of course in order to analytically continue numerical data
c	the covariance matrix is necessary.  The following takes
c	the data in meas, and the associated vector of averages
c	ave(index) and calculate and writes to standard output
c 	the covariance matrix.

c	Now calculate the covariance using 
c	covar(meas,nm1,n1,n2,ave,nm2)
c
c	Where 
c	meas(nm1,nm2) contains the measurements
c	nm1 declared dimension of meas
c	nm2 declared dimension of meas
c	n1 dimension of meas and ave
c	n2 dimension of meas
c	ave average over runs measurments
c
c	Thus, with partilce-hole symmetry, we could take
c	n1=nlg/2 in order to save space

c	call covar(meas,indexmax,nl,nbins,ave,binsmax)


      SUBROUTINE MOMENT(DATA,N,AVE,ADEV,SDEV,VAR,SKEW,CURT) 
      DIMENSION DATA(N) 
      IF(N.LE.1) then
	write(6,*) 'N must be at least 2' 
	return
      end if
      S=0. 
      DO 11 J=1,N 
        S=S+DATA(J) 
11    CONTINUE 
      AVE=S/N 
      ADEV=0. 
      VAR=0. 
      SKEW=0. 
      CURT=0. 
      DO 12 J=1,N 
        S=DATA(J)-AVE 
        ADEV=ADEV+ABS(S) 
        P=S*S 
        VAR=VAR+P 
        P=P*S 
        SKEW=SKEW+P 
        P=P*S 
        CURT=CURT+P 
12    CONTINUE 
      ADEV=ADEV/N 
      VAR=VAR/(N-1) 
      SDEV=SQRT(VAR) 
      IF(VAR.NE.0.)THEN 
        SKEW=SKEW/(N*SDEV**3) 
        CURT=CURT/(N*VAR**2)-3.   
      ELSE 
        write(6,*) 'no skew or kurtosis when zero variance' 
	return
      ENDIF 
      RETURN 
      END 

      SUBROUTINE FTEST(DATA1,N1,DATA2,N2,F,PROB)
      DIMENSION DATA1(N1),DATA2(N2)
      CALL AVEVAR(DATA1,N1,AVE1,VAR1)
      CALL AVEVAR(DATA2,N2,AVE2,VAR2)
      IF(VAR1.GT.VAR2)THEN
         F=VAR1/VAR2
         DF1=N1-1
         DF2=N2-1
      ELSE
         F=VAR2/VAR1
         DF1=N2-1
         DF2=N1-1
      ENDIF
      PROB = BETAI(0.5*DF2,0.5*DF1,DF2/(DF2+DF1*F))
     *    +(1.-BETAI(0.5*DF1,0.5*DF2,DF1/(DF1+DF2/F)))
      RETURN
      END


      SUBROUTINE AVEVAR(DATA,N,AVE,VAR)
      DIMENSION DATA(N)
      AVE=0.0
      VAR=0.0
      DO 11 J=1,N
         AVE=AVE+DATA(J)
 11   CONTINUE
      AVE=AVE/N
      DO 12 J=1,N
         S=DATA(J)-AVE
         VAR=VAR+S*S
 12   CONTINUE
      VAR=VAR/(N-1)
      RETURN
      END
      
      FUNCTION BETAI(A,B,X)
      IF(X.LT.0..OR.X.GT.1.)then
	write(6,*) 'bad argument X in BETAI'
	return
      end if
      IF(X.EQ.0..OR.X.EQ.1.)THEN
         BT=0.
      ELSE
         BT=EXP(GAMMLN(A+B)-GAMMLN(A)-GAMMLN(B)
     *        +A*ALOG(X)+B*ALOG(1.-X))
      ENDIF
      IF(X.LT.(A+1.)/(A+B+2.))THEN
         BETAI=BT*BETACF(A,B,X)/A
         RETURN
      ELSE
         BETAI=1.-BT*BETACF(B,A,1.-X)/B
         RETURN
      ENDIF
      END
      
      
      FUNCTION BETACF(A,B,X)
      PARAMETER (ITMAX=100,EPS=3.E-7)
      AM=1.
      BM=1.
      AZ=1.
      QAB=A+B
      QAP=A+1.
      QAM=A-1.
      BZ=1.-QAB*X/QAP
      DO 11 M=1,ITMAX
         EM=M
         TEM=EM+EM
         D=EM*(B-M)*X/((QAM+TEM)*(A+TEM))
         AP=AZ+D*AM
         BP=BZ+D*BM
         D=-(A+EM)*(QAB+EM)*X/((A+TEM)*(QAP+TEM))
         APP=AP+D*AZ
         BPP=BP+D*BZ
         AOLD=AZ
         AM=AP/BPP
         BM=BP/BPP
         AZ=APP/BPP
         BZ=1.
         IF(ABS(AZ-AOLD).LT.EPS*ABS(AZ))GO TO 1
 11   CONTINUE
      write(6,*) 'A or B too big, or ITMAX too small'
 1    BETACF=AZ
      RETURN
      END
      
      

      subroutine histo(data,runmax,run,avg,var,sdev,chsq,prob,iq)
      integer run,runmax,iq
      real hist(50), curve(50)
      real data(runmax),width,sdev
      
      do 10 i = 1, 50
         curve(i) = 0.0
         hist(i) = 0.0
 10   continue
      
      nx=20
      width=6.0
      xbeg=avg-0.5*width*sdev
      delx=width*sdev/float(nx)
      xend = xbeg + nx*delx
      do 40 i = 1, run
         x=data(i)
         if(x.lt.xbeg) goto 40
         if(xend.lt.x) goto 40
         j = (x-xbeg)/delx + 1
         hist(j) = hist(j) + 1.
 40   continue
      
      call gcurve(xbeg,delx,nx,avg,var,run,curve)
      
      call chsone(hist,curve,nx,2,df,chsq,prob)
      if(iq.eq.1)call phist(xbeg,delx,nx,hist,curve)
      
      return
      end
      
      
      subroutine gcurve(xbeg,delx,nx,xm,var,ndata,curve)
      dimension curve(*)
      
      sigma = sqrt(var)
      
      do 10 j = 1, nx
         x = xbeg + (j-0.5)*delx
         c = exp(-1.0*(x-xm)**2/(2.*var))/(2.5066*sigma)
         curve(j) = c*ndata*delx
 10   continue
      
      return
      end
      
      SUBROUTINE CHSONE(BINS,EBINS,NBINS,KNSTRN,DF,CHSQ,PROB)
      DIMENSION BINS(NBINS),EBINS(NBINS)
      DF=NBINS-1-KNSTRN
      CHSQ=0.
      DO 11 J=1,NBINS
         IF(EBINS(J).LE.0.) goto 11
         CHSQ=CHSQ+(BINS(J)-EBINS(J))**2/EBINS(J)
 11   CONTINUE
      PROB=GAMMQ(0.5*DF,0.5*CHSQ)
      RETURN
      END

      FUNCTION GAMMQ(A,X)
      IF(X.LT.0..OR.A.LE.0.)return
      IF(X.LT.A+1.)THEN
        CALL GSER(GAMSER,A,X,GLN)
        GAMMQ=1.-GAMSER
      ELSE
        CALL GCF(GAMMCF,A,X,GLN)
        GAMMQ=GAMMCF
      ENDIF
      RETURN
      END
      
      
      SUBROUTINE GCF(GAMMCF,A,X,GLN)
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A)
      GOLD=0.
      A0=1.
      A1=X
      B0=0.
      B1=1.
      FAC=1.
      DO 11 N=1,ITMAX
         AN=FLOAT(N)
         ANA=AN-A
         A0=(A1+A0*ANA)*FAC
         B0=(B1+B0*ANA)*FAC
         ANF=AN*FAC
         A1=X*A0+ANF*A1
         B1=X*B0+ANF*B1
         IF(A1.NE.0.)THEN
            FAC=1./A1
            G=B1*FAC
            IF(ABS((G-GOLD)/G).LT.EPS)GO TO 1
            GOLD=G
         ENDIF
 11   CONTINUE
      write(6,*) 'A too large, ITMAX too small'
 1    GAMMCF=EXP(-X+A*ALOG(X)-GLN)*G
      RETURN
      END
      
      
      SUBROUTINE GSER(GAMSER,A,X,GLN)
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A)
      IF(X.LE.0.)THEN
         IF(X.LT.0.)return
         GAMSER=0.
         RETURN
      ENDIF
      AP=A
      SUM=1./A
      DEL=SUM
      DO 11 N=1,ITMAX
         AP=AP+1.
         DEL=DEL*X/AP
         SUM=SUM+DEL
         IF(ABS(DEL).LT.ABS(SUM)*EPS)GO TO 1
 11   CONTINUE
      write(6,*) 'A too large, ITMAX too small'
 1    GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)
      RETURN
      END
      
      
      FUNCTION GAMMLN(XX)
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *     -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
         X=X+ONE
         SER=SER+COF(J)/X
 11   CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END

      subroutine phist(xbeg,delx,nx,hist,curve)
      character*1 blank, cross, ast
      character*50 zline
      dimension hist(*), curve(*)
      data blank, cross, ast/' ', 'X', '*'/
      
      scfact = 1.0
      hmax = 0.
      do 2 j = 1, nx
         if (hmax.gt.hist(j)) goto 2
         hmax = hist(j)
 2    continue
      
 3    if (hmax.le.50.) goto 5
      hmax = hmax/2.
      scfact = scfact/2.
      goto 3
      
 5    do 40 j = 1, nx
         x = xbeg + (j-0.5)*delx
         do 10 i = 1, 50
            zline(i:i) = blank
 10      continue
         k = hist(j)*scfact + 0.5
         if (k.le.0) go to 25
         do 20 i = 1, k
            zline(i:i) = cross
 20      continue
 25      k = curve(j)*scfact + 0.5
         if(k.gt.50) goto 40
         zline(k:k) = ast
         write(6,1000) x, nint(hist(j)), zline
 40   continue
      
 1000 format(1h ,f7.3,i5,2h  ,(a))
      return
      end

	subroutine covar(data,i1,i2,n1,n2,ave,lun)
	integer n1,n2,i,j,k,i1,i2,lun
	real data(i1,i2),ave(i1),cij,cij2,sigg
c	thic code calculates the cavaraince of data(n1,n2) where
c	n1 is the data index and n2 is the bin index.  ave(n1)
c	contains the average values of the measurement.  The covariance,
c	and its error are printed out to logical unit 9.

	do 50 i=1,n1
	do 50 j=i,n1
	  cij=0.0
	  cij2=0.0
	  do 55 k=1,n2
	   cij=cij + (data(i,k)-ave(i))*(data(j,k)-ave(j))
           cij2=cij2+((data(i,k)-ave(i))*(data(j,k)-ave(j)))**2
 55	  continue
	  cij=cij/float(n2)
	  cij2=cij2/float(n2)
	  sigg=(cij2-cij**2)/float(n2-1)
c	  if(sigg.lt.0.0) write(6,*) ' - variance of cij ',sigg
	  sigg=sqrt(abs(sigg))
	  write(lun,90) i,j,cij,sigg
 50	continue
 90	format(1x,i4,2x,i4,2x,e14.7,2x,e14.7)

	return
	end

