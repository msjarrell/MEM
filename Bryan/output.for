      subroutine output

      include 'defs.h'

      integer i, l
      double precision co, residual(ndmax)

c     Write a new input file
c      open (unit=4,status='unknown',file='new.input')
c      write (4,*) 'niter, ftr, aflag'
c      write (4,4)  niter, ftr, aflag
c      write (4,*) 'alpha, precision, step,range,nregions'
c      write (4,*)  alphainit, precision, step,range,nregions
c      write (4,*) 'testf,optionm,aim'
c      write (4,*)  testf,optionm,aimaim
c 4    format(i4,x,f10.5,x,i4,x,i4,x,i4)

c     Form chi-squared and the residual of the image
      open(unit=29,file='residual.xvgr',status='unknown')
      write(29,*) '@    default font 0'
      write(29,*) '@    yaxis  label "residual"'
      write(29,*) '@    xaxis  label "l"'
      write(29,*) '@    s0 symbol 2'
      write(29,*) '@    s0 linestyle 0'
      write(29,*) '@    s0 symbol size 0.5'
      write(29,*) '@TYPE xy'
      co=0.0
c     Form the reconstructed data Tf=T*f
      call MXVA (T,1,ndmax,f,1,Tf,1,ndata,nf)
      do l=1,ndata
        co=co+((Tf(l)-data(l))/error(l))**2
        residual(l)=(Tf(l)-data(l))/error(l)
	write(29,*) l,residual(l)
      end do
      aim=co/float(ndata)
      write(29,*) '&'
      write(29,*) '@TYPE xy'
      write(29,*) '0.0   0.0 '
      write(29,*) ndata,'   0.0 '
      write(29,*) '&'


c	Renormalize the image for printing
	if(optionm.eq.0) then
	  do i=1,nf
	    u(i)=f(i)
	    f(i)=u(i)/dw(i)
	    model(i)=mr(i)
	    u(i)=f(i)/model(i)
          end do
	else
	  do i=1,nf
	    u(i)=f(i)
	    f(i)=u(i)*mr(i)
	    model(i)=mr(i)
          end do
	end if

      if(aflag.eq.2.or.aflag.eq.0) weight=2.0*weight
c     Write the present image
      if(ixvgr.eq.1) then
      write(7,*) '#weight: ',weight,' aflag:',aflag
      write(7,*) '# Ngood: ',ngood
c      if(aflag.eq.2.or.aflag.eq.0) write(7,*) '#alpha: ',alpha
      write(7,*) '#alpha: ',alpha
c      if(aflag.eq.2.or.aflag.eq.0) write(7,*) '#sigma: ',sig,' aim:',aim
      write(7,*) '#sigma: ',sig,' aim:',aim
	write(7,*) '@    s2 symbol 2'
	write(7,*) '@    s2 linestyle 0'
	write(7,*) '@    s2 symbol size 0.5'
	write(7,*) '@    legend on'
	write(7,*) '@    legend x1 0.66792'
	write(7,*) '@    legend y1 0.775444'
	write(7,*) '@    legend string 0 "image"'
	write(7,*) '@    legend string 1 "model"'
	write(7,*) '@    legend string 2 "errors"'
	write(7,*) '@TYPE xy'
      end if
      do i=1,nf
         write (7,12) w(i),f(i),dw(i)
         write (99,12) w(i),max(1.0e-3,f(i)),dw(i)
         write (98,12) w(i),max(1.0e-7,f(i)),dw(i)
      end do
 12   format(f14.8,2x,e16.9,2x,f14.8)
      if(ixvgr.eq.1) then
	write(7,*) '&'
	write(7,*) '@TYPE xy'
      else
	write(7,*) '   '
      end if
	 do i=1,nf
	  write(7,*) w(i),model(i)
	end do
      if(ixvgr.eq.1) then
	write(7,*) '&'
	write(7,*) '@TYPE xydxdy'
	do i=1,nregions
	  write(7,55) 0.5*(low(i)+high(i)),Gvalue(i),
     &             width,dGvalue(i),
     &             low(i),high(i)
	end do
      end if
 55   format(1x,f7.3,f8.4,2x,f7.3,2x,e10.4,2x,f7.3,2x,f7.3,2x)

      if(aflag.eq.1.or.aflag.eq.3) then
         open (unit=8,status='unknown',file='posterior.xvgr')
         write(8,3000) 
 3000	 format('@    default font 0')
         write(8,3001) 
 3001    format('@    yaxis  label "P(alpha|D,m)"')
         write(8,3002) 
 3002    format('@    xaxis  label "alpha"')
         do i = 1,iter
            if(padm(i).gt.1.0e-29) write (8,*) alpt(i),padm(i)
         end do
      endif

      write(6,*) ' '
      if(iuflow.eq.1) write(6,*) 'An image pixel fell below e^-60 '
      write(6,*) 'Alpha: ', alpha
      write(6,*) 'Sigma: ', sig
      write(6,*) 'Ngood: ', ngood
      write (6,*) 'Chisq: ', co
      write(6,*) 'aim=',aim
      write(6,*) 'weight: ',weight
      write(22,*) weight
c
      return
      end
