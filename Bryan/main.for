      program main

c     This is the main block of the MaxEnt code
c     modified to do both automatic error, alpha scaling, and alpha 
c     integration. Implemented by Mark Jarrell and J.E. Gubernatis, 1990
c     following J.Skilling and R.K. Bryan, Mon. Not. R. Ast. Soc., 1984, 
c     211, 111-124 and R.K. Bryan, Eur. Biophys, J, 1990, 18, 165-174.
c     Intrinsic correlation function were added about 1/5/91 (deleted)
c     Historic maxent was added in 1991.
c
c     The main function of the algorithm is controlled by aflag
c
c	aflag	function
c	0	Classic		with Jeffrey Prior for tighter fit.
c	1	Bryan		with Jeffrey Prior for tighter fit.
c	2	Classic
c	3	Bryan		 
c	4	Historic MEM 	image fit to chi^2 = ndata


      include 'defs.h'
      double precision test
      integer i
      character*30 string


c     First read the data and then initialize derived quantities
      call init

      if (niter.gt.npmax) then
         write(6,*) 'main: niter > npmax'
         stop
      endif
      if(iprint.gt.0.and.aflag.ne.4) then
	write(6,*) '  '
	write(6,*) ' iteration #     Alpha '
      end if
 15   format(5x,i4,3x,f13.5)
      if(aflag.ne.4) then
c     Then iterate the image
c        find the value of alpha that maximizes the posterior probability
c        for the image
      do 10 i=1,niter
         iter = i
         alpha1=alpha
         if(iprint.gt.0) write(6,15) i,alpha
         call update
         write(string,*) 'alpha= ',alpha
         test = (alpha-alpha1)/(alpha+alpha1)
         if (alpha.ne.alpha1 .and. abs(test).le.0.0005) then
            if(iprint.gt.0) write(6,*) 'Alpha converged'
	    cflag=1
            goto 20
         endif
 10   continue
 20   alpha1 = alpha
c	If we are to integrate over the posterior, we should perform the
c	error propagation now, one we have the most likely image
	 if(aflag.eq.1.or.aflag.eq.3) then
	   call errcalc
	   if(iasmc.eq.1) call Av_spec
	end if

c     Alpha is assumed to be initialized to a large value.  As alpha 
c     converges, values of padm(alpha) are thus generated for large
c     values of alpha. This section of the code computes padm for 
c     small values of alpha.  The 16, 21 and 20.0 are arbitrary. 

      alphamax = alpha
      if (aflag.eq.1.or.aflag.eq.3) then
c        find the posterior probability of the image for other values 
c        of alpha and integrate the image over it.
         if (niter+19.gt.npmax) then
            write(6,*) 'main: niter+16 > npmax'
            stop
         endif
         do 30 i = 1, 19
            iter = iter + 1
            alpha1 = alpha
            alpha = alphamax*float(21-i)/20.0
            call update
            write(string,*) 'alpha= ',alpha
            if(iprint.gt.0) write(6,15) iter,alpha
 30      continue
         do 40 i = 1, nf
            f(i) = fbar(i)/weight
 40      continue
         do 50 i = 1, iter
            padm(i) = padm(i)/weight
 50      continue
         alpha = alphabar/weight
         sig = sigbar/weight
         ngood = ngbar/weight
      endif
      else
c     Use historic MAXENT
c     First set cflag=1 so that update does not change alpha.
      cflag=1
c     now find alphamin and alphamax
      alphamin=0.0
 100  call update
      write(string,*) 'alpha= ',alpha
      if(aim.lt.aimaim)then
	alphamin=alpha
	alpha=2.0*alpha
	goto 100
      else
	alphamax=alpha
      end if
c     Now that we have an acceptable range for alpha, we must iterate
c     to find alpha(aimaim).
      alpha=0.5*(alphamin+alphamax)
      do 110 i=1,niter
	call update
        write(string,*) 'alpha= ',alpha
        if(abs(aim-aimaim)/(aim+aimaim).lt.0.0001) goto 120
	if(aim.gt.aimaim)then
	  call dchop(alpha,alphamin,alphamax)
   	else
	  call uchop(alpha,alphamin,alphamax)
	end if
 110  continue
 120  continue
      end if

c     Next calculate the error 
      if(aflag.eq.0.or.aflag.eq.2.or.aflag.eq.4) then
	   call errcalc
	   if(iasmc.eq.1) call Av_spec
	end if
      
c     Finally rewrite the image 
      call output
      
      stop
      end


      subroutine mxva(a,ira,ica,b,irb,c,irc,na,nb)
c     This should emulate the SCILIB routine MXVA.  It should form the
c     matrix-vector product C=AB, where NA is the number of rows in A and the
c     number of elements in C, NB is the number of elements in B and columns
c     in A.  To get between two elements in the same row but different columns,
c     add IC. To get between two elements in the same column but different
c     rows, add IR. 
      
      implicit double precision (a-h,o-z)
      double precision a(*),b(*),c(*)
      
      iaddc=1
      do 100 irow=1,na
         iadda=1+irow*ira-ira
         iaddb=1
         xxx=0.
         do 50 icol=1,nb
            xxx=xxx+a(iadda)*b(iaddb)
            iadda=iadda+ica
            iaddb=iaddb+irb
 50      continue
         c(iaddc)=xxx
         iaddc=iaddc+irc
 100  continue
      
      return
      end

