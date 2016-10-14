c     This is the global definitions block for the MaxEnt
c     code.  All global variables must be declared here.

c       alpha       Lagrange multiplier which relates C and S
c       alpha1      The value of alpha from the previous iterate
c	alpha0      The value of alpha which maximizes the posterior
c	aflag       Flag that determines how the image is found
c       data        raw data to be deconvolved
c	post        current value (unnormalized ) of the posterior
c       du          increment in the solution
c       dw          frequency step
c       error       standard error bars associated data
c       f           image array
c	f_real	    real image (f/dw)
c	width	    The smoothing width
c       ftr         alphanew=ftr*alphanew + (1-ftr)*alphaold
c       lo          current value of L
c       model       default model values of f
c       mu          step-size cut-off parameter
c       ndata       number of data points
c	nbins	    number of bins of data used in the course-graining
c       ndmax       maximum value of ndata
c       nf          number of image pixels
c       nfmax       maximum value of nf
c       niter       total number of iterations of alpha
c       ns          number of nf values used 
c       precision   assumed numerical precision of the code
c       sig         error scaling parameter
c       so          current value of S
c       step        multiple of sum(f) used to define max step los
c       t           kernel which relates f and data (data=t*f)
c       test1       test constant as defined by alpha*S+G
c       testf       defines an acceptable value for test
c       w           image frequencies
c	iuflow	    set to 1 if the image falls below e^-60
c	iasmc	    if 1, then an average spectrum MC is done upon convergence
c	iprint	    print level	    
c	ixvgr	    if 1, then print out xvgr graphics directives

      implicit none
      integer ndmax, nfmax, npmax
      parameter (ndmax=200, nfmax=500, npmax=300)

      double precision t(ndmax,nfmax), data(ndmax), error(ndmax),
     .     w(nfmax), dw(nfmax), Tf(ndmax)
      double precision f(nfmax), f_real(nfmax), model(nfmax), u(nfmax),
     .     Umat(nfmax,nfmax), s(nfmax), v(ndmax,ndmax), 
     .     xm(nfmax,nfmax),mr(nfmax)
      double precision ftr, precision, step, testf, range
      integer iter, niter, ndata, nf, ns, aflag,   
     .         cflag, nregions, nbins, iuflow

      double precision RHS(nfmax), YIRHS(nfmax), du(nfmax), 
     .                 lambda(nfmax), mu
      double precision alpha, alpha0, sig, sigbar, alphabar, alpha1,
     .     weight, padm(npmax), fbar(nfmax), lo, so,alpt(npmax),
     .     post,sbar,dssbar,ngood,ngbar,alphainit


      common /qmcdata/ t, data, error,w, dw
      common /image/ f, f_real, model, u,Umat, s, v, xm, 
     .               mr,Tf,
     &               nregions 
      common /parms/ ftr, precision, step, testf, range, ndata, 
     . nf, ns, nbins, iuflow, niter
      common /contrl/ RHS, YIRHS, du, lambda, mu
      common /varibles/ alpha, alpha0, sig, sigbar, alphabar,
     &       weight, padm, fbar, lo, so, alpha1,
     .       alpt,post,sbar,dssbar, alphainit,
     .       ngood, ngbar, aflag, iter, cflag


	double precision pi
	integer optiono,optionm
	parameter(pi=3.1415927)
	common/probdep/ optiono,optionm

c       parameters necessary for "historic MAXENT" (aflag=5)
	double precision aim,aimaim,alphamax,alphamin
	common/hist/ aim,aimaim,alphamax,alphamin

c	parameters for the average spectrum monte carlo
	integer iasmc,nsweeps,nwarms,nruns,idummy,iprint,ixvgr
	double precision fave(nfmax),faves(nfmax),fsdev(nfmax),
     &       alphaave,alphaaves,alphasdev,delu,delalpha
	common/asmc/ alphaave,alphaaves,alphasdev,fave,faves,fsdev,
     &               delu,delalpha,
     &               iasmc,nsweeps,nwarms,nruns,idummy,iprint,ixvgr

c	parameters and arrays which need to be used in both update
c	and errcalc, or errcalc and output

      	double precision P(nfmax,nfmax), Z(nfmax), R(nfmax,nfmax)
        double precision Gvalue(300),Dgvalue(300),low(300),
     .                   high(300),width

	common/both/P,Z,R,Gvalue,Dgvalue,low,high,width
