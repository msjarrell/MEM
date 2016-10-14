// This code produces a toy problem to test the bryan code for the 
// analytic continuation problem defined by
//
//          /    A(w) exp(-tau*w)
// G(tau) = | dw -----------------
//          /     1 + exp(-beta*w)
//
// or in matrix form
//
// G=KA
// 
// Generated are data with relative error fixed by sigma, the kernel,
// and default model.  
//
// The model is assmed to be a gaussian
//
// A is assumed to be a simple lorentzian located at w0, with width 
// gamma and normalized to one.
//
// to compile gcc Toy_Generator.c -lm
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//
// some fundamental parameters
//
#define PI 3.1415927

double rang(void);

int main()
{
	int iseed,i,j,nf,ntau;
	double r1,r2,r3,w0,wmin,dw,beta,gamma,dtau,sigma,modelWidth;
	
// declare some pointers
	FILE *ModelFile, *DataFile, *InitFile, *ToyfFile;
	double *w,*kernel,*model,*initimage,*A,*G,*tau;
	
	ModelFile= fopen("model","w");
	DataFile= fopen("kerneldata","w");
	InitFile= fopen("initimage","w");
	ToyfFile= fopen("toyA","w");

	printf("Enter a random number seed (int): ");
	scanf("%d",&iseed);
	srand(iseed); // seed the random number generator
	printf("Enter w0 the mean for the toy: ");
	scanf("%lf",&w0);
	printf("Enter gamma, the toy width: ");
	scanf("%lf",&gamma);
	printf("Enter the dw, the frequency step: ");
	scanf("%lf",&dw);
	printf("Enter the nf, the number pixels in A(odd int): ");
	scanf("%d",&nf);
	printf("Enter beta, the inverse temperature: ");
	scanf("%lf",&beta);
	printf("Enter ntau, the number of time steps: ");
	scanf("%d",&ntau);
	printf("Enter sigma, the absolute error of data: ");
	scanf("%lf",&sigma);
	printf("Enter modelWidth, the width of the gaussain model: ");
	scanf("%lf",&modelWidth);
	
	dtau=beta/(double)ntau;
	wmin=-dw*(nf/2);
	printf("dtau= %lf \n",dtau);
	printf("wmin= %lf \n",wmin);
	
// allocate some memory

	w = 		(double *) malloc(sizeof(double)*nf);
	kernel = 	(double *) malloc(sizeof(double)*nf*ntau);
	model = 	(double *) malloc(sizeof(double)*nf);
	initimage = 	(double *) malloc(sizeof(double)*nf);
	A = 		(double *) malloc(sizeof(double)*nf);
	G = 		(double *) malloc(sizeof(double)*ntau);
	tau = 	(double *) malloc(sizeof(double)*ntau);
	
//calculate the Matsubara time grid
	for(j=0;j<ntau;j++)
	{
		tau[j]=dtau*(double)j;
	}
	
// Calculate the Lorentzian toy spectrum
	for (i=0;i<nf;i++)
	{
		w[i]=i*dw+wmin;
		A[i]=dw*(gamma/PI)/(pow((w[i]-w0),2)+gamma*gamma);
		fprintf(ToyfFile,"%lf  %lf\n",w[i],A[i]/dw);
	}
	
	
// calculate the kernel
	for(i=0;i<nf;i++)
	{
		for(j=0;j<ntau;j++)
		{
		if(w[i]>0.0)
		{
			kernel[i+j*nf]=exp(-tau[j]*w[i])/(1.0+exp(-beta*w[i]));
		}
		else
		{
			kernel[i+j*nf]=exp((beta-tau[j])*w[i])/(1.0+exp(beta*w[i]));
		}
		}
	}
	
// calculate the pure (eror free) data from f and the kernel
	for(j=0;j<ntau;j++)
	{
		G[j]=0.0;
		for(i=0;i<nf;i++)
		{
			G[j]+=A[i]*kernel[i+j*nf];
		}
	}

//
// Now start writing thigs to the files
//
// first the header
	i=ntau; 
	j=nf;
	fprintf(DataFile,"%d  %d  1\n",i,j);
//
// now the data and its error
//
	for(j=0;j<ntau;j++)
	{
		r3=rang();
		r1=G[j] + sigma*r3;
		r2=sigma;
		fprintf(DataFile,"%lf  %lf \n",r1,r2);
//		printf("gaussian rn %lf\n",r3);
	}
//
//now write the kernel
//
	for(j=0;j<ntau;j++)
	{
		for(i=0;i<nf;i++)
		{
			fprintf(DataFile,"%lf\n",kernel[i+j*nf]);
		}
	}

//
//now form and write the model
//
	r1=0.0;
	for (i=0;i<100;i++)
	{
		r1+=pow(rang(),2);
	}
	r1=r1/100.0;
	printf("test of gaussian RNG=%lf\n",r1);
	
	for(i=0;i<nf;i++)
	{
		model[i]=(1.0/(modelWidth*sqrt(PI)))*exp(-pow((w[i]/modelWidth),2));
		fprintf(ModelFile,"%lf  %lf  %lf\n",w[i],model[i],dw);
		fprintf(InitFile,"%lf  %lf\n ",w[i],model[i]);
	}

// free the allocated memory
	free(w 	);
	free(kernel );
	free(model 	);
	free(initimage );
	free(A );
	free(G );
	free(tau );
return 0;
} 

double rang(void)
{
// we assume that the random number seed is set (srand)
//this code generates a gaussianly distributed random number with mean zero
//and a standard deviation of one.  http://www.taygeta.com/random/gaussian.html
         double x1, x2, w, y1, y2;
 
         do {
                 x1 = 2.0 * ((double)rand()/(double)RAND_MAX )- 1.0;
                 x2 = 2.0 * ((double)rand()/(double)RAND_MAX ) - 1.0;
                 w = x1 * x1 + x2 * x2;
         } while ( w >= 1.0 );

         w = sqrt( (-2.0 * log( w ) ) / w );
         y1 = x1 * w;
         y2 = x2 * w;
	 return y1;
}
