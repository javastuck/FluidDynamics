/*compile with gcc -lm vib01.c */
/*uses file vars.txt to read in parameter values */

/** twod01.c is a 2d continuum simulation of   	**/
/** the vibrated granular layer experiment    	**/
/** using forward time, centered differencing	**/
/** on a uniform grid.				**/
/** it is basically the matlab program 		**/
/** "fullslip00.m" converted to c.		**/

/**threed03.c extends twod01.c into 3d			**/
/**threed04.c saves less frequently			**/
/**threed04.c allows you to load init conds from file	**/
/**threed05.c fixes smoothing				**/
/**threed06.c comments out smoothing			**/


/** WARNING: threed07.c actually has worse smoothing    **/
/** than threed06.c.  Please use threed06.c		**/

/**psoned01.c is threed07.c, except it			**/
/**saves out at regular time intervals			**/
/**uses var numsnaps as # of snapshots saved per period **/
/**also saves out and loads in f*t (time in units of    **/
/**period), instead of t				**/

/**psoned02.c is psoned01.c except it uses better 	**/
/**smoothing from threed06.c				**/

/**psoned03.c gets more info from the variable file, 	**/
/**making it easier to change gridsize and size of box	**/

/**psoned04.c improves initial conditions and makes it  **/
/**easier to change the depth of the layer		**/

/**shocks01.c is the same as psoned04.c except writes	**/
/** data in binary form rather than in ascii		**/

/**bcfix01.c is the same as shocks01.c			**/

/** depthfix01.c modifies the depth to agree with the   **/
/** depth at random close packing.  This is the version **/
/** used for the papers in my dissertation              **/

/** memfix01.c is an attempt to fix a memory leak, but  **/
/** should be operationally identical to depthfix01.c   **/
/** memfix02.c is a further iteration, and is the last  **/
/** stable code from UT                                 **/

/**memfix03.c tries to free the write buffer in a       **/
/**attempt to fix memory issues.                        **/

/**memfix04a.c just changes variable names              **/

/**memfix04b.c explicitly calculates stress tensor      **/
/**and heat flux terms.  Should be equivalent to old    **/
/**versions but slightly different order of calc.       **/

/**student01.c is the same as memfix04b.c               **/
/**student01.c quits if density changes too much        **/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <signal.h>  
#include <omp.h>
#include <time.h>
/*#include "dthree.h"*/
/*#include "dgasdevrand.h"*/
/*#include "dm.h"*/

#define NUM_THREADS 4

/**The following functions create ways to access memory      **/
/**addresses as arrays rather than doing pointer arithmetic  **/
/**similar to numerical recipes in c nrutil but modified     **/
/**for three dimensional arrays and double precision         **/
void nrerror(error_text)
char error_text[];
{
	void exit();

	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}



double *vector(nl,nh)
int nl,nh;
{
	double *v;

	v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl;
}

int *ivector(nl,nh)
int nl,nh;
{
	int *v;

	v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl;
}

double *dvector(nl,nh)
int nl,nh;
{
	double *v;

	v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl;
}



double **matrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
	int i;
	double **m;

	m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
		if (!m[i]) nrerror("allocation failure 2 in matrix()");
		m[i] -= ncl;
	}
	return m;
}

double **dmatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
	int i;
	double **m;

	m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	if (!m) nrerror("allocation failure 1 in dmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
		if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}


double ***dthreearray(nrl,nrh,ncl,nch,nthreel,nthreeh)
     int nrl,nrh,ncl,nch,nthreel,nthreeh;
{
  int i,j;
  double ***a;


  a=(double ***) malloc((unsigned) (nrh-nrl+1)*sizeof(double**));
  if (!a) nrerror("allocation failure 1 in dmatrix()");
  a -= nrl;

  for(i=nrl;i<=nrh;i++) {
    a[i]=(double **) malloc((unsigned) (nch-ncl+1)*sizeof(double*));
    if (!a[i]) nrerror("allocation failure 2 in dmatrix()");
    a[i] -= ncl;
  }
	
  for(i=nrl;i<=nrh;i++) {
    for(j=ncl;j<=nch;j++){
      a[i][j]=(double *) malloc((unsigned) (nthreeh-nthreel+1)*sizeof(double));
      if (!a[i][j]) nrerror("allocation failure 3 in dmatrix()");
      a[i][j] -= ncl;
    }
  }
  return a;
}


int **imatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
	int i,**m;

	m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
	if (!m) nrerror("allocation failure 1 in imatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
		if (!m[i]) nrerror("allocation failure 2 in imatrix()");
		m[i] -= ncl;
	}
	return m;
}



double **submatrix(a,oldrl,oldrh,oldcl,oldch,newrl,newcl)
double **a;
int oldrl,oldrh,oldcl,oldch,newrl,newcl;
{
	int i,j;
	double **m;

	m=(double **) malloc((unsigned) (oldrh-oldrl+1)*sizeof(double*));
	if (!m) nrerror("allocation failure in submatrix()");
	m -= newrl;

	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+oldcl-newcl;

	return m;
}



void free_vector(v,nl,nh)
double *v;
int nl,nh;
{
	free((char*) (v+nl));
}

void free_ivector(v,nl,nh)
int *v,nl,nh;
{
	free((char*) (v+nl));
}


void free_dvector(v,nl,nh)
double *v;
int nl,nh;
{
	free((char*) (v+nl));
}



void free_matrix(m,nrl,nrh,ncl,nch)
double **m;
int nrl,nrh,ncl,nch;
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_dmatrix(m,nrl,nrh,ncl,nch)
double **m;
int nrl,nrh,ncl,nch;
{
	int i;
	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_dthreearray(m,nrl,nrh,ncl,nch,nthreel,nthreeh)
     double ***m;
     int nrl,nrh,ncl,nch,nthreel,nthreeh;
{
  int i,j;
  for(i=nrh;i>=nrl;i--){
    for(j=nch;j>=ncl;j--){
      free((char*) (m[i][j]+nthreel));
    }
  }
  free((char*) (m+nrl));
}


void free_imatrix(m,nrl,nrh,ncl,nch)
int **m;
int nrl,nrh,ncl,nch;
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}



void free_submatrix(b,nrl,nrh,ncl,nch)
double **b;
int nrl,nrh,ncl,nch;
{
	free((char*) (b+nrl));
}



double **convert_matrix(a,nrl,nrh,ncl,nch)
double *a;
int nrl,nrh,ncl,nch;
{
	int i,j,nrow,ncol;
	double **m;

	nrow=nrh-nrl+1;
	ncol=nch-ncl+1;
	m = (double **) malloc((unsigned) (nrow)*sizeof(double*));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m -= nrl;
	for(i=0,j=nrl;i<=nrow-1;i++,j++) m[j]=a+ncol*i-ncl;
	return m;
}



void free_convert_matrix(b,nrl,nrh,ncl,nch)
double **b;
int nrl,nrh,ncl,nch;
{
	free((char*) (b+nrl));
}


/*  This function creates random numbers with a normal distribution */
/*  And unit variance, just like matlab's randn		            */

#include <math.h>

double gasdevrand()
    

{
  static int iset = 0;
  static double gset;
  double fac, r, v1,v2;

  if(iset == 0){
    do{
      v1=2.0*rand()/(1.0*RAND_MAX+1.0)-1.0;
      v2=2.0*rand()/(1.0*RAND_MAX+1.0)-1.0;
      r=v1*v1+v2*v2;
    }while(r>=1.0);
    fac=sqrt(-2.0*log(r)/r);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  }else{
    iset=0;
    return gset;
  }
}



static int checkforsignaltoend;


/** catch ctl-c and make it clear memory before quitting **/

/* here is the signal handler */
void catch_int(int sig_num)

{
    /* re-set the signal handler again to catch_int, for next time */
    signal(SIGINT, catch_int);
    /* set the check to zero, telling while to exit */
    checkforsignaltoend=0;

    printf("quitting, please wait. . . \n");
}


/***********     Start Main Program ******************/

int main(int argc, char **argv) {
  /** Input/Output **/
  omp_set_num_threads(4);
  time_t start, end;
  start = clock();
  FILE *output1;
  FILE *input;
  char file_name[255],base[255];
  
  /** Flow variables **/

  double ***ph,***T,***v,***w,***u,***P,***mu,***K;
  double nu,***G,***lam,***gam;
  
  /**derivatives **/

  double  ***phy,***vy,***wy,***uy,***Ty,***Py;
  double  ***phz,***vz,***wz,***uz,***Tz,***Pz;
  double  ***phx,***vx,***wx,***ux,***Tx,***Px;
  double  ***divU, ***stressxx, ***stressyy, ***stresszz;
  double  ***stressyz,***stressxy,***stressxz,***stressyzDz,***stressxyDx;
  double  ***stressxyDy,***stressxzDz,***stressyzDy,***stressxzDx;
  double  ***stressxxDx,***stressyyDy,***stresszzDz,tone;
  double  ***divflux,***heatfluxx,***heatfluxy,***heatfluxz;

  
  /** Parameters **/

  double GAM,fstar,e,mass,f,A,To,pho,numax,gr,pi;

  /** Grid variables **/

  double ggx,ggy,h,*x,*y,*z,dx,dy,dz,sizex,sizey;
  int Nx,Ny,Nz;
    
  /** Iteration variables **/

  double dt,C,s,tau,tmp,cmass,temp,ttmp,t,time,pl,maxdt;
  double mindt,dtstep, ev, ed,dtone,gamo,onepluse;
  double ***dph,***dw,***dv,***du,***dT;
  int	n,nstart,Nf,i,j,k,m,savecount,numsnaps,tsnap,numsmooth,tsmooth;

  /** Smoothing and superviscosity variables **/

  double supervis, ***fv,***fw,***fu,***fT,***fph;
  int 	presmooth, preload;

  /** fix mass to correspond with correct layer depth **/

  double masscorr,depth;

  /** uselesss junk variable to make malloc really use free **/

  char *junkvariable;
junkvariable=(char *) malloc(1);
 free(junkvariable);




  /****** DEFINE PI *****************/
  
    pi=1.0*acos(-1.0);

    /** Get the appropriate parameter values from file vars.txt**/
    
  input=fopen("vars.txt","r");
  
  fscanf(input,"%lf",&GAM);
  printf("Gamma is:");
  printf("%lf\n",GAM);

  fscanf(input,"%lf",&fstar);
  printf("fstar is:");
  printf("%lf\n",fstar);

  fscanf(input,"%lf",&e);
  printf("coefficient of restitution is: ");
  printf("%lf\n",e);

  fscanf(input,"%lf",&depth);
  printf("depth of the layer is: ");
  printf("%lf\n",depth);

  masscorr=0.9;
  mass=depth/masscorr;

  fscanf(input,"%d",&Nx);
  printf("number of gridpoints in the x direction: ");
  printf("%d\n",Nx);

  fscanf(input,"%d",&Ny);
  printf("number of gridpoints in the y direction: ");
  printf("%d\n",Ny);

  fscanf(input,"%d",&Nz);
  printf("number of gridpoints in the z direction: ");
  printf("%d\n",Nz);

  fscanf(input,"%lf",&h);
  printf("height of the box (z-direction): ");
  printf("%lf\n",h);

  fscanf(input,"%lf",&sizex);
  printf("size of box in x-direction: ");
  printf("%lf\n",sizex);

  ggx=sizex/h; /*aspect ratio*/
  
  fscanf(input,"%lf",&sizey);
  printf("size of box in y-direction: ");
  printf("%lf\n",sizey);

  ggy=sizey/h; /*aspect ratio*/
  
  fscanf(input,"%d",&numsnaps);
  printf("number of snapshots to save per period: ");
  printf("%d\n",numsnaps);
  
  /** Get save filename from file **/

  fscanf(input,"%s",base);
  printf("variables are saved in file:");
  printf("%s\n",base);
  
  fscanf(input,"%d",&savecount);
  printf("start iterations at file number:");
  printf("%d\n",savecount);

  fscanf(input, "%d",&preload);
  printf("preload?");
  printf("%d\n",preload);
 fclose(input);


  /** Define size of grid in y and z directions **/

  /*
  Nz=201;
  Ny=7;
  Nx=7;

  */
  /*savecount=0;*/

  supervis=200.0;		/*constant used in superviscous smoothing*/
  presmooth=0;		/*how many times to presmooth*/
  numsmooth=0;			/*haven't smoothed yet at this point*/
  
  
  /*define the largest and smallest timesteps that will be permitted*/
  /*also define how much bigger the timestep can get each  step */
  maxdt=.002;
  mindt=.00001;
  dtstep=.000005;

  /*set dt to mindt */
  dt=mindt;

  /*define variables to check timestep size */ 

  C=.02;
  s=1.5;
  tau=.3;

  /** decide how many timesteps to do. 	*/
  /* starts at n=nstart+1, goes to n=Nf  	*/

  //Nf=10000000; number of iterations
  Nf = 500;
  
  nstart=1;

  /** define properties of the system **/
  
  /*  ggx=3.0*Nx/59.0*201.0/(1.0*Nz);*/	/*aspect ratio of cell */
  /* ggy=3.0*Ny/59.0*201.0/(1.0*Nz);*/	/*aspect ratio of cell */
  /*h=40.0;*/
				/*height of cell */
  numax=.65;			/*max volume fraction  */
  
  /*mass=6.0;*/			/*mass=integrate(ph(y,z)dy dz/L), */
				/*  L=gg*h.  mass=bad name*/

  pho=.9*6.0/pi*numax;		/* base value for initial conditions*/
  To=.01;


  f=fstar/sqrt(mass);		/* frequency of oscillation */ 
  A=GAM/((2.0*pi*f)*(2.0*pi*f));	/*amplitude of oscillation */

  ed=12.0*(1.0-e*e)/sqrt(pi);
  gamo=ed;

  gr=1.0;				/*gravitational constant */

  onepluse=2.0;
  
  /** define the matrices for the flow variables **/
  
  ph = dthreearray(1,Nx,1,Ny,1,Nz);
  T=dthreearray(1,Nx,1,Ny,1,Nz);
  u=dthreearray(1,Nx,1,Ny,1,Nz);
  v=dthreearray(1,Nx,1,Ny,1,Nz);
  w=dthreearray(1,Nx,1,Ny,1,Nz);
  P = dthreearray(1,Nx,1,Ny,1,Nz);
  G = dthreearray(1,Nx,1,Ny,1,Nz);

  fph = dthreearray(1,Nx,1,Ny,1,Nz);
  fT=dthreearray(1,Nx,1,Ny,1,Nz);
  fu=dthreearray(1,Nx,1,Ny,1,Nz);
  fv=dthreearray(1,Nx,1,Ny,1,Nz);
  fw=dthreearray(1,Nx,1,Ny,1,Nz);

  dph=dthreearray(1,Nx,1,Ny,1,Nz);
  du=dthreearray(1,Nx,1,Ny,1,Nz);
  dv=dthreearray(1,Nx,1,Ny,1,Nz);
  dw=dthreearray(1,Nx,1,Ny,1,Nz);
  dT=dthreearray(1,Nx,1,Ny,1,Nz);

  
  phx=dthreearray(1,Nx,1,Ny,1,Nz);
  ux=dthreearray(1,Nx,1,Ny,1,Nz);
  vx=dthreearray(1,Nx,1,Ny,1,Nz);
  wx=dthreearray(1,Nx,1,Ny,1,Nz);
  Tx=dthreearray(1,Nx,1,Ny,1,Nz);
  Px=dthreearray(1,Nx,1,Ny,1,Nz);
  phy=dthreearray(1,Nx,1,Ny,1,Nz);
  uy=dthreearray(1,Nx,1,Ny,1,Nz);
  vy=dthreearray(1,Nx,1,Ny,1,Nz);
  wy=dthreearray(1,Nx,1,Ny,1,Nz);
  Ty=dthreearray(1,Nx,1,Ny,1,Nz);
  Py=dthreearray(1,Nx,1,Ny,1,Nz);
  phz=dthreearray(1,Nx,1,Ny,1,Nz);
  uz=dthreearray(1,Nx,1,Ny,1,Nz);
  vz=dthreearray(1,Nx,1,Ny,1,Nz);
  wz=dthreearray(1,Nx,1,Ny,1,Nz);
  Tz=dthreearray(1,Nx,1,Ny,1,Nz);
  Pz = dthreearray(1,Nx,1,Ny,1,Nz);

  divU = dthreearray(1,Nx,1,Ny,1,Nz);
  mu=dthreearray(1,Nx,1,Ny,1,Nz);
  K=dthreearray(1,Nx,1,Ny,1,Nz);
  lam=dthreearray(1,Nx,1,Ny,1,Nz);
  gam=dthreearray(1,Nx,1,Ny,1,Nz);
  
  stressxx = dthreearray(1,Nx,1,Ny,1,Nz);
  stressyy = dthreearray(1,Nx,1,Ny,1,Nz);
  stresszz = dthreearray(1,Nx,1,Ny,1,Nz);
  stressyz = dthreearray(1,Nx,1,Ny,1,Nz);
  stressxy = dthreearray(1,Nx,1,Ny,1,Nz);
  stressxz = dthreearray(1,Nx,1,Ny,1,Nz);
  stressxxDx=dthreearray(1,Nx,1,Ny,1,Nz);
  stressyyDy=dthreearray(1,Nx,1,Ny,1,Nz);
  stresszzDz=dthreearray(1,Nx,1,Ny,1,Nz);
  stressxyDy=dthreearray(1,Nx,1,Ny,1,Nz);
  stressxzDz=dthreearray(1,Nx,1,Ny,1,Nz);
  stressyzDy=dthreearray(1,Nx,1,Ny,1,Nz);
  stressxzDx=dthreearray(1,Nx,1,Ny,1,Nz);
  stressyzDz=dthreearray(1,Nx,1,Ny,1,Nz);
  stressxyDx=dthreearray(1,Nx,1,Ny,1,Nz);

  divflux=dthreearray(1,Nx,1,Ny,1,Nz);
  heatfluxx=dthreearray(1,Nx,1,Ny,1,Nz);
  heatfluxy=dthreearray(1,Nx,1,Ny,1,Nz);
  heatfluxz=dthreearray(1,Nx,1,Ny,1,Nz);


  printf("ok\n");




  /**initialize while loop to keep running **/
  
  checkforsignaltoend=1;

/* set the INT (Ctrl-C) signal handler to 'catch_int' */
signal(SIGINT, catch_int);

  /*define vectors representing height and horizontal distance*/

  x=dvector(1,Nx);
  y=dvector(1,Ny);
  z=dvector(1,Nz);

  
  dz=h/(1.0*(Nz-1));
  dy=ggy*h/(1.0*(Ny-1));
  dx=ggx*h/(1.0*(Nx-1));

//Change01 8/2/2016
for(k=1;k<=Nz;k++){
	z[k]=dz*(k-1);} 
for(j=1;j<=Ny;j++){
	y[j]=dy*(j-1);}
for(i=1;i<=Nx;i++){
	x[i]=dx*(i-1);}

  
//end change01





  /*define initial conditions */ 
  for(i=1;i<=Nx;i++){
	for(j=1;j<=Ny;j++){
	  /*
	temp=2.0*rand()/(1.0*RAND_MAX+1.0)-1.0;
	ttmp=2.0*rand()/(1.0*RAND_MAX+1.0)-1.0;
	  */
	
	  temp=gasdevrand();
	  ttmp=gasdevrand();
	
	  for(k=1;k<=Nz;k++){
	tmp=(.99+.01*rand()/(1.0*RAND_MAX+1.0))*(1.0-tanh(z[k]-(1.07+mass+.8*temp)))/2.0;	
	ph[i][j][k]=tmp*(.99+.01*rand()/(1.0*RAND_MAX+1.0))*(1.0-tanh((1.0+.8*ttmp)-z[k]))/2.0+ .000001;
	

	/*	tmp=(1.0-tanh(z[k]-(1.07+mass+.8*temp)))/2.0;
	ph[i][j][k]=tmp*(1.0-tanh((1.0+.8*ttmp)-z[k]))/2.0;
	*/
	  }
	}
  }


  
  /*Now set the average layer depth to mass -- correct if neccessary	*/
//start change02
  
    cmass=0;

    for(i=1;i<=Nx;i++){  
      for(j=1;j<=Ny;j++){	
	for(k=1;k<=Nz;k++){
	  cmass = cmass + ph[i][j][k];
	}
      }
    }
    printf("%lf %lf %lf\n",dx,dz, dy);
    cmass=cmass*h/(1.0*Nx*Ny*Nz);
    printf("cmass %lf\n",cmass);
      
    temp=0.0;

    for(i=1;i<=Nx;i++){
      for(j=1;j<=Ny;j++){
	for(k=1;k<=Nz;k++){
	  fph[i][j][k] = ph[i][j][k]*exp(-100.0*(ph[i][j][k]-numax*6.0/pi/2.0)*(ph[i][j][k]-numax*6.0/pi/2.0));
	  temp=temp + fph[i][j][k];
	}
      }
    }
    for(i=1;i<=Nx;i++){
      for(j=1;j<=Ny;j++){
	for(k=1;k<=Nz;k++){
	  ph[i][j][k]=ph[i][j][k] + fph[i][j][k]*(mass-cmass)/
	    (h*temp/(1.0*Nx*Ny*Nz));
	  if(ph[i][j][k]<0){
	    printf("negph \n");}
	}
      }
    }
//end change02 

//start change03
  /*use values found above for density, and initialize u,v,w,t*/
    for(i=1;i<=Nx;i++){
      for(j=1;j<=Ny;j++){
	for(k=1;k<=Nz;k++){
	  u[i][j][k]=0.0;
	  v[i][j][k]=0.0;
	  w[i][j][k]=-.3;
	  T[i][j][k]=2.0*To;     
	}
      }
    }
//end change03
  /*initial time is one quarter period before t=0*/
  
  t=((acos(0.0)+pi)-2.0*pi)/2.0/pi/f;
  time=((acos(0.0)+pi)-2.0*pi)/2.0/pi;

  
  /*find values for G, P, etc (not really necessary */
  /*
  for(j=1;j<=Ny;j++){
    for(k=1;k<=Nz;k++){
    nu=pi/6.0*ph[j][k];
    G[j][k]=nu/(1.0-pow(nu/numax,4.0*numax/3.0));
    P[j][k]=(1.0+2.0*onepluse*G[j][k])*T[j][k]*ph[j][k];
    tone=ph[j][k]*sqrt(T[j][k]);
    lam[j][k]=8.0/3.0/sqrt(pi)*tone*G[j][k];
    mu[j][k]=sqrt(pi)/6.0*tone*(1.0+4.0/5.0*(1.0+12.0/pi)*G[j][k]);
    K[j][k]=15.0*sqrt(pi)/16.0*tone*(1.0+6.0/5.0*(1.0+32.0/9.0/pi)*G[j][k]);
    ev=1.0-(1.0-e)*pow((sqrt(T[j][k])/1.0),.75);
    if(T[j][k]>1.0){
      ev=e;}
    ed=12.0*(1.0-ev*ev)/sqrt(pi);
    gam[j][k]=ed*G[j][k]*T[j][k]*tone;
    }
  }
  */
  /***** matlab code finds derivatives here *********/

  
  /*done initializing */

  /**********load in variables from old file if so desired **********/
  if(preload==1){
    
     sprintf(file_name,"%s%04d.dat",base,savecount); 


     input = fopen(file_name,"rb");

 
     fread(&tmp,sizeof(double),1,input);  	/* fstar from file 	*/
     fread(&tmp,sizeof(double),1,input);	/*GAM from file		*/
     
     /*    fscanf(input,"%lf %lf \n",&tmp, &tmp);*/
 
    printf("%s\n", file_name);

    fread(&Nx,sizeof(int),1,input);
    fread(&Ny,sizeof(int),1,input);
    fread(&Nz,sizeof(int),1,input);
    
    /*    fscanf(input,"%d %d %d",&Nx, &Ny, &Nz);*/
      printf("%d %d %d\n",Nx, Ny, Nz);
      
      fread(&n,sizeof(int),1,input);
      fread(&time,sizeof(double),1,input);
      fread(&numsmooth,sizeof(int),1,input);
      
      /*fscanf(input,"%d %lf %d\n",&n,&t,&numsmooth);*/

    t=time/f; /*what is really saved in this position is t*f, not t*/

    fread(&sizex,sizeof(double),1,input);
    fread(&sizey,sizeof(double),1,input);
    fread(&h,sizeof(double),1,input);
    fread(&depth,sizeof(double),1,input);

    mass=depth/masscorr;

    ggx=sizex/h;
    ggy=sizey/h;
    
    /*    fscanf(input," %lf %lf %lf %lf\n", &ggx, &ggy,&h,&mass);*/

    /*
    fread(ph,sizeof(double),Nx*Ny*Nz,input);
    fread(u,sizeof(double),Nx*Ny*Nz,input);
    fread(v,sizeof(double),Nx*Ny*Nz,input);
    fread(w,sizeof(double),Nx*Ny*Nz,input);
    fread(T,sizeof(double),Nx*Ny*Nz,input);
    */
    

    
    for(i=1;i<=Nx;i++){
      for(j=1;j<=Ny;j++){
	for(k=1;k<=Nz;k++){
	  fread(&ph[i][j][k],sizeof(double),1,input);
	}
      }
    }    
    for(i=1;i<=Nx;i++){
      for(j=1;j<=Ny;j++){
	for(k=1;k<=Nz;k++){
	fread(&u[i][j][k],sizeof(double),1,input);
	}
      }
    }

   for(i=1;i<=Nx;i++){
      for(j=1;j<=Ny;j++){
        for(k=1;k<=Nz;k++){
          fread(&v[i][j][k],sizeof(double),1,input);
        }
      }
    }
    
    for(i=1;i<=Nx;i++){
      for(j=1;j<=Ny;j++){
        for(k=1;k<=Nz;k++){
          fread(&w[i][j][k],sizeof(double),1,input);
        }
      }
    }

    for(i=1;i<=Nx;i++){    
      for(j=1;j<=Ny;j++){
        for(k=1;k<=Nz;k++){
          fread(&T[i][j][k],sizeof(double),1,input);
        }
      }
    }

    
    
    /*close data file*/
      
    fclose(input);
  }

  
  /*smooth out initial conditions from 1 to presmooth */
  
  for(n=1;n<=presmooth;n++){
    for(i=1;i<=Nx;i++){
      for(j=1;j<=Ny;j++){
	for(k=2;k<=Nz-1;k++){
	  fph[i][j][k]=(.5*ph[i][j][k-1] + ph[i][j][k] + .5*ph[i][j][k+1])/2.0;
	  fu[i][j][k]=(.5*u[i][j][k-1] + u[i][j][k] + .5*u[i][j][k+1])/2.0;
	  fv[i][j][k]=(.5*v[i][j][k-1] + v[i][j][k] + .5*v[i][j][k+1])/2.0;
	  fw[i][j][k]=(.5*w[i][j][k-1] + w[i][j][k] + .5*w[i][j][k+1])/2.0;
	  fT[i][j][k]=(.5*T[i][j][k-1] + T[i][j][k] + .5*T[i][j][k+1])/2.0;
	}
	
	fph[i][j][1]=(.75*ph[i][j][1] + .25*ph[i][j][2]);
	fu[i][j][1]=(.75*u[i][j][1] + .25*u[i][j][2]);
	fv[i][j][1]=(.75*v[i][j][1] + .25*v[i][j][2]);

		/*fv[j][1]=(.75*v[j][1] + .25*v[j][2]);*/
	fw[i][j][1]=(.75*w[i][j][1] + .25*w[i][j][2]);
	fT[i][j][1]=(.75*T[i][j][1] + .25*T[i][j][2]);

	fph[i][j][Nz]=(.75*ph[i][j][Nz] + .25*ph[i][j][Nz-1]);
	fu[i][j][Nz]=(.75*u[i][j][Nz] + .25*u[i][j][Nz-1]);
	fv[i][j][Nz]=(.75*v[i][j][Nz] + .25*v[i][j][Nz-1]);
	/*fv[j][Nz]=(.75*v[j][Nz] + .25*v[j][Nz-1]);*/
	fw[i][j][Nz]=(.75*w[i][j][Nz] + .25*w[i][j][Nz-1]);
	fT[i][j][Nz]=(.75*T[i][j][Nz] + .25*T[i][j][Nz-1]);
	
      }
    }

    for(i=1;i<=Nx;i++){
      for(k=1;k<=Nz;k++){
	for(j=2;j<Ny;j++){
	  fph[i][j][k]=(.1*fph[i][j-1][k] + fph[i][j][k]
		       + .1*fph[i][j+1][k])/1.2;
	  fu[i][j][k]=(.1*fu[i][j-1][k] + fu[i][j][k] + .1*fu[i][j+1][k])/1.2;
	  fv[i][j][k]=(.1*fv[i][j-1][k] + fv[i][j][k] + .1*fv[i][j+1][k])/1.2;
	  fw[i][j][k]=(.1*fw[i][j-1][k] + fw[i][j][k] + .1*fw[i][j+1][k])/1.2;
	  fT[i][j][k]=(.1*fT[i][j-1][k] + fT[i][j][k] + .1*fT[i][j+1][k])/1.2;
	}
	fph[i][1][k]=(.1*fph[i][Ny][k] + fph[i][1][k] + .1*fph[i][2][k])/1.2;
	fu[i][1][k]=(.1*fu[i][Ny][k] + fu[i][1][k] + .1*fu[i][2][k])/1.2;
	fv[i][1][k]=(.1*fv[i][Ny][k] + fv[i][1][k] + .1*fv[i][2][k])/1.2;
	fw[i][1][k]=(.1*fw[i][Ny][k] + fw[i][1][k] + .1*fw[i][2][k])/1.2;
	fT[i][1][k]=(.1*fT[i][Ny][k] + fT[i][1][k] + .1*fT[i][2][k])/1.2;

	fph[i][Ny][k]=(.1*fph[i][Ny-1][k] + fph[i][Ny][k]
		      + .1*fph[i][1][k])/1.2;
	fu[i][Ny][k]=(.1*fu[i][Ny-1][k] + fu[i][Ny][k] + .1*fu[i][1][k])/1.2;
	fv[i][Ny][k]=(.1*fv[i][Ny-1][k] + fv[i][Ny][k] + .1*fv[i][1][k])/1.2;
	fw[i][Ny][k]=(.1*fw[i][Ny-1][k] + fw[i][Ny][k] + .1*fw[i][1][k])/1.2;
	fT[i][Ny][k]=(.1*fT[i][Ny-1][k] + fT[i][Ny][k] + .1*fT[i][1][k])/1.2;
      }
    }
  
    for(j=1;j<=Ny;j++){
      for(k=1;k<=Nz;k++){
	for(i=2;i<Nx;i++){
	  ph[i][j][k]=(.1*fph[i-1][j][k] + fph[i][j][k] +
		       .1*fph[i+1][j][k])/1.2;
	  u[i][j][k]=(.1*fu[i-1][j][k] + fu[i][j][k] + .1*fu[i+1][j][k])/1.2;
	  v[i][j][k]=(.1*fv[i-1][j][k] + fv[i][j][k] + .1*fv[i+1][j][k])/1.2;
	  w[i][j][k]=(.1*fw[i-1][j][k] + fw[i][j][k] + .1*fw[i+1][j][k])/1.2;
	  T[i][j][k]=(.1*fT[i-1][j][k] + fT[i][j][k] + .1*fT[i+1][j][k])/1.2;
	}
	ph[1][j][k]=(.1*fph[Nx][j][k] + fph[1][j][k] + .1*fph[2][j][k])/1.2;
	u[1][j][k]=(.1*fu[Nx][j][k] + fu[1][j][k] + .1*fu[2][j][k])/1.2;
	v[1][j][k]=(.1*fv[Nx][j][k] + fv[1][j][k] + .1*fv[2][j][k])/1.2;
	w[1][j][k]=(.1*fw[Nx][j][k] + fw[1][j][k] + .1*fw[2][j][k])/1.2;
	T[1][j][k]=(.1*fT[Nx][j][k] + fT[1][j][k] + .1*fT[2][j][k])/1.2;

	ph[Nx][j][k]=(.1*fph[Nx-1][j][k] + fph[Nx][j][k] +
		      .1*fph[1][j][k])/1.2;
	u[Nx][j][k]=(.1*fu[Nx-1][j][k] + fu[Nx][j][k] + .1*fu[1][j][k])/1.2;
	v[Nx][j][k]=(.1*fv[Nx-1][j][k] + fv[Nx][j][k] + .1*fv[1][j][k])/1.2;
	w[Nx][j][k]=(.1*fw[Nx-1][j][k] + fw[Nx][j][k] + .1*fw[1][j][k])/1.2;
	T[Nx][j][k]=(.1*fT[Nx-1][j][k] + fT[Nx][j][k] + .1*fT[1][j][k])/1.2;
      }
    }

  
  }


  numsmooth += presmooth;

  if(preload==0){
    tsnap=floor(time*numsnaps-1);} /*set tsnap so there is a save on 1st step*/
  else{
    tsnap=floor(time*numsnaps);} /*set tsnap so there is no save on 1st step*/
  /* now start timestepping */

  n=nstart+1;
  //while(checkforsignaltoend==1){
  for(n=nstart;n<=Nf;n++){
    //n++; 
      /* find plate position */
    
      pl=A*(1.0-cos(2.0*pi*time));
    
      /*Use superviscosity for low density regions	*/
      /*Smooth to prevent blowups in regions where */
      /*there are physically no particles		*/

      for(i=1;i<=Nx;i++){	
	for(j=1;j<=Ny;j++){
	  for(k=2;k<Nz;k++){  
	    fu[i][j][k]=(u[i][j][k-1] + u[i][j][k] + u[i][j][k+1])/3.0;
	    fv[i][j][k]=(v[i][j][k-1] + v[i][j][k] + v[i][j][k+1])/3.0;
	    fw[i][j][k]=(w[i][j][k-1] + w[i][j][k] + w[i][j][k+1])/3.0;
	    fT[i][j][k]=(T[i][j][k-1] + T[i][j][k] + T[i][j][k+1])/3.0;
	  }
	  
	  fu[i][j][1]=(2.0*u[i][j][1] + u[i][j][2])/3.0;
	  fv[i][j][1]=(2.0*v[i][j][1] + v[i][j][2])/3.0;
	  fw[i][j][1]=(2.0*w[i][j][1] + w[i][j][2])/3.0;
	  fT[i][j][1]=(2.0*T[i][j][1] + T[i][j][2])/3.0;

	  fu[i][j][Nz]=(2.0*u[i][j][Nz] + u[i][j][Nz-1])/3.0;
	  fv[i][j][Nz]=(2.0*v[i][j][Nz] + v[i][j][Nz-1])/3.0;
	  fw[i][j][Nz]=(2.0*w[i][j][Nz] + w[i][j][Nz-1])/3.0;
	  fT[i][j][Nz]=(2.0*T[i][j][Nz] + T[i][j][Nz-1])/3.0;
	
	}
      }

      for(i=1;i<=Nx;i++){      
	for(k=1;k<=Nz;k++){
	  for(j=2;j<Ny;j++){
	    fu[i][j][k]=(fu[i][j-1][k] + fu[i][j][k] + fu[i][j+1][k])/3.0;
	    fv[i][j][k]=(fv[i][j-1][k] + fv[i][j][k] + fv[i][j+1][k])/3.0;
	    fw[i][j][k]=(fw[i][j-1][k] + fw[i][j][k] + fw[i][j+1][k])/3.0;
	    fT[i][j][k]=(fT[i][j-1][k] + fT[i][j][k] + fT[i][j+1][k])/3.0;
	  }
	  fu[i][1][k]=(fu[i][Ny][k] + fu[i][1][k] + fu[i][2][k])/3.0;
	  fv[i][1][k]=(fv[i][Ny][k] + fv[i][1][k] + fv[i][2][k])/3.0;
	  fw[i][1][k]=(fw[i][Ny][k] + fw[i][1][k] + fw[i][2][k])/3.0;
	  fT[i][1][k]=(fT[i][Ny][k] + fT[i][1][k] + fT[i][2][k])/3.0;

	  fu[i][Ny][k]=(fu[i][Ny-1][k] + fu[i][Ny][k] + fu[i][1][k])/3.0;
	  fv[i][Ny][k]=(fv[i][Ny-1][k] + fv[i][Ny][k] + fv[i][1][k])/3.0;
	  fw[i][Ny][k]=(fw[i][Ny-1][k] + fw[i][Ny][k] + fw[i][1][k])/3.0;
	  fT[i][Ny][k]=(fT[i][Ny-1][k] + fT[i][Ny][k] + fT[i][1][k])/3.0;
	}
      }

    for(k=1;k<=Nz;k++){
	for(j=1;j<=Ny;j++){
	  for(i=2;i<Nx;i++){
	    fu[i][j][k]=(fu[i-1][j][k] + fu[i][j][k] + fu[i+1][j][k])/3.0;
	    fv[i][j][k]=(fv[i-1][j][k] + fv[i][j][k] + fv[i+1][j][k])/3.0;
	    fw[i][j][k]=(fw[i-1][j][k] + fw[i][j][k] + fw[i+1][j][k])/3.0;
	    fT[i][j][k]=(fT[i-1][j][k] + fT[i][j][k] + fT[i+1][j][k])/3.0;
	  }

	  fu[1][j][k]=(fu[Nx][j][k] + fu[1][j][k] + fu[2][j][k])/3.0;
	  fv[1][j][k]=(fv[Nx][j][k] + fv[1][j][k] + fv[2][j][k])/3.0;
	  fw[1][j][k]=(fw[Nx][j][k] + fw[1][j][k] + fw[2][j][k])/3.0;
	  fT[1][j][k]=(fT[Nx][j][k] + fT[1][j][k] + fT[2][j][k])/3.0;

	  
	  fu[Nx][j][k]=(fu[Nx-1][j][k] + fu[Nx][j][k] + fu[1][j][k])/3.0;
	  fv[Nx][j][k]=(fv[Nx-1][j][k] + fv[Nx][j][k] + fv[1][j][k])/3.0;
	  fw[Nx][j][k]=(fw[Nx-1][j][k] + fw[Nx][j][k] + fw[1][j][k])/3.0;
	  fT[Nx][j][k]=(fT[Nx-1][j][k] + fT[Nx][j][k] + fT[1][j][k])/3.0;
	
	}
    }

     
      
    for(i=1;i<=Nx;i++){
      for(j=1;j<=Ny;j++){
	for(k=1;k<=Nz;k++){
	  temp=exp(-ph[i][j][k]*supervis);
	  u[i][j][k]=temp*fu[i][j][k] + (1.0-temp)*u[i][j][k];
	  v[i][j][k]=temp*fv[i][j][k] + (1.0-temp)*v[i][j][k];
	  w[i][j][k]=temp*fw[i][j][k] + (1.0-temp)*w[i][j][k];
	  T[i][j][k]=temp*fT[i][j][k] + (1.0-temp)*T[i][j][k];
	}
      }
    }

    printf("%d \n", n);

    /** Make sure T is nonnegative everywhere 	**/
    /** And check for conservation of mass, so	**/
    /** that avg layer deth is mass		**/
    
//change04      
    cmass=0.0;
    tsmooth=0;
    for(i=1;i<=Nx;i++){
      for(j=1;j<=Ny;j++){
	for(k=1;k<=Nz;k++){
	  if(T[i][j][k]<0.0){
	    printf("negative T");
	    tsmooth+=1;
	    T[i][j][k]=pow(10.0,-10.0);
	  }
	  cmass = cmass + ph[i][j][k];
	}
      }
    }
    cmass=cmass*h/(1.0*Nx*Ny*Nz);
    
    temp=0.0;
    for(i=1;i<=Nx;i++){
      for(j=1;j<=Ny;j++){
	for(k=1;k<=Nz;k++){
	  fph[i][j][k] = ph[i][j][k] *
	    exp(-100.0*(ph[i][j][k]-numax*6.0/pi/2.0) *
		(ph[i][j][k]-numax*6.0/pi/2.0));
	  temp+=fph[i][j][k];
	}
      }
    }
    for(i=1;i<=Nx;i++){
      for(j=1;j<=Ny;j++){
	for(k=1;k<=Nz;k++){
	  ph[i][j][k]=ph[i][j][k]
	    +fph[i][j][k]*(mass-cmass)/(h*temp/(Nx*Ny*Nz));
	}
      }
    }
//end change04

    printf("cmass %lf\n", cmass);

    /** Now set boundary conditions		**/
    for(i=1;i<=Nx;i++){
      for(j=1;j<=Ny;j++){
	w[i][j][1]=0.0;			/* impenetrable walls 	*/
	w[i][j][Nz]=0.0;		/* impenetrable walls	*/


      /* fullslip b.c. for horizontal velocities */

	/*****  FIX THIS ***********************/
	
	u[i][j][1]=(4.0*u[i][j][2]-u[i][j][3])/3.0;	
	u[i][j][Nz]=(4.0*u[i][j][Nz-1]-u[i][j][Nz-2])/3.0;

	
	v[i][j][1]=(4.0*v[i][j][2]-v[i][j][3])/3.0;	
	v[i][j][Nz]=(4.0*v[i][j][Nz-1]-v[i][j][Nz-2])/3.0;

      /* no temp flux through the wall 		*/
      
      T[i][j][1]=(4.0*T[i][j][2]-T[i][j][3])/3.0;
      T[i][j][Nz]=(4.0*T[i][j][Nz-1]-T[i][j][Nz-2])/3.0;

      if(T[i][j][1]<0.0){
	T[i][j][1]=0.0;}
      if(T[i][j][Nz]<0.0){
	T[i][j][Nz]=0.0;}
      }
    }
    

    /*******  Now calculate variables for iteration of equations ****/

    /******* Calculate first derivatives of the flow variables ***/
     
    /*First in x */
      
    for(k=1;k<=Nz;k++){
      for(j=1;j<=Ny;j++){
	for(i=2;i<Nx;i++){
	  phx[i][j][k]=(ph[i+1][j][k]-ph[i-1][j][k])/(2.0 * dx);
	  ux[i][j][k]=(u[i+1][j][k]-u[i-1][j][k])/(2.0 * dx);
	  vx[i][j][k]=(v[i+1][j][k]-v[i-1][j][k])/(2.0 * dx);
	  wx[i][j][k]=(w[i+1][j][k]-w[i-1][j][k])/(2.0 * dx);
	  Tx[i][j][k]=(T[i+1][j][k]-T[i-1][j][k])/(2.0 * dx);
	  Px[i][j][k]=(P[i+1][j][k]-P[i-1][j][k])/(2.0 * dx);
	}
	phx[1][j][k]=(ph[2][j][k]-ph[Nx][j][k])/(2.0*dx);
	ux[1][j][k]=(u[2][j][k]-u[Nx][j][k])/(2.0*dx);
	vx[1][j][k]=(v[2][j][k]-v[Nx][j][k])/(2.0*dx);
	wx[1][j][k]=(w[2][j][k]-w[Nx][j][k])/(2.0*dx);
	Tx[1][j][k]=(T[2][j][k]-T[Nx][j][k])/(2.0*dx);
	Px[1][j][k]=(P[2][j][k]-P[Nx][j][k])/(2.0*dx);
      
	phx[Nx][j][k]=(ph[1][j][k]-ph[Nx-1][j][k])/(2.0*dx);
	ux[Nx][j][k]=(u[1][j][k]-u[Nx-1][j][k])/(2.0*dx);
	vx[Nx][j][k]=(v[1][j][k]-v[Nx-1][j][k])/(2.0*dx);
	wx[Nx][j][k]=(w[1][j][k]-w[Nx-1][j][k])/(2.0*dx);
	Tx[Nx][j][k]=(T[1][j][k]-T[Nx-1][j][k])/(2.0*dx);
	Px[Nx][j][k]=(P[1][j][k]-P[Nx-1][j][k])/(2.0*dx);    
      }

    }
    /*Then in y */

     
    for(i=1;i<=Nx;i++){
      for(k=1;k<=Nz;k++){
	for(j=2;j<Ny;j++){
	  phy[i][j][k]=(ph[i][j+1][k]-ph[i][j-1][k])/(2.0 * dy);
	  uy[i][j][k]=(u[i][j+1][k]-u[i][j-1][k])/(2.0 * dy);
	  vy[i][j][k]=(v[i][j+1][k]-v[i][j-1][k])/(2.0 * dy);
	  wy[i][j][k]=(w[i][j+1][k]-w[i][j-1][k])/(2.0 * dy);
	  Ty[i][j][k]=(T[i][j+1][k]-T[i][j-1][k])/(2.0 * dy);
	  Py[i][j][k]=(P[i][j+1][k]-P[i][j-1][k])/(2.0 * dy);
	}
	phy[i][1][k]=(ph[i][2][k]-ph[i][Ny][k])/(2.0*dy);
	uy[i][1][k]=(u[i][2][k]-u[i][Ny][k])/(2.0*dy);
	vy[i][1][k]=(v[i][2][k]-v[i][Ny][k])/(2.0*dy);
	wy[i][1][k]=(w[i][2][k]-w[i][Ny][k])/(2.0*dy);
	Ty[i][1][k]=(T[i][2][k]-T[i][Ny][k])/(2.0*dy);
	Py[i][1][k]=(P[i][2][k]-P[i][Ny][k])/(2.0*dy);
      
	phy[i][Ny][k]=(ph[i][1][k]-ph[i][Ny-1][k])/(2.0*dy);
	uy[i][Ny][k]=(u[i][1][k]-u[i][Ny-1][k])/(2.0*dy);
	vy[i][Ny][k]=(v[i][1][k]-v[i][Ny-1][k])/(2.0*dy);
	wy[i][Ny][k]=(w[i][1][k]-w[i][Ny-1][k])/(2.0*dy);
	Ty[i][Ny][k]=(T[i][1][k]-T[i][Ny-1][k])/(2.0*dy);
	Py[i][Ny][k]=(P[i][1][k]-P[i][Ny-1][k])/(2.0*dy);
      }
    }
	

     /* finally in z */
    for(i=1;i<=Nx;i++){
      for(j=1;j<=Ny;j++){      
	for(k=2;k<Nz;k++){
	  phz[i][j][k]=(ph[i][j][k+1]-ph[i][j][k-1])/(2.0*dz);
	  uz[i][j][k]=(u[i][j][k+1]-u[i][j][k-1])/(2.0*dz);
	  vz[i][j][k]=(v[i][j][k+1]-v[i][j][k-1])/(2.0*dz);
	  wz[i][j][k]=(w[i][j][k+1]-w[i][j][k-1])/(2.0*dz);
	  Tz[i][j][k]=(T[i][j][k+1]-T[i][j][k-1])/(2.0*dz);
	  Pz[i][j][k]=(P[i][j][k+1]-P[i][j][k-1])/(2.0*dz);

      }


	
      phz[i][j][1]=(-ph[i][j][3] + 4.0*ph[i][j][2]-3.0*ph[i][j][1])/(2.0*dz);
      uz[i][j][1]=(-u[i][j][3] + 4.0*u[i][j][2]-3.0*u[i][j][1])/(2.0*dz);
      vz[i][j][1]=(-v[i][j][3] + 4.0*v[i][j][2]-3.0*v[i][j][1])/(2.0*dz);
      wz[i][j][1]=(-w[i][j][3] + 4.0*w[i][j][2]-3.0*w[i][j][1])/(2.0*dz);
      Tz[i][j][1]=(-T[i][j][3] + 4.0*T[i][j][2]-3.0*T[i][j][1])/(2.0*dz);
      Pz[i][j][1]=(-P[i][j][3] + 4.0*P[i][j][2]-3.0*P[i][j][1])/(2.0*dz);

      phz[i][j][Nz]=(ph[i][j][Nz-2] - 4.0*ph[i][j][Nz-1]+
		     3.0*ph[i][j][Nz])/(2.0*dz);
      uz[i][j][Nz]=(u[i][j][Nz-2] - 4.0*u[i][j][Nz-1]+
		    3.0*u[i][j][Nz])/(2.0*dz);
      vz[i][j][Nz]=(v[i][j][Nz-2] - 4.0*v[i][j][Nz-1]+
		    3.0*v[i][j][Nz])/(2.0*dz);
      wz[i][j][Nz]=(w[i][j][Nz-2] - 4.0*w[i][j][Nz-1]+
		    3.0*w[i][j][Nz])/(2.0*dz);
      Tz[i][j][Nz]=(T[i][j][Nz-2] - 4.0*T[i][j][Nz-1]+
		    3.0*T[i][j][Nz])/(2.0*dz);
      Pz[i][j][Nz]=(P[i][j][Nz-2] - 4.0*P[i][j][Nz-1]+
		    3.0*P[i][j][Nz])/(2.0*dz);
      
      }
    }


    
    /**** Now calculate combinations of the variables and derivatives****/
//change05
    // for shared(G,P,lam,mu,K,T,gam,divU,stressxx,ux,wz,vy,ph,stressyy,stresszz,stressxy,stressxz,stressyz,heatfluxx,heatfluxy,heatfluxz,vx,uy,wx,uz,wy,vz,Tx,Ty,Tz) private(i,j,k) num_threads(4)
    for(i=1;i<=Nx;i++){
      for(j=1;j<=Ny;j++){
	for(k=1;k<=Nz;k++){
    /**** First calculate transport coefficients *******/


	  nu = ph[i][j][k]*pi/6.0;
	  G[i][j][k]=nu/(1.0-pow((nu/numax),(4.0*numax/3.0)));
	  P[i][j][k]=(1.0+2.0*(onepluse)*G[i][j][k])*T[i][j][k]*ph[i][j][k];
	  if(P[i][j][k]<0.0){
		  P[i][j][k]=0.0;}
	  tone=ph[i][j][k]*sqrt(T[i][j][k]);
	  lam[i][j][k] = 8.0/3.0/sqrt(pi)*tone*G[i][j][k];
	  mu[i][j][k] = sqrt(pi) / 6.0 * tone *
	    (1.0 +4.0/5.0*(1.0+12.0/pi)*G[i][j][k]);
	  K[i][j][k] = 15.0 * sqrt(pi)/16.0 * tone *
	    (1.0 + 6.0/5.0 * (1.0 + 32.0/9.0/pi)*G[i][j][k]);
	  ev = 1.0 - (1.0 -e) * pow((sqrt(T[i][j][k])/1.0),.75);
	  if(T[i][j][k]>1){
	    ev = e;}
	  ed = 12.0*(1.0-ev*ev)/sqrt(pi);
	  gam[i][j][k] = ed*G[i][j][k]*T[i][j][k]*tone; 


	  /**now the components of the stress tensor**/        

	  divU[i][j][k]=ux[i][j][k]+wz[i][j][k]+vy[i][j][k];
	  stressxx[i][j][k] = 2.0*mu[i][j][k]*ux[i][j][k] + 
	    (lam[i][j][k] - 2.0/3.0*mu[i][j][k])*divU[i][j][k];
	    /*sqrt(2.0*T[i][j][k]*(4.0/3.0*mu[i][j][k]+lam[i][j][k])) * 
	      gasdevrand();*/
	  stressyy[i][j][k] = 2.0*mu[i][j][k]*vy[i][j][k] + 
	    (lam[i][j][k] - 2.0/3.0*mu[i][j][k])*divU[i][j][k]; 
	  /*sqrt(2.0*T[i][j][k]*(4.0/3.0*mu[i][j][k]+lam[i][j][k])) * 
	    gasdevrand();*/		 
	  stresszz[i][j][k] = 2.0*mu[i][j][k]*wz[i][j][k] + 
	    (lam[i][j][k] - 2.0/3.0*mu[i][j][k])*divU[i][j][k]; 
	  /*sqrt(2.0*T[i][j][k]*(4.0/3.0*mu[i][j][k]+lam[i][j][k])) * 
	    gasdevrand();*/
	  stressxy[i][j][k] = mu[i][j][k] * (vx[i][j][k]+uy[i][j][k]); 
	  /* +  sqrt(2.0*T[i][j][k]*mu[i][j][k])*gasdevrand();*/
	  stressxz[i][j][k] = mu[i][j][k] * (wx[i][j][k]+uz[i][j][k]);  
	  /*+   sqrt(2.0*T[i][j][k]*mu[i][j][k])*gasdevrand();*/
	  stressyz[i][j][k] = mu[i][j][k] * (wy[i][j][k]+vz[i][j][k]); 
	  /*+  sqrt(2.0*T[i][j][k]*mu[i][j][k])*gasdevrand();*/
	  heatfluxx[i][j][k]=-K[i][j][k]*Tx[i][j][k]; 
	  /* +sqrt(2.0*K[i][j][k])*T[i][j][k]*gasdevrand(); */
	  heatfluxy[i][j][k]=-K[i][j][k]*Ty[i][j][k];
	    /* +   sqrt(2.0*K[i][j][k])*T[i][j][k]*gasdevrand();*/
	  heatfluxz[i][j][k]=-K[i][j][k]*Tz[i][j][k]; 
	  /*+   sqrt(2.0*K[i][j][k])*T[i][j][k]*gasdevrand();*/
	}
      }
    }
//end change05
  
 
    /* and higher order derivatives and combinations */

    /* first x derivatives */

    for(j=1;j<=Ny;j++){
      for(k=1;k<=Nz;k++){
	for(i=2;i<Nx;i++){
	  divflux[i][j][k]=(heatfluxx[i+1][j][k]-
			      heatfluxx[i-1][j][k])/(2.0*dx);
	  stressxxDx[i][j][k]=(stressxx[i+1][j][k]-
			       stressxx[i-1][j][k])/(2.0*dx); 
	  stressxzDx[i][j][k]=(stressxz[i+1][j][k]-
			   stressxz[i-1][j][k])/(2.0*dx);
	  stressxyDx[i][j][k]=(stressxy[i+1][j][k]-
			   stressxy[i-1][j][k])/(2.0*dx);
	}

	divflux[1][j][k]=(heatfluxx[2][j][k] -
			    heatfluxx[Nx][j][k])/(2.0*dx);
	stressxxDx[1][j][k]=(stressxx[2][j][k] -
			     stressxx[Nx][j][k])/(2.0*dx);
	stressxzDx[1][j][k]=(stressxz[2][j][k] -
			 stressxz[Nx][j][k])/(2.0*dx);
	stressxyDx[1][j][k]=(stressxy[2][j][k] -
			 stressxy[Nx][j][k])/(2.0*dx);
	
	divflux[Nx][j][k]=(heatfluxx[1][j][k] -
			    heatfluxx[Nx-1][j][k])/(2.0*dx);
	stressxxDx[Nx][j][k]=(stressxx[1][j][k] -
			      stressxx[Nx-1][j][k])/(2.0*dx);
	stressxzDx[Nx][j][k]=(stressxz[1][j][k] -
			  stressxz[Nx-1][j][k])/(2.0*dx);
	stressxyDx[Nx][j][k]=(stressxy[1][j][k] -
			  stressxy[Nx-1][j][k])/(2.0*dx);
	
      }
    }
    /* now y derivatives */
    
    for(i=1;i<=Nx;i++){
      for(k=1;k<=Nz;k++){
	for(j=2;j<Ny;j++){
	  divflux[i][j][k] += (heatfluxy[i][j+1][k]
				 -heatfluxy[i][j-1][k])/(2.0 * dy);
	  stressyyDy[i][j][k]=(stressyy[i][j+1][k] - 
			       stressyy[i][j-1][k])/(2.0*dy);
	  stressxyDy[i][j][k]=(stressxy[i][j+1][k] - 
			       stressxy[i][j-1][k])/(2.0 * dy);
	  stressyzDy[i][j][k]=(stressyz[i][j+1][k] - 
			       stressyz[i][j-1][k])/(2.0 * dy);


      }

      divflux[i][1][k]+=(heatfluxy[i][2][k]-
			   heatfluxy[i][Ny][k])/(2.0*dy);
      stressyyDy[i][1][k]=(stressyy[i][2][k]-stressyy[i][Ny][k])/(2.0*dy);
      stressxyDy[i][1][k]=(stressxy[i][2][k]-stressxy[i][Ny][k])/(2.0*dy);
      stressyzDy[i][1][k]=(stressyz[i][2][k]-stressyz[i][Ny][k])/(2.0*dy);

      
      divflux[i][Ny][k]+=(heatfluxy[i][1][k]-
			   heatfluxy[i][Ny-1][k])/(2.0*dy);
      stressyyDy[i][Ny][k]=(stressyy[i][1][k]-stressyy[i][Ny-1][k])/(2.0*dy);
      stressxyDy[i][Ny][k]=(stressxy[i][1][k]-stressxy[i][Ny-1][k])/(2.0*dy);
      stressyzDy[i][Ny][k]=(stressyz[i][1][k]-stressyz[i][Ny-1][k])/(2.0*dy);
      }
    }

      /* then z derivatives */

    for(i=1;i<=Nx;i++){    
      for(j=1;j<=Ny;j++){
	for(k=2;k<Nz;k++){
	  divflux[i][j][k] +=
	    (heatfluxz[i][j][k+1]-heatfluxz[i][j][k-1])/
	    (2.0*dz);
	  stresszzDz[i][j][k]=(stresszz[i][j][k+1] - 
			       stresszz[i][j][k-1])/(2.0*dz);
	  stressxzDz[i][j][k]=(stressxz[i][j][k+1] - 
			       stressxz[i][j][k-1])/(2.0*dz);
	  stressyzDz[i][j][k]=(stressyz[i][j][k+1] - 
			       stressyz[i][j][k-1])/(2.0*dz);
      }

	divflux[i][j][1] +=
	  (-heatfluxz[i][j][3] +
	   4.0*heatfluxz[i][j][2] -
	   3.0*heatfluxz[i][j][1])/
	  (2.0*dz);
	stresszzDz[i][j][1]=
	  (-stresszz[i][j][3] +
	   4.0*stresszz[i][j][2]-
	   3.0*stresszz[i][j][1])/
	  (2.0*dz);
	stressxzDz[i][j][1]=
	  (-stressxz[i][j][3] +
	   4.0*stressxz[i][j][2]-
	   3.0*stressxz[i][j][1])/
	  (2.0*dz);
	stressyzDz[i][j][1]=
	  (-stressyz[i][j][3] +
	   4.0*stressyz[i][j][2]-
	   3.0*stressyz[i][j][1])/
	  (2.0*dz);
      
	divflux[i][j][Nz] +=
	  (heatfluxz[i][j][Nz-2] -
	   4.0*heatfluxz[i][j][Nz-1]
	   +3.0*heatfluxz[i][j][Nz])/
	  (2.0*dz);
	stresszzDz[i][j][Nz]=
	  (stresszz[i][j][Nz-2] -
	   4.0*stresszz[i][j][Nz-1] +
	   3.0*stresszz[i][j][Nz])/(2.0*dz);
	stressxzDz[i][j][Nz]=
	  (stressxz[i][j][Nz-2] -
	   4.0*stressxz[i][j][Nz-1] +
	   3.0*stressxz[i][j][Nz])/(2.0*dz);
	stressyzDz[i][j][Nz]=
	  (stressyz[i][j][Nz-2] -
	   4.0*stressyz[i][j][Nz-1] +
	   3.0*stressyz[i][j][Nz])/(2.0*dz);

      }
    }
    

    /********* Finally, calculate dph, dv, dw, dT *****************/
    for(i=1;i<=Nx;i++){
      for(j=1;j<=Ny;j++){
	for(k=1;k<=Nz;k++){
	  dph[i][j][k] = -(u[i][j][k] * phx[i][j][k] +
			   v[i][j][k] * phy[i][j][k] +
			   w[i][j][k] * phz[i][j][k] +
			   ph[i][j][k] * divU[i][j][k]);

	  dw[i][j][k] =
	    -u[i][j][k]*wx[i][j][k]
	    -v[i][j][k]*wy[i][j][k]
	    -w[i][j][k]*wz[i][j][k]
	    + (-Pz[i][j][k] +
	       stresszzDz[i][j][k] +
	       stressxzDx[i][j][k] +
	       stressyzDy[i][j][k]) /ph[i][j][k]
	    - gr - GAM*cos(2.0*pi*f*t);

	  du[i][j][k]=
	    -u[i][j][k]*ux[i][j][k]
	    -v[i][j][k]*uy[i][j][k]
	    -w[i][j][k]*uz[i][j][k]
	    + (-Px[i][j][k] +
	       stressxxDx[i][j][k] +
	       stressxyDy[i][j][k] +
	       stressxzDz[i][j][k])/ph[i][j][k]; 
	       
	    
	  dv[i][j][k] =
	    -u[i][j][k]*vx[i][j][k]
	    -v[i][j][k]*vy[i][j][k]
	    -w[i][j][k]*vz[i][j][k]
	    + (-Py[i][j][k]+
	       stressyyDy[i][j][k]+
	       stressxyDx[i][j][k]+
	       stressyzDz[i][j][k])/ph[i][j][k];

	  dT[i][j][k] =
	    -u[i][j][k]*Tx[i][j][k]
	    -v[i][j][k]*Ty[i][j][k]
	    -w[i][j][k]*Tz[i][j][k]
	    + (-P[i][j][k] * divU[i][j][k] +
	       stressxx[i][j][k]*ux[i][j][k] +
	       stressyy[i][j][k]*vy[i][j][k] +
	       stresszz[i][j][k]*wz[i][j][k] +
	       stressxy[i][j][k]*(uy[i][j][k]+vx[i][j][k])+
	       stressyz[i][j][k]*(vz[i][j][k]+wy[i][j][k])+
	       stressxz[i][j][k]*(uz[i][j][k]+wx[i][j][k])-
	       divflux[i][j][k] -gam[i][j][k])/
	    (1.5*ph[i][j][k]);
	}
      }
    }

    
      /* Write data to file */
   
    /*    if(n%1000==2){	*/
    
    printf("%lf\n",t*f*numsnaps-tsnap); 
    
    /*if time elapsed since last save is >=1/numsnaps, then save*/

    if((time*numsnaps-tsnap)>=1.0){

      if(cmass>mass/2.0 && cmass<mass*2.0){
      
      
      tsnap=floor(t*f*numsnaps); /*(set counter at last multiple of 1/numsnaps*/

      /*now save out the variables */

      printf("time=%lf\n",time);
    
      savecount +=1;
      printf("savecount= %d\n",savecount);
      
    sprintf(file_name,"%s%04d.dat",base,savecount); 

    output1 = fopen(file_name,"wb");

    /** zero size buffer.  Slows down write, but doesn't use RAM for buffer**/
    setvbuf(output1,NULL,_IONBF,0);


    fwrite(&fstar,sizeof(double),1,output1);
    fwrite(&GAM,sizeof(double),1,output1);
    fwrite(&Nx,sizeof(int),1,output1);
    fwrite(&Ny,sizeof(int),1,output1);
    fwrite(&Nz,sizeof(int),1,output1);
    fwrite(&n,sizeof(int),1,output1);
    fwrite(&time,sizeof(double),1,output1);
    fwrite(&numsmooth,sizeof(int),1,output1);
    fwrite(&sizex,sizeof(double),1,output1);
    fwrite(&sizey,sizeof(double),1,output1);
    fwrite(&h,sizeof(double),1,output1);
    fwrite(&depth,sizeof(double),1,output1);

    /*    
    fprintf(output,"%lf %lf \n",fstar, GAM);
    
    fprintf(output,"%d %d %d \n",Nx, Ny, Nz);

    fprintf(output,"%d %lf %d\n",n,t*f,numsmooth);

    fprintf(output," %lf %lf %lf %lf\n", ggx, ggy,h,mass);
    */
    /*
    fwrite(ph,sizeof(double),Nx*Ny*Nz,output1);
    fwrite(u,sizeof(double),Nx*Ny*Nz,output1);
    fwrite(v,sizeof(double),Nx*Ny*Nz,output1);
    fwrite(w,sizeof(double),Nx*Ny*Nz,output1);
    fwrite(T,sizeof(double),Nx*Ny*Nz,output1);
    */

    
    for(i=1;i<=Nx;i++){
      for(j=1;j<=Ny;j++){
	for(k=1;k<=Nz;k++){
	 fwrite(&ph[i][j][k],sizeof(double),1,output1); 
	}
      }
    }
    
    for(i=1;i<=Nx;i++){
      for(j=1;j<=Ny;j++){
	for(k=1;k<=Nz;k++){
	  fwrite(&u[i][j][k],sizeof(double),1,output1); 
	}
      }
    }
    
    for(i=1;i<=Nx;i++){
      for(j=1;j<=Ny;j++){
	for(k=1;k<=Nz;k++){
	  fwrite(&v[i][j][k],sizeof(double),1,output1); 
	}
      }
    }
    
    for(i=1;i<=Nx;i++){
      for(j=1;j<=Ny;j++){
	for(k=1;k<=Nz;k++){
	  fwrite(&w[i][j][k],sizeof(double),1,output1); 
	}
      }
    }

  
    for(i=1;i<=Nx;i++){    
      for(j=1;j<=Ny;j++){
	for(k=1;k<=Nz;k++){
	  fwrite(&T[i][j][k],sizeof(double),1,output1);
	}
      }
    }

    
    /*close data file*/

    fclose(output1);
    output1=0;
      }
     else{
	printf("mass is outside bounds\n");
        checkforsignaltoend=0;
      }
    }
    /* Use additional smoothing to keep stable */
    /*smooth if there was a negative temperature previous timestep*/
    
    if(tsmooth>0){	

      printf("smoothing \n");

      numsmooth+=1; /*counts the number of times smoothing happens during run*/
      
      for(i=1;i<=Nx;i++){
	for(j=1;j<=Ny;j++){
	  for(k=2;k<Nz;k++){
	    fph[i][j][k]=(.5*ph[i][j][k-1] + ph[i][j][k] +
			  .5*ph[i][j][k+1])/2.0;
	    fu[i][j][k]=(.5*u[i][j][k-1] + u[i][j][k] +
			 .5*u[i][j][k+1])/2.0;
	    fv[i][j][k]=(.5*v[i][j][k-1] + v[i][j][k] +
			 .5*v[i][j][k+1])/2.0;
	    fw[i][j][k]=(.5*w[i][j][k-1] + w[i][j][k] +
			 .5*w[i][j][k+1])/2.0;
	    fT[i][j][k]=(.5*T[i][j][k-1] + T[i][j][k] +
			 .5*T[i][j][k+1])/2.0;
	  }
	
	  fph[i][j][1]=(.75*ph[i][j][1] + .25*ph[i][j][2]);
	  fu[i][j][1]=(.75*u[i][j][1] + .25*u[i][j][2]);
	  fv[i][j][1]=(.75*v[i][j][1] + .25*v[i][j][2]);
	  fw[i][j][1]=(.75*w[i][j][1] + .25*w[i][j][2]);
	  fT[i][j][1]=(.75*T[i][j][1] + .25*T[i][j][2]);

	  fph[i][j][Nz]=(.75*ph[i][j][Nz] + .25*ph[i][j][Nz-1]);
	  fu[i][j][Nz]=(.75*u[i][j][Nz] + .25*u[i][j][Nz-1]);
	  fv[i][j][Nz]=(.75*v[i][j][Nz] + .25*v[i][j][Nz-1]);
	  fw[i][j][Nz]=(.75*w[i][j][Nz] + .25*w[i][j][Nz-1]);
	  
	  fT[i][j][Nz]=(.75*T[i][j][Nz] + .25*T[i][j][Nz-1]);
	}
      }
    
    for(i=1;i<=Nx;i++){
      for(k=1;k<=Nz;k++){
	for(j=2;j<Ny;j++){
	  fph[i][j][k]=(.01*fph[i][j-1][k] + fph[i][j][k] +
		       .01*fph[i][j+1][k])/1.02;
	  fu[i][j][k]=(.01*fu[i][j-1][k] + fu[i][j][k] +
		      .01*fu[i][j+1][k])/1.02;
	  fv[i][j][k]=(.01*fv[i][j-1][k] + fv[i][j][k] +
		      .01*fv[i][j+1][k])/1.02;
	  fw[i][j][k]=(.01*fw[i][j-1][k] + fw[i][j][k] +
		      .01*fw[i][j+1][k])/1.02;
	  fT[i][j][k]=(.01*fT[i][j-1][k] + fT[i][j][k] +
		      .01*fT[i][j+1][k])/1.02;
	}
	fph[i][1][k]=(.01*fph[i][Ny][k] + fph[i][1][k] +
		     .01*fph[i][2][k])/1.02;
	fu[i][1][k]=(.01*fu[i][Ny][k] + fu[i][1][k] +
		    .01*fu[i][2][k])/1.02;
	fv[i][1][k]=(.01*fv[i][Ny][k] + fv[i][1][k] +
		    .01*fv[i][2][k])/1.02;
	fw[i][1][k]=(.01*fw[i][Ny][k] + fw[i][1][k] +
		    .01*fw[i][2][k])/1.02;
	fT[i][1][k]=(.01*fT[i][Ny][k] + fT[i][1][k] +
		    .01*fT[i][2][k])/1.02;
	
	fph[i][Ny][k]=(.01*fph[i][Ny-1][k] + fph[i][Ny][k] +
		      .01*fph[i][1][k])/1.02;
	fu[i][Ny][k]=(.01*fu[i][Ny-1][k] + fu[i][Ny][k] +
		       .01*fu[i][1][k])/1.02;
	fv[i][Ny][k]=(.01*fv[i][Ny-1][k] + fv[i][Ny][k] +
		     .01*fv[i][1][k])/1.02;
	fw[i][Ny][k]=(.01*fw[i][Ny-1][k] + fw[i][Ny][k] +
		     .01*fw[i][1][k])/1.02;
	fT[i][Ny][k]=(.01*fT[i][Ny-1][k] + fT[i][Ny][k] +
		     .01*fT[i][1][k])/1.02;
	
      }
    }

    for(j=1;j<=Ny;j++){
      for(k=1;k<=Nz;k++){
	for(i=2;i<Nx;i++){
	  
	  ph[i][j][k]=(.01*fph[i-1][j][k] + fph[i][j][k] +
		       .01*fph[i+1][j][k])/1.02;
	  u[i][j][k]=(.01*fu[i-1][j][k] + fu[i][j][k] +
			.01*fu[i+1][j][k])/1.02;
	  v[i][j][k]=(.01*fv[i-1][j][k] + fv[i][j][k] +
			.01*fv[i+1][j][k])/1.02;
	  w[i][j][k]=(.01*fw[i-1][j][k] + fw[i][j][k] +
		      .01*fw[i+1][j][k])/1.02;
	  T[i][j][k]=(.01*fT[i-1][j][k] + fT[i][j][k] +
		      .01*fT[i+1][j][k])/1.02;

	}

	  ph[1][j][k]=(.01*fph[Nx][j][k] + fph[1][j][k] +
		       .01*fph[2][j][k])/1.02;
	  u[1][j][k]=(.01*fu[Nx][j][k] + fu[1][j][k] +
		      .01*fu[2][j][k])/1.02;
	  v[1][j][k]=(.01*fv[Nx][j][k] + fv[1][j][k] +
		      .01*fv[2][j][k])/1.02;
	  w[1][j][k]=(.01*fw[Nx][j][k] + fw[1][j][k] +
		      .01*fw[2][j][k])/1.02;
	  T[1][j][k]=(.01*fT[Nx][j][k] + fT[1][j][k] +
		      .01*fT[2][j][k])/1.02;
	
	  ph[Nx][j][k]=(.01*fph[Nx-1][j][k] + fph[Nx][j][k]+
			.01*fph[1][j][k])/1.02;
	  u[Nx][j][k]=(.01*fu[Nx-1][j][k] + fu[Nx][j][k] +
		       .01*fu[1][j][k])/1.02;
	  v[Nx][j][k]=(.01*fv[Nx-1][j][k] + fv[Nx][j][k] +
		       .01*fv[1][j][k])/1.02;
	  w[Nx][j][k]=(.01*fw[Nx-1][j][k] + fw[Nx][j][k] +
		       .01*fw[1][j][k])/1.02;
	  T[Nx][j][k]=(.01*fT[Nx-1][j][k] + fT[Nx][j][k] +
		       .01*fT[1][j][k])/1.02;
	
	}

      }
    }
    
    /***** find length of timestep *******/
	    
    temp=dz;
    if(dy<temp){
      temp=dy;
    }
    if(dx<temp){
      temp=dx;
    }
    
    tmp=0.0;
    for(i=1;i<=Nx;i++){
      for(j=1;j<=Ny;j++){
	for(k=1;k<=Nz;k++){
	  if(K[i][j][k]/ph[i][j][k]>tmp){
	    tmp=K[i][j][k]/ph[i][j][k];}
	  if(lam[i][j][k]/ph[i][j][k]>tmp){
	    tmp=lam[i][j][k]/ph[i][j][k];}
	  if(mu[i][j][k]/ph[i][j][k]>tmp){
	    tmp=mu[i][j][k]/ph[i][j][k];}
	}
      }
    }
    dtone=s * temp * temp / tmp;

    tmp=0.0;
    for(i=1;i<=Nx;i++){
      for(j=1;j<=Ny;j++){
	for(k=1;k<=Nz;k++){
	  if(fabs(w[i][j][k])>tmp){
	    tmp=fabs(w[i][j][k]);}
	  if(fabs(v[i][j][k])>tmp){
	    tmp=fabs(v[i][j][k]);}
	   if(fabs(u[i][j][k])>tmp){
	    tmp=fabs(u[i][j][k]);}
	}
      }
    }

     if(C*temp/tmp < dtone){
       dtone=C*temp/tmp;}

     tmp=0.0;
     for(i=1;i<=Nx;i++){
       for(j=1;j<=Ny;j++){
	 for(k=1;k<=Nz;k++){
	   if(G[i][j][k]*sqrt(T[i][j][k])>tmp){
	     tmp=(G[i][j][k]*sqrt(T[i][j][k]));}
	 }
       }
     }
       if(tau/gamo/tmp<dtone){
	 dtone=tau/gamo/tmp;}

       if(mindt>dtone){
	 dtone=mindt;}

       if(dt>dtone){
	 dt=dtone;}

       if(dt<dtone){
	 dt += dtstep/2.0;}

       if(dt>maxdt){
	 dt=maxdt;}
       
       for(i=1;i<=Nx;i++){
	 for(j=1;j<=Ny;j++){
	   for(k=1;k<=Nz;k++){
	     ph[i][j][k]=ph[i][j][k]+dph[i][j][k]*dt;
	     u[i][j][k]=u[i][j][k]+du[i][j][k]*dt;
	     v[i][j][k]=v[i][j][k]+dv[i][j][k]*dt;
	     w[i][j][k]=w[i][j][k]+dw[i][j][k]*dt;
	     T[i][j][k]=T[i][j][k]+dT[i][j][k]*dt;
	     if(ph[i][j][k]<pow(10.0,-5.0)){
	       ph[i][j][k]=pow(10.0,-5.0);}
	   }	
	 }
       }
       t=t+dt;
       time=t*f;
    /*********************END TIMESTEPPING ****************/

       
  }
  end=clock();
  float tottime=((float)(end-start))/CLOCKS_PER_SEC;
  printf("Time elapsed: %f",tottime);

  printf("freeing memory\n");
       
  /* FREE MEMORY ALLOCATED PREVIOUSLY */

  free_dvector(x,1,Nx);
  free_dvector(y,1,Ny);
  free_dvector(z,1,Nz);

  free_dthreearray(ph,1,Nx,1,Ny,1,Nz);
  free_dthreearray(T,1,Nx,1,Ny,1,Nz);
  free_dthreearray(u,1,Nx,1,Ny,1,Nz);
  free_dthreearray(v,1,Nx,1,Ny,1,Nz);
  free_dthreearray(w,1,Nx,1,Ny,1,Nz);
  free_dthreearray(P,1,Nx,1,Ny,1,Nz);
  free_dthreearray(G,1,Nx,1,Ny,1,Nz);

  free_dthreearray(fph,1,Nx,1,Ny,1,Nz);
  free_dthreearray(fT,1,Nx,1,Ny,1,Nz);
  free_dthreearray(fu,1,Nx,1,Ny,1,Nz);
  free_dthreearray(fv,1,Nx,1,Ny,1,Nz);
  free_dthreearray(fw,1,Nx,1,Ny,1,Nz);

  free_dthreearray(dph,1,Nx,1,Ny,1,Nz);
  free_dthreearray(dT,1,Nx,1,Ny,1,Nz);
  free_dthreearray(du,1,Nx,1,Ny,1,Nz);
  free_dthreearray(dv,1,Nx,1,Ny,1,Nz);
  free_dthreearray(dw,1,Nx,1,Ny,1,Nz);

  free_dthreearray(phx,1,Nx,1,Ny,1,Nz);
  free_dthreearray(Tx,1,Nx,1,Ny,1,Nz);
  free_dthreearray(ux,1,Nx,1,Ny,1,Nz);
  free_dthreearray(vx,1,Nx,1,Ny,1,Nz);
  free_dthreearray(wx,1,Nx,1,Ny,1,Nz);
  free_dthreearray(Px,1,Nx,1,Ny,1,Nz);

  free_dthreearray(phy,1,Nx,1,Ny,1,Nz);
  free_dthreearray(Ty,1,Nx,1,Ny,1,Nz);
  free_dthreearray(uy,1,Nx,1,Ny,1,Nz);
  free_dthreearray(vy,1,Nx,1,Ny,1,Nz);
  free_dthreearray(wy,1,Nx,1,Ny,1,Nz);
  free_dthreearray(Py,1,Nx,1,Ny,1,Nz);
  
  free_dthreearray(phz,1,Nx,1,Ny,1,Nz);
  free_dthreearray(Tz,1,Nx,1,Ny,1,Nz);
  free_dthreearray(uz,1,Nx,1,Ny,1,Nz);
  free_dthreearray(vz,1,Nx,1,Ny,1,Nz);
  free_dthreearray(wz,1,Nx,1,Ny,1,Nz);
  free_dthreearray(Pz,1,Nx,1,Ny,1,Nz);


  free_dthreearray(divU,1,Nx,1,Ny,1,Nz);
  free_dthreearray(mu,1,Nx,1,Ny,1,Nz);
  free_dthreearray(K,1,Nx,1,Ny,1,Nz);
  free_dthreearray(lam,1,Nx,1,Ny,1,Nz);
  free_dthreearray(gam,1,Nx,1,Ny,1,Nz);


  free_dthreearray(stressxx,1,Nx,1,Ny,1,Nz);
  free_dthreearray(stressyy,1,Nx,1,Ny,1,Nz);
  free_dthreearray(stresszz,1,Nx,1,Ny,1,Nz);
  free_dthreearray(stressyz,1,Nx,1,Ny,1,Nz);
  free_dthreearray(stressxy,1,Nx,1,Ny,1,Nz);
  free_dthreearray(stressxz,1,Nx,1,Ny,1,Nz);
  free_dthreearray(stressxxDx,1,Nx,1,Ny,1,Nz);
  free_dthreearray(stressyyDy,1,Nx,1,Ny,1,Nz);
  free_dthreearray(stresszzDz,1,Nx,1,Ny,1,Nz);


  free_dthreearray(stressxyDy,1,Nx,1,Ny,1,Nz);
  free_dthreearray(stressxzDz,1,Nx,1,Ny,1,Nz);
  free_dthreearray(stressyzDy,1,Nx,1,Ny,1,Nz);
  free_dthreearray(stressxzDx,1,Nx,1,Ny,1,Nz);
  free_dthreearray(stressyzDz,1,Nx,1,Ny,1,Nz);
  free_dthreearray(stressxyDx,1,Nx,1,Ny,1,Nz);

  free_dthreearray(divflux,1,Nx,1,Ny,1,Nz);
  free_dthreearray(heatfluxx,1,Nx,1,Ny,1,Nz);
  free_dthreearray(heatfluxy,1,Nx,1,Ny,1,Nz);
  free_dthreearray(heatfluxz,1,Nx,1,Ny,1,Nz);

  
  printf("done freeing memory\n");


  /*  while(1){  break;}*/

  /* END; returning zero because successful */
  return(0);

}
 










