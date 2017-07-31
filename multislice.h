
/* Constants */
#define planck 6.626176e-34 /*Plancks constant */
#define mass 9.10938356e-31 /* Electron Mass */
#define charge 1.60217e-19 /*Electron Charge */
#define pi 3.14159265359  /* Constant pi */ 
#define lightspeed 2.99792e8 /* Speed of Light */
#define lightspeed2 8.98755e16 /* Speed of Light Squared */ 
#define hartree 4.35974e-18 /* How mant hartrees in a Joule */
#define angstrom 1e-10 /* Angstrom in metres */
#define abberation 1.8e-3 /*Coefficient of Spherical Abberation (Microscope constant) */ 
#define df 94.25e-9 /* Defocus */


/* Variable declarations */ 
int i,j,loop1, height,width,depth,size,wavefunctionsize,fy, tx, ty, fx; 
double a,b,c,v,total,total1,dx,dy,dz,*x,*y,*V,abbfunction,*k_x2,*k_y2,wavelength,sigma;
double complex *T,*P,*wavefunction,*PSF,var,*P2,tl, tr, bl, br;
char rec;

/* Write width and height to an output file */
int outputparameters(void){
  FILE* output2;
  output2 = fopen("size.txt", "w");
  fprintf(output2, "%d \n%d\n", width,height);
  fclose(output2);
}

/* Spatial Frequency Calculation (x direction) */ 
double calck_x2(int width, double dx,double*x){
  for (i=0; i<wavefunctionsize;i++){
    k_x2[i] = pow((x[i]) / (width*pow(dx,2)),2);
  } 
}   

/* Spatial frequency Calculation (y direction) */
double calck_y2(int height, double dy,double*y){
   for (i=0; i<wavefunctionsize;i++){
      k_y2[i] = pow((y[i]) / (height*pow(dy,2)),2);
   }
}


/* Relativistic De broglie wavelength */ 
double calcwavelength(double v){
  wavelength = (planck*lightspeed) / sqrt(charge*v*(2*mass*lightspeed2+charge*v));
}



/* Interaction parameter */ 
double calcsigma(double wavelength, double v){
  sigma = ((2*pi)/(wavelength*v))*((mass*lightspeed2+charge*v)/(2*mass*lightspeed2+charge*v)); 
}




