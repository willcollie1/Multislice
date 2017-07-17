
/* Constants */
#define planck 6.63e-34 /*Plancks constant */
#define mass 9.11e-31 /* Electron Mass */
#define charge 1.69e-19 /*Electron Charge */
#define pi 3.14159265359  /* Constant pi */ 
#define lightspeed 2.99792e8 /* Speed of Light */
#define lightspeed2 8.98755e16 /* Speed of Light Squared */ 
#define hartree 4.359e-18 /* Units of the potential file - this many joules */
#define angstrom 1e-10 /* Angstrom in metres */
#define abberation 0.00001 /* Coefficient of Spherical Abberation (Microscope constant) */ 
#define df 180e-9 /* Defocus */ 



/* Variables used for looping */ 
int i,j,z,r = 0; 
double a,b,c,v;
int height,width,depth, size;




int outputparameters(void){
  FILE* output2;
  output2 = fopen("size.txt", "w");
  fprintf(output2, "%d \n%d\n", width,height);
  fclose(output2);
}



/* Calculates the wavelength given the accelerating voltage */ 
double calcwavelength(double v)
{
  double wavelength = (planck*lightspeed) / sqrt(charge*v*(2*mass*lightspeed2+charge*v));
}



/* Calculates the interaction parameter */ 
double calcsigma(double wavelength){
  double sigma = (2*pi*mass*wavelength*charge) / pow(planck, 2); 
}



