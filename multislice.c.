/* This program uses the "Multislice Algorithm" to simulate the wavefunction for a TEM electron as it passes through a sample.
The multislice method reduces to a succession of transmission and propagation operations with a Fast Fourier Transform (FFT) in 
between each. FFTs are calculated using the FFTW libary. The final wavefunction is outputted to the file called "output.txt" 
which is used by the python script "program.py" to create an image of the specimin */ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>    
#include <complex.h> /* Makes the use of complex numbers much easier */
#include <fftw3.h>  /* Used for FFTs */
#include "multislice.h" /* Contains definitions and some function declarations */ 


int main()
{
  
  /* Get input parameters from  the user */ 
  printf("Lattice Parameter (a) in Angstroms:");
  scanf("%lf",&a);

  printf("Lattice Parameter (b) in Angstroms:");
  scanf("%lf",&b);

  printf("Lattice parameter (c) in Angstroms:");
  scanf("%lf",&c);

  printf("Microscope Voltage (volts):");
  scanf("%lf",&v);

  printf("Number of Voxels in x direction:");
  scanf("%d",&width);

  printf("Number of Voxels in y direction:");
  scanf("%d",&height);

  printf("Number of Voxels in z direction:");
  scanf("%d",&depth);
  

  /* <a b c> Lattice Parameters in Angstroms */
  a = a * angstrom; /* 2.54 */
  b = b * angstrom; /* 4.36 */
  c = c * angstrom; /* 28.00 */
  

  /* Size of an individual "voxel" */
  double dx = (a / width);
  double dy = (b / height);
  double dz = (c / depth);  /* Slicethickness */ 
  int size = height*width*depth; /* The number of potential readings in file */
  int wavefunctionsize = height*width; /* Number of points wavefunction is evaluated at */

      
  /* Output the wavefunction parameters to the file "size.text" (used in python script) */ 
  outputparameters();

  
  /*Allocating memory to store data from potential file in three seperate arrays */ 
  double *x, *y, *V; 
  x  = malloc((size)*sizeof(double));
  y  = malloc((size)*sizeof(double));
  V  = malloc((size)*sizeof(double));

 
  /* Reading in the correctly ordered potential file "output1.txt" */ 
    FILE *potential = fopen("output1.txt", "r");
  for (i=0; i<size; i++){
    fscanf(potential, "%lf %lf %*f  %lf", &x[i], &y[i], &V[i]);
  }
  fclose(potential);

  
  /* Converting potentials from Hartrees to Volts (energy to volltage) and positions into metres */ 
  for (i=0; i<size; i++){
    V[i] = V[i] * hartree/charge;  
    x[i] = x[i] * dx;
    y[i] = y[i] * dy;
  }
  
    /* If the potential is concentrated at edge of the unit cell then we can skip to where the atoms are  */ 
  for(i = 0; i<size; i++){
    total = abs(V[i]) + total; /* Sum of all the potentials */
  }
  for(i = 0.9*size; i<size; i++){
    total1 = abs(V[i]) + total1;  /* Sum of all potentials in final 10% */
  }
  

  /* Spatial frequency squared (k^2) */ 
  double *k_x2,*k_y2;
  k_x2 = malloc((size)*sizeof(double));
  k_y2 = malloc((size)*sizeof(double));
  
  for (i=0; i<size;i++) {
    if (x[i] < (width/2)) {
      k_x2[i] = pow((x[i]) / (width*pow(dx,2)),2);
    }
    else{
       k_x2[i] = pow((x[i] - width) / (width*pow(dx,2)),2);
    }
  }
  
   for (i=0; i<size;i++) {
    if (y[i] < (height/2)) {
      k_y2[i] = pow((y[i]) / (height*pow(dy,2)),2);
    }
    else{
       k_y2[i] = pow((y[i] - height) / (height*pow(dy,2)),2);
    }
   }

  /* Relativistic de broglie wavelength and Interaction parameter (sigma) */ 
  double wavelength = calcwavelength(v);
  double sigma = calcsigma(wavelength,v);
  
  /* Transmission Function */
  double complex *T; 
  T  = malloc((size)*sizeof(double complex));
  for (i=0; i<size; i++){
    T[i] = cexp(I*sigma*V[i]*dz);  
  }
  

  /* Propogation Function in k space */
  double complex *P; 
  P  = malloc((size)*sizeof(double complex));
  for (i=0; i<size; i++){
    P[i] = cexp(-I*pi*dz*wavelength*(k_x2[i]+k_y2[i]));
  }

  /* Initialising and assigning memory for wavefunction and wavefunction_next */
  double complex *wavefunction; 
  wavefunction = malloc((wavefunctionsize)*sizeof(double complex));

  /* Set initial wavefunction equal to 1 */
  for (i=0; i<wavefunctionsize; i++){
    wavefunction[i] = 1.;
  }
  
  /* Defining j loop boundaries (depends where the atoms are located) */ 
  if (total1 > 0.1*total){
    loop1 = 0.6*depth;
  }
  else{
    loop1 = 0;
  }
  
  
  /* Initialising and allocating memory for FFT function */
  fftw_complex *in, *out;
  fftw_plan plan;
  fftw_plan planinverse;
  int n[2];
  n[0] = width;
  n[1] = height;
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * wavefunctionsize);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * wavefunctionsize);
  plan = fftw_plan_dft_2d(height, width, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  planinverse = fftw_plan_dft_2d(height, width, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  
   
 /************************************************** Multislice Iteration loop **********************************************************/
  
  /* Looping over the number of slices */ 
  for(j=loop1; j<depth; j++){
    
    /* Multiplying the Wavefunction by the Transmission function */ 
    for (i=0; i<wavefunctionsize; i++){
      wavefunction[i] = T[i+(wavefunctionsize*j)]*wavefunction[i];
    }

    /* Taking the 2D Fourier transform of this product */ 
    for (i=0; i<wavefunctionsize; i++){
      in[i] = wavefunction[i];
    }
    
    fftw_execute(plan);

    
    /* Multiplying by the Propogation function and taking the Inverse FT  */ 
    for (i=0; i<wavefunctionsize; i++){
      wavefunction[i] = out[i]*P[i+(wavefunctionsize*j)];
    }
    
    for(i=0; i<wavefunctionsize; i++){
    in[i] = wavefunction[i];
    }
    
    fftw_execute(planinverse);
  
    for (i=0;i<wavefunctionsize;i++){
      wavefunction[i] = out[i] / wavefunctionsize ; /* Scaling factor */
    }

  } 
 
/*************************** End of Multislice loop **************************************/ 

  
    /* Fourier Transform of Exit Wavefunction */ 
    for (i=0; i<wavefunctionsize; i++){
      in[i] = wavefunction[i];
    }
    
    fftw_execute(plan); 
  
    /* Calculating the point spread function (PSF) of the objective lens */
    double complex *PSF;
    double *susceptability;
    susceptability = malloc((wavefunctionsize)*sizeof(double));
    PSF = malloc((wavefunctionsize)*sizeof(double complex));

    for (i=0; i<wavefunctionsize; i++) {
      susceptability[i] = (2*pi/wavelength)*((0.25*abberation*pow(wavelength,4)*(pow(pow((k_x2[i]+k_y2[i]),0.5),4)))-(0.5*df*pow(wavelength,2)*(k_x2[i]+k_y2[i])));
      PSF[i] = cexp(-I*susceptability[i]);
    }
    

    /* FT of exit wavefunction multiplied by the PSF */ 
    for (i=0; i<wavefunctionsize; i++) {
      wavefunction[i] = out[i]*PSF[i];
    }
    
    /* Inverse FT and calculate square modulus to get final image intensity */
    for (i=0; i<wavefunctionsize; i++){
      in[i] = wavefunction[i];
    }
    
    fftw_execute(planinverse);

    for (i=0; i<wavefunctionsize; i++){
      wavefunction[i] = out[i] / wavefunctionsize;  /* Scaling factor */ 
    }
  

    
  /* Output the Intensity at each point to a file called "output.txt" */ 
  FILE* output;
  output = fopen("output.txt", "w");
  for (i=0; i<wavefunctionsize; i++){   
    fprintf(output, "%.20lf \n", pow(cabs(wavefunction[i]),2)); /* Output the square modulus of the wavefunction */ 
  }
  fclose(output);
   

   /* Free the dynamic arrays */
   fftw_destroy_plan(plan);
   fftw_destroy_plan(planinverse);
   fftw_free(in); fftw_free(out); 
   free(T), free(P), free(wavefunction),free(x), free(y), free(V), free(k_x2), free(k_y2);

}