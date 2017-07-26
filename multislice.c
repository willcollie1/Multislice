/* This program uses the "Multislice Algorithm" to simulate the wavefunction for a TEM electron as it passes through a sample.
The multislice method reduces to a succession of transmission and propagation operations with a Fast Fourier Transform (FFT) in 
between each. FFTs are calculated using the FFTW libary. The square modulus of the exit wavefunction is (once abberations are accounted for) outputted to the file called "output.txt", which is used by the python script "program.py" to create an image of the specimin */ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>    /* Use of pow function */ 
#include <complex.h> /* Makes the use of complex numbers much easier */
#include <fftw3.h>  /* Used for FFTs */
#include "multislice.h" /* Contains definitions and some function declarations */ 


int main()
{
  
  /* Lattice Parameters */ 
  a = 2.54;
  b = 4.36;
  c = 28.00;
  
  /* Microscope voltage */ 
  printf("Microscope Voltage (volts):");
  scanf("%lf",&v);
  
  /* Number of voxels */ 
  width = 24;
  height = 40;
  depth = 240;
  
  /* <a b c> Lattice Parameters in Angstroms */
  a = a * angstrom; /* 2.54 */
  b = b * angstrom; /* 4.36 */
  c = c * angstrom; /* 28.00 */
  
  /* Size of an individual "voxel" */
  dx = (a / width);
  dy = (b / height);
  dz = (c / depth);  /* Slicethickness */ 
  size = height*width*depth; /* The number of potential readings in file */
  wavefunctionsize = height*width; /* Number of points wavefunction is evaluated at */
  
  /* Output the wavefunction parameters to the file "size.text" (used in python script) */ 
  outputparameters();

  /*Allocating memory to store data from potential file in three seperate arrays */   
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
  
  /* Calculating the Spatial frequency squared (reciprical lattice) */ 
  k_x2 = malloc((wavefunctionsize)*sizeof(double));
  k_y2 = malloc((wavefunctionsize)*sizeof(double));
  calck_x2(width,dx,x);
  calck_y2(height,dy,y);
  
  /* Relativistic de broglie wavelength and Interaction parameter (sigma) */ 
  wavelength = calcwavelength(v);
  sigma = calcsigma(wavelength,v);

  /* Dynamic memory allocation for Transmission and Propogation function and wavefunction */ 
  T  = malloc((wavefunctionsize)*sizeof(double complex));
  P  = malloc((wavefunctionsize)*sizeof(double complex));
  wavefunction = malloc((wavefunctionsize)*sizeof(double complex));


  /* Propogation Function (k space) */
  for (i=0; i<wavefunctionsize; i++){
    P[i] = cexp(-I*pi*dz*wavelength*(k_x2[i]+k_y2[i]));
  }

  /* Initialising the wavefunction */
  for (i=0; i<wavefunctionsize; i++){
    wavefunction[i] = 1.;
  }
  
  /* Deciding where multislice loop should start */ 
  if (total1 > 0.1*total){
    loop1 = 0.6*depth;
  }
  else{
    loop1 = 0;
  }
  
  
  /* Initialising and allocating memory for FFTW function */
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
  
  
  /**** Multislice Loop ****/ 
  for(j=loop1; j<depth; j++){
    for (i=0; i<wavefunctionsize; i++){
      T[i] = cexp(I*sigma*V[i+j*wavefunctionsize]*dz);  
    }
    for (i=0; i<wavefunctionsize; i++){
      in[i] = T[i]*wavefunction[i];
    }    
    
    fftw_execute(plan);
 
    for (i=0; i<wavefunctionsize; i++){
      in[i] = out[i]*P[i];
    }
     
    fftw_execute(planinverse);
  
    for (i=0;i<wavefunctionsize;i++){
      wavefunction[i] = out[i] / wavefunctionsize ; 
    }
  } 
  /****End of Multislice Loop****/ 


  /* Accounting for Abberations */
  abbfunction = malloc((wavefunctionsize)*sizeof(double)); /* Abberation function */ 
  PSF = malloc((wavefunctionsize)*sizeof(double complex)); /* Point spread of "optical transfer function" */

  for (i=0; i<wavefunctionsize; i++){
    in[i] = wavefunction[i];
  }
    
  fftw_execute(plan); 
 
  for (i=0; i<wavefunctionsize; i++){
    abbfunction[i] = 2*pi/wavelength; /* (pi*wavelength*(k_x2[i]+k_y2[i]))*((0.5*abberation*pow(wavelength,2)*(k_x2[i]+k_y2[i]))-df); */ 
     PSF[i] = cexp(-I*abbfunction[i]);
   }

  for (i=0; i<wavefunctionsize; i++){
    in[i] = out[i]*PSF[i];
  }
    
  fftw_execute(planinverse);

  for (i=0; i<wavefunctionsize; i++){
    wavefunction[i] = out[i] / wavefunctionsize;  
  }
    
  /* Outputting the final intensity */ 
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
  free(T), free(P), free(wavefunction),free(x), free(y), free(V), free(k_x2), free(k_y2), free(abbfunction), free(PSF);

}