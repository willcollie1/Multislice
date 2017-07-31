/* This program uses the "Multislice Algorithm" to simulate the wavefunction for a TEM electron as it passes through a sample.
The multislice method reduces to a succession of transmission and propagation operations with a Fast Fourier Transform (FFT) in 
between each. FFTs are calculated using the FFTW libary. The square modulus of the exit wavefunction is (once abberations are accounted for) outputted to the file called "output.txt", which is used by the python script "program.py" to create an image of the specimin */ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>    /* Use of pow function */ 
#include <complex.h> /* Makes the use of complex numbers much easier */
#include <fftw3.h>  /* Used for FFTs */
#include "multislice.h" /* Contains definitions and some function declarations */ 
#include <limits.h>


int main()
{
  
  /* Lattice Parameters (Angstroms) */ 
  a = 14.794;  
  b = 14.794;
  c = 28.00;  
  
  /* Microscope Voltage (Volts) */ 
  printf("Microscope Voltage (volts):");
  scanf("%lf",&v);
  
  /* Number of voxels */ 
  width = 120;
  height = 120;
  depth = 225;
  
  /* <a b c> Lattice Parameters in Angstroms */
  a = a * angstrom; /* Working with S.I units throughout */  
  b = b * angstrom; 
  c = c * angstrom;
  
  /* Size of an individual "voxel" */
  dx = (a / width);
  dy = (b / height);
  dz = (c / depth);  /* Slicethickness */ 
  size = height*width*depth; /* The number of potential readings in file */
  wavefunctionsize = height*width; /* Number of points wavefunction is evaluated at */

  /* Relativistic de broglie wavelength and Interaction parameter (sigma) */ 
  wavelength = calcwavelength(v);
  sigma = calcsigma(wavelength,v);
  
  /* Output the wavefunction parameters to the file "size.text" (used in python script) */ 
  outputparameters();

 /* Dynamic Memory Allocation */   
  x  = malloc((size)*sizeof(double)); /* X postion */
  y  = malloc((size)*sizeof(double)); /* Y position */
  V  = malloc((size)*sizeof(double)); /* Potential */ 
  P2 = malloc((size)*sizeof(double complex)); /*Interpolated Potential */ 
  T  = malloc((wavefunctionsize)*sizeof(double complex)); /* Transmission Function */
  P  = malloc((wavefunctionsize)*sizeof(double complex)); /* Propogation Function */
  wavefunction = malloc((wavefunctionsize)*sizeof(double complex)); /* Wavefunction */
  PSF = malloc((wavefunctionsize)*sizeof(double complex)); /* Optical Transfer function */
  k_x2 = malloc((wavefunctionsize)*sizeof(double)); /* Spatial Frequency */
  k_y2 = malloc((wavefunctionsize)*sizeof(double)); /* Spatial Frequency */


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

  /* Deciding where multislice loop should start */ 
  if (total1 > 0.1*total){
    loop1 = 0.6*depth;
  }
  else{
    loop1 = 0;
  }

  /* Calculating the Spatial frequency squared (reciprical lattice) */ 
  calck_x2(width,dx,x);
  calck_y2(height,dy,y);

  /* Allows us to sample from the unit cell using x and y as if it were a rectangular unit cell */
  for (i=0; i<size; i++){
    var = (-y[i])*(sin(pi/6)/cos(pi/6));
    y[i] = y[i] / cos(pi/6); 
    x[i] = x[i] - var;
  
    x[i] = x[i] / (dx);
    y[i] = y[i] / (dy);
  }

  /* Interpolate to get potential (P2) for these non integer coordinates */
  for (j = loop1; j < depth; j++){
    for (i=0; i<wavefunctionsize; i++){
      fx = floor(x[i+wavefunctionsize*j]);
      fy = floor(y[i+wavefunctionsize*j]);
      tx = fx + 1;
      ty = fy + 1;

      tl = (1-(x[i + wavefunctionsize*j] - floor(x[i+wavefunctionsize*j]))) * (y[i + wavefunctionsize*j]-floor(y[i + wavefunctionsize*j]));
      tr = (x[i + wavefunctionsize*j]- floor(x[i + wavefunctionsize*j])) * (y[i + wavefunctionsize*j] - floor(y[i + wavefunctionsize*j]));
      bl = (1 - (x[i + wavefunctionsize*j] - floor(x[i + wavefunctionsize*j]))) * (1 - (y[i + wavefunctionsize*j] - floor(y[i + wavefunctionsize*j])));
      br = (x[i + wavefunctionsize*j] - floor(x[i + wavefunctionsize*j])) * (1 - (y[i + wavefunctionsize*j] - floor(y[i + wavefunctionsize*j])));

      P2[i+wavefunctionsize*j] = (tl * V[wavefunctionsize*j + ty*width + fx ]) + (tr * V[wavefunctionsize*j + ty*width + tx ]) + (bl * V[wavefunctionsize*j + width*fy + fx ]) + (br * V[wavefunctionsize*j + fy*width + tx]);

  } 
}

  /* Propogation Function (k space) */
  for (i=0; i<wavefunctionsize; i++){
    P[i] = cexp(-I*pi*dz*wavelength*(k_x2[i]+k_y2[i]));
  }

  /* Initialising the wavefunction */
  for (i=0; i<wavefunctionsize; i++){
    wavefunction[i] = 1.;
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
      T[i] = cexp(-I*sigma*P2[i+j*wavefunctionsize]*dz);  
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


  /* Accounting for Abberations (spherical and defocus) */
  for (i=0; i<wavefunctionsize; i++){
    in[i] = wavefunction[i];
  }
    
  fftw_execute(plan); 
 
  for (i=0; i<wavefunctionsize; i++){
    abbfunction = 2*pi/wavelength; /* (pi*wavelength*(k_x2[i]+k_y2[i]))*((0.5*abberation*pow(wavelength,2)*(k_x2[i]+k_y2[i]))-df); */ 
    PSF[i] = cexp(-I*abbfunction);
   }

  for (i=0; i<wavefunctionsize; i++){
    in[i] = out[i]*PSF[i];
  }
    
  fftw_execute(planinverse);

  for (i=0; i<wavefunctionsize; i++){
    wavefunction[i] = out[i] / wavefunctionsize;  
  }

    
  /* Outputting the final image intensity (square modulus of wavefunction) */ 
  FILE* output;
  output = fopen("output.txt", "w");
  for (i=0; i<wavefunctionsize; i++){   
    fprintf(output, "%.20lf \n", pow(cabs(wavefunction[i]),2)); 
  }
  fclose(output);
   

  /* Free the dynamic arrays */
  fftw_destroy_plan(plan);
  fftw_destroy_plan(planinverse);
  fftw_free(in); fftw_free(out); 
  free(T), free(P), free(wavefunction),free(x), free(y), free(V), free(k_x2), free(k_y2), free(PSF), free(P2); 

}