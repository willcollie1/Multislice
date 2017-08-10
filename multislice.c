/* This program uses a multislice technique to simulate the wavefunction for electrons in a TEM as they pass through a thin sample.
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
 
  /* Gets the required input parameters from the user (and corrects units) */  
  getuserinput();
   
  /* Converting to S.I units */ 
  a = a * angstrom;   
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
  allocatememory();

  /* Reading in the correctly ordered potential file "ordered.txt" */ 
  readinputfile();

  /* Converting potentials from Hartrees to Volts (energy to volltage) and positions into metres */ 
  for (i=0; i<size; i++){
    V[i] = V[i] * hartree/charge;  
    x[i] = x[i] * dx;
    y[i] = y[i] * dy;
  }

  /* Working out where the atoms are located (where to start the j loop)  */ 
  whereareatoms();

  /* Calculating the Spatial frequency squared (reciprical lattice) */ 
  calck_x2(width,dx,x);
  calck_y2(height,dy,y);

  for (i=0; i<wavefunctionsize; i++){
    k[i] = pow((k_x2[i]+k_y2[i]),0.5);
  }
  kmax = pow((pow((width/(2*a)),2) + pow((height/(2*b)),2)),0.5);

  /* If the unit cell is Hexagonal, follow this step of the loop */ 
  if(unitcellcode == 2){
   hexagonalsampling();
   interpolation();
  }

  /* If the unit cell is rectangular */ 
  if(unitcellcode == 1){
    for(i=0;i<size;i++){
      P2[i] = V[i];
    }
  }

  /* Propogation Function (k space) No specimin tilt*/
  for (i=0; i<wavefunctionsize; i++){
    P[i] = cexp(-I*pi*dz*wavelength*(k_x2[i]+k_y2[i]));
  }


  /* Allowing for Specimin tilt */ 
  if(tilt == 1){
    for (i=0; i<wavefunctionsize; i++){
    P[i] = cexp((-I*pi*dz*wavelength*(k_x2[i]+k_y2[i]))+(2*pi*I*dz*((pow(k_x2[i],0.5)*tan(tiltx))+k_y2[i]*tan(tilty))));
  }
  }
  
  /* Calculating the point spread function for the objective lens */ 
  df = pow(1.5*abberation*wavelength,0.5); /* Defocus */ 
  aperture = pow(((4*wavelength)/abberation),0.25); /* Aperture size */

  


  for (i=0; i<wavefunctionsize; i++){
    if(k[i] <= aperture){
      abbfunction = (pi*wavelength*(k_x2[i]+k_y2[i]))*((0.5*abberation*pow(wavelength,2)*(k_x2[i]+k_y2[i]))-df); 
      PSF[i] = cexp(-I*abbfunction);
    }
    else{
      abbfunction = 0; 
      PSF[i] = cexp(-I*abbfunction);
    }
  }

  
  /* Initialising FFTW */
  n[2];
  n[0] = width;
  n[1] = height;
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * wavefunctionsize);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * wavefunctionsize);
  plan = fftw_plan_dft_2d(height, width, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  planinverse = fftw_plan_dft_2d(height, width, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);


    /* Initialising the wavefunction */
  for (i=0; i<wavefunctionsize; i++){
    wavefunction[i] = 1.;
  }
  

  /* Multislice loop */ 
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

  /* Convolution with Objective lens */
  for (i=0; i<wavefunctionsize; i++){
    in[i] = wavefunction[i];
  }  
  
  fftw_execute(plan);
 
  for (i=0; i<wavefunctionsize; i++){
    in[i] = out[i]*PSF[i];
  }
  
  fftw_execute(planinverse);

  for (i=0; i<wavefunctionsize; i++){
    wavefunction[i] = out[i] / wavefunctionsize;  
  }



  /* Convergence Test */ 
  for(i=0;i<wavefunctionsize;i++){
    tester = pow(cabs(wavefunction[i]),2) + tester;
  }

  if(tester/wavefunctionsize < 0.9){
    printf("The total integrated intensity is %lf which is less than 0.9. There is likely to be a problem with this simulation \n", tester/wavefunctionsize);
  }

  if(tester/wavefunctionsize > 1.0){
    printf("The total integrated intensity is %lf which is greater than 1. There is likely to be a problem with this simulation \n", tester/wavefunctionsize);
  }

  else{
    printf("The total integrated intensity is: %lf - Test Passed \n", tester/wavefunctionsize);
  }

 
  /* Outputting the final image intensity (square modulus of wavefunction) */ 
  FILE* output;
  output = fopen("output.txt", "w");
  for (i=0; i<wavefunctionsize; i++){   
    fprintf(output, "%.20lf \n", pow(cabs(wavefunction[i]),2)); 
  }
  fclose(output);
   

  /* Free the dynamic arrays */
  destroymem(); 

}