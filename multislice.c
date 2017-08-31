#include <stdio.h>
#include <stdlib.h>
#include <math.h>    /* Use of pow function */ 
#include <complex.h> /* Makes the use of complex numbers much easier */
#include <fftw3.h>  /* Used for FFTs */
#include "multislice.h" /* Contains definitions and some function declarations */ 
 
int main()
{

  /* Gets the required input parameters from the user (in angstroms) */  
  getuserinput();

  /* Converting lattice parameters to S.I units */ 
  a = a * angstrom;   
  b = b * angstrom; 
  c = c * angstrom; 

  /* Swapping a with c and width with depth, consistent with a side view */
  if(view == 1){
  temp = a;
  a = c;
  c = temp;

  temp1 = width;
  width = depth;
  depth = temp1;
  }

  /* Swapping a with c and width with depth, consistent with a side view */
  if(view == 2){
  temp = b;
  b = c;
  c = temp;

  temp1 = height;
  height = depth ;
  depth = temp1;
  }
    
  /* Size of an individual voxel */
  dx = (a / width);
  dy = (b / height);
  dz = (c / depth);  /* Slicethickness */ 
  size = height*width*depth; /* The number of potential readings in file */
  wavefunctionsize = height*width; /* Number of points wavefunction is evaluated at */

  /* Outputs the wavefunction parameters to the file "size.text" to be used by the python script */ 
  outputparameters();

  /* Calculate the Relativistic De Broglie wavelength and Electron Interaction parameter */ 
  wavelength = calcwavelength(v);
  sigma = calcsigma(wavelength,v);
  
  /* Dynamically Allocates Memory for the arrays used */   
  allocatememory();

  /* Reads in the correctly ordered potential file "ordered.txt" */ 
  if(view == 0){
    readinputfile(); 
  }

  if(view == 1){
    readinputfilesideview();
  }

  if(view == 2){
    readinputfiletopview();
  }

  /* Converts potentials from Hartrees to Volts (energy to volltage),and positions into metres */ 
  for (i=0; i<size; i++){
    V[i] = V[i] * hartree/charge;  
    x[i] = x[i] * dx;
    y[i] = y[i] * dy;
  }

  /* Working out where the atoms are located (maximising efficiency) */ 
  whereareatoms();

  /* Calculating the Spatial Frequency squared (reciprical lattice) */ 
  calck_x2(width,dx,x);
  calck_y2(height,dy,y);

  /* If the unit cell is Hexagonal, follow this step of the loop */ 
  if(unitcellcode == 1){
   calck_x2H(width,dx,x);
   calck_y2H(height,dy,y);
   hexagonalsampling();
   interpolation();
  }

  /* If the unit cell is rectangular then no remapping/interpolation is needed */ 
  if(unitcellcode == 0){
    for(i=0;i<size;i++){
      P2[i] = V[i];
    }
  }

  /* Propogation Function */
  for (i=0; i<wavefunctionsize; i++){
    P[i] = cexp(-I*pi*dz*wavelength*(k_x2[i]+k_y2[i]));
  }

  /* Propogation Function including specimin tilt */ 
  if(tilt == 1){
    for (i=0; i<wavefunctionsize; i++){
      P[i] = cexp((-I*pi*dz*wavelength*(k_x2[i]+k_y2[i]))+(2*pi*I*dz*((pow(k_x2[i],0.5)*tan(tiltx))+pow(k_y2[i],0.5)*tan(tilty))));
    }
  }
  
  /* Calculating the Transfer Function of the objective Lens */ 
  df = pow(1.5*abberation*wavelength,0.5); /* Optimum defocus */ 
  aperture = pow(((4*wavelength)/abberation),0.25); /* Optimum Aperture size */
  for (i=0; i<wavefunctionsize; i++){
    k = pow((k_x2[i]+k_y2[i]),0.5);
    if(k <= aperture){
      abbfunction = (pi*wavelength*(k_x2[i]+k_y2[i]))*((0.5*abberation*pow(wavelength,2)*(k_x2[i]+k_y2[i]))-df); 
      PSF[i] = cexp(-I*abbfunction);
    }
    else{
      abbfunction = 0; 
      PSF[i] = cexp(-I*abbfunction);
    }
  }

  /* Initialising FFTW arrays */
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

  if(tester/wavefunctionsize > 1.00001){
    printf("The total integrated intensity is %lf which is greater than 1. There is likely to be a problem with this simulation \n", tester/wavefunctionsize);
  }

  else{
    printf("The total integrated intensity is: %lf - Test Passed \n", tester/wavefunctionsize);
  }

  /* Output the image intensity to "output.txt" */ 
  FILE* output;
  output = fopen("output.txt", "w");
  for (i=0; i<wavefunctionsize; i++){   
    fprintf(output, "%.20lf \n", pow(cabs(wavefunction[i]),2)); 
  }
  fclose(output);

  /* Free the dynamic arrays */
  destroymem(); 

}