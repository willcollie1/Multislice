/* This is the header file for "Multilsice.c" containing all the definitions and function declarations */

#define planck 6.626176e-34 /*Plancks constant */
#define mass 9.10938356e-31 /* Electron Mass */
#define charge 1.60217e-19 /*Electron Charge */
#define pi 3.14159265359  /* Constant pi */ 
#define lightspeed 2.99792e8 /* Speed of Light */
#define lightspeed2 8.98755e16 /* Speed of Light Squared */ 
#define hartree 4.35974e-18 /* How mant hartrees in a Joule */
#define angstrom 1e-10 /* Angstrom in metres */
#define abberation 1e-6 /*Coefficient of Spherical Abberation (Microscope constant) */ 

int i,j,loop1,tilt, height,width,depth,size,wavefunctionsize,fy, tx, ty, fx,unitcellcode,n[2],temp1,sideview; 
double a,b,c,v,total,total1,dx,dy,dz,*x,*y,*V,abbfunction,*k_x2,*k_y2,wavelength,sigma,df,aperture,tiltx,tilty,*k,tester,kmax,temp;
double complex *T,*P,*wavefunction,*PSF,var,*P2,tl, tr, bl, br;
char rec;

/* FFTW declarations */ 
fftw_complex *in, *out;
fftw_plan plan;
fftw_plan planinverse;


  /* Getting the user input */ 
  void getuserinput(void){
    printf("Unit cell code (1-Rectangular, 2-Hexagonal):");
    scanf("%d",&unitcellcode); 
    printf("Lattice parameter a (Angstroms):");
    scanf("%lf",&a);
    printf("Lattice parameter b (Angstroms):");
    scanf("%lf",&b); 
    printf("Lattice Parameter c (Angstroms):");
    scanf("%lf",&c); 
    printf("Width:");
    scanf("%d",&width);  
    printf("Height:");
    scanf("%d",&height); 
    printf("Depth:");
    scanf("%d",&depth); 
    printf("Microscope Voltage (volts):");
    scanf("%lf",&v);
  }

  /* Output a file containing the array width and height */
  int outputparameters(void){
    FILE* output2;
    output2 = fopen("size.txt", "w");
    fprintf(output2, "%d \n%d\n", width,height);
    fclose(output2);
  }


  /* Relativistic De broglie wavelength */ 
  double calcwavelength(double v){
    wavelength = (planck*lightspeed) / sqrt(charge*v*(2*mass*lightspeed2+charge*v));
  }

  /* Interaction parameter */ 
  double calcsigma(double wavelength, double v){
    sigma = ((2*pi)/(wavelength*v))*((mass*lightspeed2+charge*v)/(2*mass*lightspeed2+charge*v)); 
  }


  /* Dynamically allocating memory */
  void allocatememory(){
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
    k = malloc((wavefunctionsize)*sizeof(double));
  }


  /* Reading in the ordered potential file for standard view */ 
  int readinputfile(){
    FILE *potential = fopen("orderedfile.txt", "r");
      for (i=0; i<size; i++){
        fscanf(potential, "%lf %lf %*f  %lf", &x[i], &y[i], &V[i]);
       }
    fclose(potential);
  }

   /* Reading in the ordered potential file for side view*/ 
  int readinputfilesideview(){
    FILE *potential = fopen("orderedfile.txt", "r");
      for (i=0; i<size; i++){
        fscanf(potential, "%*f %lf %lf  %lf", &y[i], &x[i], &V[i]);
       }
    fclose(potential);
  }

  /* Calculating where to start the multislice loop (optimisation) */ 
  double whereareatoms(){
    for(i = 0; i<size; i++){
      total = abs(V[i]) + total; 
    }
    for(i = 0.9*size; i<size; i++){
      total1 = abs(V[i]) + total1;  
    }

    if (total1 > 0.1*total){
      loop1 = 0.6*depth;
    }
    else{
      loop1 = 0;
    }
  }



  /* Spatial Frequency Calculation (x direction) */ 
  double calck_x2(int width, double dx,double*x){
    for (i=0; i<wavefunctionsize;i++){
      if((x[i]/dx) > (width/2)){
          k_x2[i] = pow((((x[i]/dx)-(width+1)) / (width*dx)),2);
      }
      else{
        k_x2[i] = pow((x[i]/dx) / (width*dx),2);
      }
    }
  } 



  /* Spatial frequency Calculation (y direction) */
  double calck_y2(int height, double dy,double*y){
    for (i=0; i<wavefunctionsize;i++){
      if((y[i]/dy) > (height/2)){
         k_y2[i] = pow((((y[i]/dy)- (height+1)) / (height*dy)),2);
      }
      
      else{
          k_y2[i] = pow((y[i]/dy) / (height*dy),2);
      }
    }
  }

  /* Hexagonal lattice SPX */
  double calck_x2H(int width, double dx,double*x){
   for(i=0; i<wavefunctionsize;i++){
     k_x2[i] = pow((x[i]/dx) / (width*dx),2);
   }
 }

 /* Heagonal lattice SPY */ 
 double calck_y2H(int height, double dy,double*y){
   for(i=0; i<wavefunctionsize;i++){
     k_y2[i] = pow((y[i]/dy) / (height*dy),2);
   }
 }


  /* Allows us to sample from the unit cell using x and y as if it were a rectangular unit cell */ 
  double hexagonalsampling(){

    for (i=0; i<size; i++){
      var = (-y[i])*(sin(pi/6)/cos(pi/6));
      y[i] = y[i] / cos(pi/6); 
      x[i] = x[i] - var;
  
      x[i] = x[i] / (dx);
      y[i] = y[i] / (dy);
    }
  }

  /* Interpolate to get potential (P2) for these non integer coordinates */
  double interpolation(){
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
  }

  /* Freeing memory */ 
  void destroymem(void){
   fftw_destroy_plan(plan);
   fftw_destroy_plan(planinverse);
   fftw_free(in); fftw_free(out); 
   free(T), free(P), free(wavefunction),free(x), free(y), free(V), free(k_x2), free(k_y2), free(PSF), free(P2); 
  }






