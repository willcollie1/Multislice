Simulating Transmission Electron Microscope Micrographs using the Multislice technique. 

This software requires a four column text file describing the potential throughout a thin specimen that is to be imaged. To obtain a description of the potential throughout the specimen, DFT software such as CASTEP can be used.

multislice.py: This is the main controller of the software. Run this and it will open the GUI, take user input, compile the C program, run the C program and do all the image processing!

multislice.c: This is where the bulk of the calculation takes place. This program reads in the potential data, undergoes a series of calcululations including mutiple fast Fourier transforms before finally outputting a file containing the absolute square of the wave function at each position.

multislice.h: Header file for multislice.c. All functions used in C code are defined here. 

GUI.py: Script for controlling the GUI that taking user input. 
