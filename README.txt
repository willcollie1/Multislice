Simulating Transmission Electron Microscope Micrographs using the Multislice technique. 

This software requires a four column text file describing the potential throughout the specimin that is to be imaged. The specimin must be a thin sample for the program to work correctly. To obtain a description of the potential within the specimin DFT software such as Castep can be used.



multislice.py: This is the main controller of the software. Run this and it will open the GUI, take user input, compile the C program, run the c program and do all the image processing!


multislice.c: This is where the bulk of the calculation takes place. This program reads in the potential data, undergoes a series of calcululations including mutiple fast Fourier transforms before finally outputting a file containing the absolute square of the wave function at each position.

multislice.h: Header file for multislice.c. All functions used in c code are defined here. 


GUI.py: Program that controls the GUI taking user input. 
