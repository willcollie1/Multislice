#This program uses the output file from "multislice.c" to create a greyscale image of the unit cell.

import subprocess
import os #Run linux commands from Python program
from PIL import Image #Image Processing
import numpy as np #Poisson Distribution


#Unzip the potential file (if it isn't already unzipped)
unzipping = os.system("bzip2 -dk 'Example Potential file .pot_fmt'.bz2")

#Uses Linux sort to reorder the input file to row major format.
reordering = os.system("sort -n -k3,3 -k2,2 -k1,1 < graphene.txt > output1.txt")  

#Runs the multislice C program
var = os.system("./multislice")

#Reads in array parameters from multislice
parameters = open("size.txt").readlines()
height = int(parameters[0])
width = int(parameters[1]) 
size = height*int(width)

#Read in intensity values from multislice.c and create a list
List = open("output.txt").readlines()

edose = 152665 

#Drawing a random number from a Poission distribution centred on edose
for j in range (0,size):
  List[j] = np.random.poisson(float(List[j])*edose) #generalise this

#Normalise intensity values so that they fill the range 0-255 (8 bit binary)
maximum = max(List)
minimum = min(List)
for i in range (0,size):
    List[i] =  255*(float(List[i])-float(minimum))/(float(maximum) - float(minimum))
    List[i] = int(List[i])

#Create greyscale image with the given intensity values.
newim = Image.new('L', (height,width))
pixels1 = newim.load()
for i in range (0,(height-1)):
    for j in range(0,(width-1)):
        l = List[i+height*j] 
 
        pixels1[i,j] = (l) 

newim.save("newim", "JPEG")