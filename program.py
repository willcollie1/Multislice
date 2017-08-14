#This program uses the output file from "multislice.c" to create a greyscale image of the unit cell.

import subprocess
import os #Run linux commands from Python program
from PIL import Image #Image Processing
import numpy as np #Poisson Distribution


compiling = os.system("gcc multislice.c -o multislice -lfftw3 -lm")

filename = raw_input('Enter a filename:' )
sideview = raw_input(' Normal view (0), Side view(1):' )

#Uses Linux sort to reorder the input file to row major format.
string = "sort -n -k3,3 -k2,2 -k1,1 < "
string1 = "> orderedfile.txt"
string2 = "sort -n -k1,1 -k2,2 -k3,3 < "


if(int(sideview) == 0):
    cmd = string + filename + string1
    reordering = os.system(cmd)
#else:
 #   cmd = string2 + filename + string1
  #  reordering = os.system(cmd)


#Runs the multislice C program
var = os.system("./multislice")

#Reads in array parameters from multislice
parameters = open("size.txt").readlines()
height = int(parameters[0])
width = int(parameters[1]) 
size = height*int(width)

#Read in intensity values from multislice.c and create a list
List = open("output.txt").readlines()

edose = 200000

#Drawing a random number from a Poission distribution centred on edose
for j in range (0,size):
  List[j] = np.random.poisson(float(List[j])*edose) #generalise this

#Normalise intensity values so that they fill the range 0-255 (8 bit binary)
maximum = max(List)
minimum = min(List)
for i in range (0,size):
    List[i] =  255*(float(List[i])-float(minimum))/(float(maximum) - float(minimum))
    List[i] = int(List[i])


#######################REMAPPINGCODE FOR NANOWIRE1 SIDE VIEW ################################
#number1 = int(0.7*size)
#number2 = int(0.3*size)
#
#for i in range (number1,size):
#   List[i-number2] = List[i]

#for i in range (0,number2):
 #   List[i+ number1] = List[i] 
  #  List[i] = 0 


#########################REMAPPINGCODE FOR NANOWIRE2 SIDE VIEW #####################################
#New = [None] * size
#for i in range (0,(width/2)):
 #   for j in range(0,(height/2)):
  #    New[(size/2+width/2)+i+(j*width)] = List[i+j*width]
   #   New[i+(j*width)] = List[(size/2+width/2)+i+(j*width)]

#for i in range ((width/2),(width)):
 #   for j in range(0,(height/2)):
  #    New[(size/2)+(i-(width/2))+(j*width)] = List[i+j*width]
   #   New[i+(j*width)] = List[(size/2)+(i-(width/2))+(j*width)]

#for i in range (0,size):
 #   List[i] = New[i]




#Create greyscale image with the given intensity values.
newim = Image.new('L', (height,width))
pixels1 = newim.load()
for i in range (0,(height-1)):
    for j in range(0,(width-1)):
        l = List[i+height*j] 
 
        pixels1[i,j] = (l) 

newim.save("newim", "JPEG")