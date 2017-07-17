#This program takes the output file from "multislice.c" containing intensity values and turns them into a greyscale image. 

import subprocess #Subprocess used to run multislice program from python script
import os 
from PIL import Image #Python Image Libary used to create greyscale image


# Running the multislice program
var = os.system("./multislice")
print(var)


# Reading in parameters from c program output
parameters = open("size.txt").readlines()
height = int(parameters[0])
width = int(parameters[1])
size = height*int(width)

# Open output file from multislice and read in the intensity values
List = open("output.txt").readlines()

# Find the maximum and minimum intensity value
maximum = max(List)
minimum = min(List)

# Normalise the data so that the intensity values fall between 0 and 255 (8 bits).
for i in range (0,size):
    List[i] =  255*(float(List[i])-float(minimum))/(float(maximum) - float(minimum))
    List[i] = int(List[i])
    print(List[i])



newim = Image.new('L', (height,width))
pixels1 = newim.load()

# Set pixels equal to their corresponding intensity values. 
for i in range (0,(height-1)):
    for j in range(0,(width-1)):
        l = List[i+height*j]
 
        pixels1[i,j] = (l)

newim.save("newim", "JPEG")