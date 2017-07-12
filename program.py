#This program uses the output file from "multislice.c" which is a text file named output.txt holding intensity values of each x and y position.
#From this output file an image is created.

import subprocess
import os



var = os.system("./multislice")
print(var)

# Image creation from TEM data (Converts Intensity values to integers between 0 and 255 (8 bits) to create an image)
from PIL import Image

parameters = open("size.txt").readlines()
height = int(parameters[0])
width = int(parameters[1])

size = height*int(width)

List = open("output.txt").readlines()

maximum = max(List)
minimum = min(List)


for i in range (0,size):
    List[i] =  255*(float(List[i])-float(minimum))/(float(maximum) - float(minimum))
    List[i] = int(List[i])




newim = Image.new('L', (height,width))
pixels1 = newim.load()


for i in range (0,(height-1)):
    for j in range(0,(width-1)):
        l = List[i+height*j]
 
        pixels1[i,j] = (l)

newim.save("newim", "JPEG")
