'''
--------------------------------------------------
--- This scripts reads a frcmod containing only DIHEDRAL data (from paramfit)
--- This script is not robust: the input file has to be formatted as paramfit output , columswise and headerwise
--- The order for the dihedral is set up as 3, if needed to be changed change the o value after the imports. all dihedral needs to be of the same order
--------------------------------------------------  

'''
import matplotlib.pyplot as plt
import math
import sys

o=3

if len(sys.argv)  != 1  : sys.exit('need  one argument : frcmod produced by paramfit' )

def term_fourrier_serie(barriere, phase  , angle , n) :
	return barriere/2 * (1+ math.cos(n*math.radians(angle)-math.radians(phase) ))
 

input=open(sys.argv[1], 'r') 
inputcontent=input.readlines()
x=[]
y=[]


for l in range(len(inputcontent[4::o])) :  # datas begins at line 4 for paramfit output 
	del x[:]
	del y[:]
	for angle in range(-180,181):
		energy=0
		for n in range(o):    # Let's add all the terms of the fourrier serie
		  line=inputcontent[4+l*o+n]
                  # DEBUG if angle== -50 and n ==0:   print float(line[50:65])   
		  if len(line) == 1 : sys.exit('reach end of the file' )
		  energy+=term_fourrier_serie(float(line[16:30]),  float(line[32:45]), angle , float(line[50:65]))  # dirty parsing 
		x.append(angle)
                y.append(energy)
		
   	fig, ax = plt.subplots( nrows=1, ncols=1 )  # create figure 
	ax.plot(x, y )
	plt.title(inputcontent[4+l*o][0:12].replace(' ','') )
	fig.savefig('%s.png' %(inputcontent[4+l*o][0:12].replace(' ','') ) )   # save the figure to file
	print 'new file created : %s.png' %(inputcontent[4+l*o][0:12].replace(' ','') ) 
	plt.close(fig)

