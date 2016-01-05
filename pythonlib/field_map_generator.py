#######################################################################################################################
#Created by Patrik Hallsjo @ University of Glasgow
#Need automatic dating through GIT, 
#Modified on 28/10-2015
#Created on 25/9-2015
#######################################################################################################################
#General python import
#######################################################################################################################
import os
import math
#######################################################################################################################
#Importing own python files
#######################################################################################################################

#######################################################################################################################
#Class generation
#######################################################################################################################
class field_map_generator:

	def __init__(self, Bmag, height, width, npanels):
		#Assume Bmag in tesla, height, width in meters and npanels as a number.

		#Setup all the required units. Only stepSize can/should be changed.
		self.mm = 1
		self.cm = 10 * self.mm
		self.meter = 1000*self.mm
		self.meter2 = self.meter*self.meter
		self.nanosecond = 1
		self.second = self.nanosecond * 1.e+9
		self.megaelectronvolt = 1
		self.eplus = 1
		self.megavolt = self.megaelectronvolt/self.eplus
		self.volt = 1.e-6 * self.megavolt
		self.tesla = self.volt * self.second / self.meter2

		self.Bx = Bmag
		self.hmax = height #* self.cm
		self.wmax = width #* self.cm
		self.seg = npanels
		self.stepSize =0.05

################################################
	def CreateRange(self,start,stop,step):
		r = start + step
		while r <= stop:
			yield r
			r += step

	def GetFieldValue(self,x,y,z):
		# Only define the magnetic field within the detector and ignore it
  		# outside of the detector. 
  		Bx = 0.
  		By = 0.
  		Bz = 0.

  		if( x <= self.wmax/2. and y <= self.hmax/2. ):
  			if(self.seg == 2):
  			#the field is divided into two equal sections
  				if(y > 0):
  					Bx = self.Bx
  				else:
  					Bx = -self.Bx
  			#There are two outside plates half the size of the inner plates
  			elif(self.seg > 2):
  				segSize = self.hmax/2./(self.seg - 1)
  				inSeg = math.ceil(abs(y/segSize))

  				if(inSeg % 2 == 0):
  					Bx = self.Bx
  				else:
  					Bx = -self.Bx

  		return Bx, By, Bz

	def Print_field_to_file(self,filename):
		'''
		Print data to file filename.
		'''
		outfile = open(filename,'w+')
		zCoord = 0
		yRange = self.CreateRange(-self.hmax/2,self.hmax/2,self.stepSize)

		for yCoord in yRange:
			xRange = self.CreateRange(-self.wmax/2,self.wmax/2,self.stepSize)
			for xCoord in xRange:
	    		        #xCoord tab yCoord tab zCoord tab Bx tab By tab Bz
	    		        Bx, By, Bz = self.GetFieldValue(xCoord,yCoord,zCoord)
				#print '%s\t%s\t%s\t%s\t%s\t%s' % (xCoord,yCoord,zCoord,Bx,By,Bz)
				data = '    %s    %s    %s    %s    %s    %s\n' % (round(xCoord,3),round(yCoord,3),round(zCoord,3),Bx,By,Bz)
				outfile.write(data)

	        outfile.close()

#######################################################################################################################
#File specific functions
#######################################################################################################################
#if __name__ == "__main__":

    #s=field_map_generator(1.5,2.0,3.5,2)
    #s.Print_field_to_file('test.txt')








