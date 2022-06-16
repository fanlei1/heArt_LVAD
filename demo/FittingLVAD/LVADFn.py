from scipy.interpolate import LinearNDInterpolator
import numpy as np

class LVAD(object):

    def __init__(self, filename):

    	self.points = np.load(filename)["points"]
    	self.vals =  np.load(filename)["vals"]

    def Flowrate(self, pressure, speed):


    	interp = LinearNDInterpolator(self.points, self.vals)
	if(pressure < 0):
		flow = 0.0
	else:
		flow = interp(pressure,speed)
    	
    	return flow


