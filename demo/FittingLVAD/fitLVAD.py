import csv
import os
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
from scipy.interpolate import LinearNDInterpolator

def PQr(x, a, b, c, d):
    return a + b/sqrt(c**2 - d*x**2)#b*x + c*x**2 + d*x**3

def readdata(pumpspeed):

	filename = open("HeartMatePumpCharacteristics.csv", "r")
	file = csv.DictReader(filename)
	
	P = []
	Q = []
	cnt = 0
	for col in file:
		

		Q_ = col[str(pumpspeed)+'rpm_x']
		P_ = col[str(pumpspeed)+'rpm_y']

		try:
			P_ = float(P_)
			Q_ = float(Q_)

			P.append(P_)
			Q.append(Q_)

		except ValueError:
			continue;



	return P, Q


pumpsp_arr = [3000, 4000, 5000, 6000, 7000, 8000, 9000]

points = []
vals = []
maxQ = []
maxP = []
for sp in pumpsp_arr:
	print "pump speed = ", sp
	P, Q = readdata(sp)

	[points.append((p,sp)) for p in P]
	[vals.append(Q[p]) for p in np.arange(0, len(Q))]

	plt.plot(Q, P, 'o', label=str(sp))
	Pfit = np.linspace(0, max(P))
	# Q as a function of P
	Qfunc = interp1d(P, Q, kind='linear', fill_value="extrapolate")
	Qfit = Qfunc(Pfit)

	#P as a function of Q
	Pfunc = interp1d(Q, P, kind='linear', fill_value="extrapolate")

	maxQ.append(Qfunc(0))
	maxP.append(Pfunc(0))

	# Add end points
	points.append((0, sp))
	vals.append(Qfunc(0))

	points.append((Pfunc(0), sp))
	vals.append(0)
	
	
maxPfunc = interp1d(pumpsp_arr, maxP, kind='linear', fill_value="extrapolate")

# Define 2D interplation
np.savez("HeartMate.npz",
	 points = points,
	 vals = vals)

def QLVAD(pressure, speed):
	points = np.load("HeartMate.npz")["points"]
	vals =  np.load("HeartMate.npz")["vals"]
	interp = LinearNDInterpolator(points, vals)
	
	return interp(pressure,speed)


pumpsp_arr2 = np.linspace(3000,9000,13)
for sp in pumpsp_arr2:
	print "pump speed = ", sp

	Pfit = np.linspace(0, maxPfunc(sp))
	Qfit = [QLVAD(p, sp) for p in Pfit]
	plt.plot(Qfit, Pfit, '-', label=str(sp))

	

plt.figure(1)
plt.legend()
plt.xlabel("Q (L/min)")
plt.ylabel("dP (mmHg)")
plt.savefig("LVADPumpCharacteristics.png")


