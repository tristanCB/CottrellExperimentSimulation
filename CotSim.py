#!/usr/bin/env python
import matplotlib.pyplot as plt
import math
import numpy
import argparse

########### INPUT VARIABLES USED IN EXECEL SUPPLIED ###########
# D = 1.00E-06 	# Diffusion Coefficient (cm^2/s)
# A = 1			# Area of the Electrode (cm^2)
# F = 96485 	# Faraday's Constant (C/mol)
# c = 1.00E-06 	# Concentration (mol/cm3)
# n = 1			# n electron process
# tk = 1.08E+06	# Know carachteristic time of the simulation (s)

########### DERIVED VARIABLES USED IN EXCEL EXAMPLE ###########
# NOT USED IN PROGRAM BUT GIVEN ON EXCEL SHEET 
# delT = tk / l			
# delx = (D*tk/Dm/l)**0.5	
# jmax = 4.2*(l**0.5)
###############################################################

def cottrell_exp(l,Dm,Boxes,tOVERtk):
	""" 
	Electrochem Simulation. Problem B.2 Bard. Mostly a translate of Fortran program Figure B.2.1 (Bard ) 
	This function takes 4 input arguments (l, Dm, Boxes, tOVERtk).
	A simulation of the Cottrell experiment which compares simulated non dimentional current to the current as derived from the cottrell equation.
	Then calculates dimentioneless concentration FA and FB over specified amount of boxes (starting from the electrode surface extending into the solution) at a specified dimentionless time.

	########### Example variables ###########
	# l = 50			# number of iterations corresponding to the known carachteristic time pf the simulation
	# Dm = 0.40			# Model diffusion coefficient in a simulation (None)
	# Boxes = 12		# How many boxes are to be plotted 
	# tOVERtk = 0.2		# Non dimentional Time
	
	Example use:
	>>> cottrell_exp(50,0.4,12,0.2)

	"""

	# At what non dimentinal time if FA and FB to be plotted 
	NonDTime = int(l*tOVERtk)	

	# Initializing some variable
	Z = list()
	T = list()
	ZCOTT = list()
	R = list()
	FA = list()
	FB = list()
	X = list()
	ConcFAfromEq5213 = list()
	ConcFBfromEq5213 = list()

	########## INITIAL CONDTS ###########
	# Electro reactant A is uniformely distributed initially.
	# Potential step is is applied at t = 0 forcing surface conc. of A to zero by converting it faradaically to species B
	# Start by setting up arrays to represent fractional conc. of A and B in each box
	# arrays initiallized to reflect uniform conc. of A and abscence of B
	# Old and new arrays for each species.... 

	# Setting up lists to represent the fractional concentrations of A and B in each box. 
	FAOLD = [1] * l
	FANEW = [1] * l
	FBOLD = [0] * l
	FBNEW = [0] * l

	# New conc. arrays calculated from old ones according to laws of diffusion
	k = 0
	while k != l:
		# print("Iteration: ", k)
		# Diffusion beyond the first box
		# for any experiment that has proceeded for time t will alter the solution from it's bulk carachter for a distance no larger than about 6*(Dt**1/2)
		jmax = int(4.2*(k**0.5))
		for j in range(1,jmax):
			FANEW[j] = FAOLD[j]+Dm*(FAOLD[j-1]-2*FAOLD[j]+FAOLD[j+1])
			FBNEW[j] = FBOLD[j]+Dm*(FBOLD[j-1]-2*FBOLD[j]+FBOLD[j+1])

		# Diffusion into the first box
		FANEW[0] = FAOLD[0]+ Dm*(FAOLD[1]-FAOLD[0])
		FBNEW[0] = FBOLD[0]+ Dm*(FBOLD[1]-FBOLD[0])

		
		#Faradaic conversion and current flow

		# Current Z(k) is calculated from amount of A converted
		Z.append(((l/Dm)**0.5)*FANEW[0])
		# FBNEW(0) is incremented by an equal amount to reflect faradaic conversion.
		FBNEW[0] = FBNEW[0] + FANEW[0]
		# Boundary cond. require that species A be zero in first box
		FANEW[0] = 0

		# if k == l/2: #this was what was given in the bard book also l/2 was 50
		# Type out conc. arrays when t/tk = 0.2 --> this will occur at iteration l/5
		# We are chosing a carachteristic time where the concentration profiles are to be plotted
		if k == NonDTime:
			#for j in range(0,jmax):
			for j in range(0,Boxes): #In the assignment we only want the first 12 boxes
				# X.append((j-1)/((Dm*l)**0.5))
				# not sure about this one here but the following line of code gives conc. profile that starts at zero 
				X.append((j)/((Dm*l)**0.5))
				FA.append(FANEW[j])
				FB.append(FBNEW[j])

		#Set up Old arrays for next iteration
		for j in range(0,jmax):
			FAOLD[j] = FANEW[j]
			FBOLD[j] = FBNEW[j]

		k = k+1 # Add an iteration of time. This terminates the while loop when k reaches l iterations.

	# Current time curve
	for k in range(1,l):
		T.append((k-0.5)/l)
		ZCOTT = ((float(math.pi)*T[k-1])**0.5)**-1
		R.append(Z[k]/ZCOTT)

	# From trial and error we find that by replacing x by nondimentinal distance and Do*t by t/tk in equation 5.2.13 we obtain 	
	distance = numpy.linspace(0, 2, num=50)
	for x in distance:
		conc = math.erf(x/(2*(tOVERtk)**0.5))
		ConcFAfromEq5213.append(conc)
		ConcFBfromEq5213.append(1-conc)

	###################################################
	################ Plotting results #################
	###################################################

	##### First Chart #####
	plt.figure(1)
	# Set chart title.
	# plt.title("Simulation for Question B.2 from (Bard, 2001)")

	# Set x, y label text.
	plt.xlabel("t/tk")
	plt.ylabel("Z/Zcott")
	plt.plot(T,R)
	# Draw point based on above x, y axis values.
	plt.scatter(T, R, s=10)

	plt.axis([0, 1.0, 0, 2.0])

	##### Second Chart #####
	plt.figure(2)
	# Set chart title.
	# plt.title("Concentration profiles of the Simulation at t/tk = 0.2 for Question B.2 from (Bard, 2001)")

	# Set x, y label text.
	plt.xlabel("X")
	plt.ylabel("f")
	plt.axis([0, 3, 0, 1])

	# Draw point based on above x, y axis values.
	FAsim = plt.scatter(X, FA, s=10, edgecolors='none', c='blue', label="FA sim")

	# Draw point based on above x, y axis values.
	FBsim = plt.scatter(X, FB, s=10, edgecolors='none', c='orange', label="FB sim")

	# Draw a line representing the concentration as derived from equation 5.2.13 in Bard for electroreactants FA and FB
	FAeq, = plt.plot(distance,ConcFAfromEq5213, label="FA as derived from eq 5.2.13")
	FBeq, = plt.plot(distance,ConcFBfromEq5213, label="FB as derived from eq 5.2.13")

	legend = plt.legend(handles=[FAeq,FBeq,FAsim,FBsim], loc=1)

	plt.show()

if __name__ == "__main__":
	# In this case the file is being run directly
	parser = argparse.ArgumentParser()
	parser.add_argument('--l', type=int, help='number of iterations corresponding to the known carachteristic time pf the simulation', default=50)
	parser.add_argument('--Dm', type=float, help='Model diffusion coefficient in a simulation', default=0.4)
	parser.add_argument('--Boxes', type=int, help='How many boxes are to be plotted', default=12)
	parser.add_argument('--tOVERtk', type=float, help='Non dimentional Time', default=0.2)

	args = parser.parse_args()

	# Run the simulation and graph results
	cottrell_exp(args.l,args.Dm,args.Boxes,args.tOVERtk)