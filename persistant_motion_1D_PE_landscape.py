import numpy as np
import matplotlib.pyplot as plt
import time
from sys import argv 
##=============================================================================


def main(alpha, beta, which_force):
	##initial stuff (very specific, I know, but those four variables... that are like, muy importante and related to alpha and beta)
	##alpha = (tau_av)*(prop_force)/L
	##beta = k*L/(prop_force)
	k = 1.0
	prop_force_mag = 1.0
	L = beta*prop_force_mag/k
	tau_av = alpha*L/prop_force_mag
	##time variables dt = time step, tmax = max, time 
	tmax = 500.0
	dt = 0.01
	##time switch variables 
	number_of_switches = 2.0*tmax/tau_av
	
	## array of times
	tx = np.arange(0.00,tmax,dt)
	##setting switching time 
	tau = np.random.exponential(tau_av, number_of_switches).cumsum()
	##initial position
	xo = 0.00
	##initial propulsive force
	prop_force_i = prop_force_mag*np.random.choice([-1,1])
	##force 
	
	###Functions 
	prop_array = random_prop_force(tau,tx,prop_force_i)
	position_array = position(prop_array, xo, dt, L, k, which_force)
	
	#PE landscape
	# xarr = np.linspace(-L-1,L+1,200)
	# farr = np.array([force(xi, L, k, which_force) for xi in xarr])
	# plt.plot(xarr, farr)
	# plt.show()
	
	#pos vs time plot
	# plt.plot(tx, position_array)
	# plt.xlabel("Time (s)")
	# plt.ylabel("Position")
	# plt.title("Position vs Time")
	# plt.show()

	##histogram 1: rice rule 2n**1/3 - intermediate between previous two
	nbins_1 = int((tx.size)**(1.0/3.0))

	##other graph stuff?
	xmax = xaxis(k, L, prop_force_i, which_force)
	domain = np.ndarray.round(np.arange(-xmax, xmax+0.1, 0.1), 1)
	U_landscape = np.array([potential_energy_landscape(x, L, which_force) for x in domain])
	prob_eq = np.exp(-U_landscape) 
	prob_eq /= np.trapz(prob_eq, domain)
	
	##histogram plot
	plt.hist(position_array, bins=nbins_1, normed=True)
	H, bins = plt.hist(position_array, bins=nbins_1, normed=True)[:2]
	plt.plot(domain, U_landscape, label = 'Potential Landscape')
	plt.plot(domain, prob_eq, label = 'Gaussian Probability Density')
	plt.ylim([0,round(H.max(),1)+0.5])
	plt.title("Histogram of Particle Position (%s, $\\alpha = %s$, $\\beta = %s$)" %(which_force, alpha, beta))
	plt.xlabel("Position")  
	plt.legend(loc = 'upper right')
	hist_name = "PDF_1D_%s_alpha%s_beta%s_%s.png" %(which_force,tau_av, L, which_force)
	histfilename = "HISDAT_1D_lin_alpha%s_beta%s_%s" %(alpha, beta, which_force)
	np.savez(histfilename, hist=H,bins=bins,domain=domain, U_landscape=U_landscape, prob_eq=prob_eq)
	plt.savefig(hist_name)
	plt.close()

	##other 
	return
	
##=============================================================================
	
##returns an array of times and an array of velocities at those times 
def random_prop_force(tau,tx,p):
	i = 0 
	prop = np.zeros(tx.size)
	for j in range(tx.size): 
		if i<tau.size:
			if tx[j] > tau[i]: 
				p = np.random.choice([-p, p])
				i = i+1 
		prop[j] = p
	return prop 
		
##returns an array of positions based on the previously calculated velocities 
def position(prop, x0, dt, L, k, which_force):
	x = np.zeros(prop.size)
	x[0] = x0
	for i in range(prop.size-1):
		x[i+1] = x[i] + dt *(prop[i]+force(x[i],L,k,which_force))
	return x

##=============================================================================

def force(x,L,k,which_force):
	if which_force == 'linear':
		if x <= -L: 
			f = -k*(x+L)
		elif x >= L:
			f = -k*(x-L)
		else:
			f = 0.0
	if which_force == 'infinite':
		if x <= -L and x>-L-1: 
			f = -(x+L)/(1-(x+L)**2)
		elif x >= L and x<L+1:
			f = -(x-L)/(1-(x-L)**2)
		else:
			f = 0.0
	return f
##=============================================================================

def xaxis(k, L, vo, which_force):
	v = abs(vo)
	if which_force == 'linear':
		xmax = L + 4.0
	elif which_force == 'infinite':
		xmax = L + 0.9
	return xmax
	
##=============================================================================

def potential_energy_landscape(x, L, which_force):
	if which_force == 'linear': 	
		if x < -L: 
			U = ((x+L)**2)/2
		elif x >= -L and x <= L: 
			U = 0.0
		elif x >= L:
			U = ((x-L)**2)/2
	elif which_force == 'infinite': 
		G = 1-(x-L)**2
		if x <= -L and x > -L-1:
			U = (-1/2)*np.log(1-(x+L)**2)
		elif x > -L and x < L:
			U = 0.0
		elif x >= L and x < L+1:
			U = (-1/2)*np.log(1-(x-L)**2)
	return U

##=============================================================================	

	
if __name__=="__main__":
	main(float(argv[1]), float(argv[2]), argv[3])

