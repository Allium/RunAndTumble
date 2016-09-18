import numpy as np
import matplotlib.pyplot as plt
import time
from sys import argv 

##=============================================================================
##FIX HIST NAME RAD GOD DAMN IT
def main(alpha, beta):
	##initial stuff (very specific, I know, but those four variables... that are like, muy importante and related to alpha and beta)
	##alpha = (tau_av)*(prop_force)/radius
	##beta = k*radius/(prop_force)
	k = 1.0
	prop_force_i = 1.0
	r = beta*prop_force_i/k
	tau_av = alpha*r/prop_force_i
	##time variables 
	dt = 0.01
	tmax = 5000.0
	##time switch variables 
	number_of_switches = 2.0*tmax/tau_av
	##time switch variables -- don't overload memory
	while True:
		try:
			number_of_switches = 1.1*tmax/tau_av
			tau = np.random.exponential(tau_av, number_of_switches).cumsum()
			break
		except MemoryError:
			tmax /= 10.0
			print "Warning: tmax decimated.", tmax
	## array of times
	tx = np.arange(0.00,tmax,dt)
	##setting switching time 
	tau = np.random.exponential(tau_av, number_of_switches).cumsum()
	N = tau.size
	##initial position
	xo = 0.00
	yo = 0.00 
	##angle
	angle_initial = 2*np.random.rand()
	##force things
	stall_r = (prop_force_i/k)+r 

	##arrays, yo 
	xcomp, ycomp = random_angle(N)
	x_prop, y_prop = random_prop_force(tau,tx,prop_force_i,angle_initial,xcomp,ycomp)
	x_position, y_position, radius_array = position(x_prop, y_prop, xo, yo, dt, r, k)
	theta_array = radial_angle(x_position, y_position)
	##radial histogram things 
	domain_rad = np.ndarray.round(np.arange(0, stall_r+1.1, 0.1), 1)
	U_landscape = np.array([radial_PE_landscape(x, r) for x in domain_rad])
	prob_eq = (np.exp(-U_landscape))
	prob_eq /= np.trapz(prob_eq*2*np.pi*domain_rad, domain_rad)
 
	##bin sizes
	nbins = int(((tx.size*2+r*2+tmax*2)**1/2)**(1.0/3.0))
	Nrbins = 32
	Nthbins = 48
	##plot with angle shit I guess? idk man
	H_ang, bin_rad_ang, bin_ang = np.histogram2d(radius_array[:-1], theta_array, [Nrbins,Nthbins])[:3]
	extent_ang = [bin_rad_ang[0], bin_rad_ang[-1], bin_ang[0], bin_ang[-1]]
	
	
	##plot that was working two seconds ago that i broke even tho i didn't touch it
	H_2d, bin_x, bin_y = np.histogram2d(x_position, y_position, bins=nbins)
	extent_2d = [bin_x[0], bin_x[-1], bin_y[0], bin_y[-1]]
	histnamexy = "PDFxy_2D_alpha%s_beta%s_tmax%s.png" %(alpha,beta,tmax)


	"""histogram plot stuff"""
	H, bin_edges = np.histogram(radius_array, bins=nbins, normed=True)[:2]
	bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
	H_rad =  H/(2*np.pi*bin_centers)
	#width = bin_edges[1] - bin_edges[0]
	#plt.plot(bin_centers, H_rad, "b-")
	#plt.fill_between(bin_centers, 0, H_rad, color = "b", alpha = 0.1)
	
	"""other plot stuff"""
	#plt.plot(domain_rad, U_landscape, label = '$U(r)$', color = 'r')
	#plt.plot(domain_rad, prob_eq, label = '$\\rho_{\\rm eq}(r)$', color = 'g')
	#plt.ylim([0,round(H_rad.max(),1)+0.1])
	#plt.title("Radial Density ($\\alpha = %s$, $\\beta = %s$, tmax = %s)" %(alpha, beta,tmax))
	histnamerad = "PDFr_2D_alpha%s_beta%s_tmax%s.png" %(alpha,beta,tmax)
	

	'''what's going on'''
	
	""""saving data file"""
	histfilename = "HISDAT_2D_alpha%s_beta%s_tmax%s" %(alpha, beta, tmax)
	np.savez("Datafiles/"+histfilename, H_rad=H_rad,bin_center_rad=bin_centers,domain_rad=domain_rad, U_landscape=U_landscape, prob_eq=prob_eq,
	alpha=alpha, beta=beta, r=r, tmax=tmax, histnamerad=histnamerad, H_2d = H_2d, bin_x=bin_x, bin_y = bin_y, stall_r = stall_r, 
	extent_2d = extent_2d, histnamexy = histnamexy, H_ang = H_ang, bin_ang = bin_ang, bin_rad_ang = bin_rad_ang, extent_ang = extent_ang, radius = r)
	
	return

##=============================================================================
def random_angle(N):
	rangle = 2*np.pi*np.random.rand(N)
	return np.cos(rangle), np.sin(rangle)
##=============================================================================

def random_prop_force(tau,tx,prop_force,angle_initial,xcomp,ycomp):
	i = 0 
	x_prop = np.zeros(tx.size)
	y_prop = np.zeros(tx.size)
	for j in range(tx.size): 
		if i<tau.size:
			if tx[j] > tau[i]:
				x = xcomp[i]*prop_force
				y = ycomp[i]*prop_force
				i = i+1 
			elif j == 0:
				x = np.cos(angle_initial)*prop_force
				y = np.sin(angle_initial)*prop_force
		x_prop[j] = x
		y_prop[j] = y
	return x_prop, y_prop
	
##=============================================================================

def force(x, y, k, r):
	pos = (x**2+y**2)**(0.5)
	if pos >= r:
		xforce = -k*(x-r*(x/pos))
		yforce = -k*(y-r*(y/pos))
	else:
		xforce = 0.0
		yforce = 0.0
	return pos, xforce, yforce 
	
##=============================================================================

def position(propx, propy, xo, yo, dt, r, k):
	x = np.zeros(propx.size)
	y = np.zeros(propy.size)
	radius = np.zeros(propy.size)
	x[0] = xo
	y[0] = yo
	radius[0] = 0
	for i in range(propx.size-1):
		pos, xforce, yforce = force(x[i], y[i], k, r)
		x[i+1] = x[i] + dt *(propx[i]+xforce)
		y[i+1] = y[i] + dt *(propy[i]+yforce)
		radius[i+1] = (x[i+1]**2 + y[i+1]**2)**(0.5)
	return x, y, radius
	
##=============================================================================

def radial_PE_landscape(x, r):
	if x < r: 
		U = 0.0
	elif x >= r:
		U = ((x-r)**2)/2
	return U		
##=============================================================================

def radial_angle(x, y):
	delta_x = np.diff(x)
	delta_y = np.diff(y)
	ang1 = np.arctan2(y, x)
	theta = np.arctan2(delta_y, delta_x)-ang1[:-1]
	theta[theta<=-np.pi] += 2*np.pi
	theta[theta>np.pi] -= 2*np.pi
	return theta

##=============================================================================

if __name__=="__main__":
	main(float(argv[1]), float(argv[2]))



