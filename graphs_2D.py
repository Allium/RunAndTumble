import numpy as np
import matplotlib.pyplot as plt
from sys import argv

## ============================================================================	

def main(): 
	
	dat = np.load(argv[1])
	plotrad2d(dat)
	plotxy2d(dat)
	plotangrad(dat)
	plotboxwall(dat)
	plotradinout(dat)
	plotang(dat)
	return
	
## ============================================================================	

def plotrad2d(dat):
	##loading data
	H_rad = dat["H_rad"]
	bins_rad = dat["bin_center_rad"]
	domain_rad = dat["domain_rad"]
	U_landscape = dat["U_landscape"]
	prob_eq = dat["prob_eq"]
	alpha = dat["alpha"]
	beta = dat["beta"]
	tmax = dat["tmax"]
	histnamerad = dat["histnamerad"]
	##Plot
	plt.plot(bins_rad, H_rad, "b-")
	plt.fill_between(bins_rad, 0, H_rad, color="b", alpha=0.1)
	plt.plot(domain_rad, U_landscape, label  = '$U(r)$', color = 'r')
	plt.plot(domain_rad, prob_eq, label = '$\\rho_{\\rm eq}(r)$', color = 'g')
	plt.ylim([0,round(H_rad.max(),1)+0.1])
	plt.title("Radial Density ($\\alpha = %s$, $\\beta = %s$, tmax = %s)" %(alpha, beta,tmax))
	plt.xlabel("Radius")  
	plt.legend(loc = 'upper right')
	plt.savefig("Results/"+str(histnamerad))
	plt.close()
	return

	
def plotxy2d(dat):
	##loading data
	alpha = dat["alpha"]
	beta = dat["beta"]
	tmax = dat["tmax"]
	H_2d = dat["H_2d"]
	bin_x = dat["bin_x"]
	bin_y = dat["bin_y"]
	stall_r = dat["stall_r"]
	extent_2d = dat["extent_2d"]
	histnamexy = dat["histnamexy"]
	r = dat["r"]
	##the plot thickens 
	plt.imshow(H_2d, extent = extent_2d, interpolation = "None", aspect = "auto")
	plt.colorbar()
	##circly things
	ang = np.linspace(0,2*np.pi,360)
	plt.plot(r*np.cos(ang),r*np.sin(ang),"g-",lw=3)
	plt.plot(stall_r*np.cos(ang),stall_r*np.sin(ang),"y-",lw=3)
	plt.xlim([-stall_r,stall_r])
	plt.ylim([-stall_r,stall_r])
	##other things that are also plotty 
	plt.title("Density ($\\alpha = %s$, $\\beta = %s$, tmax = %s)" %(alpha, beta,tmax))
	plt.xlabel("x-position") 
	plt.ylabel("y-position")
	plt.savefig("Results/"+str(histnamexy))
	plt.close()
	return 

def plotangrad(dat):
	##loading data
	alpha = dat["alpha"]
	beta = dat["beta"]
	tmax = dat["tmax"]
	H_ang = dat["H_ang"] 
	extent_ang = dat["extent_ang"]
	##plotting
	plt.imshow(H_ang, extent = extent_ang, interpolation = "None", aspect = "auto")
	plt.colorbar()
	plt.title("Angular Density at Radii ($\\alpha = %s$, $\\beta = %s$, tmax = %s)" %(alpha, beta,tmax))
	plt.xlabel("Radius")
	plt.ylabel("Angle")  
	plt.savefig("Results/PDFangrad_2D_alpha%s_beta%s_tmax%s.png" %(alpha,beta,tmax))
	plt.close()
	return
	
def plotang(dat):
	##loading data
	alpha = dat["alpha"]
	beta = dat["beta"]
	tmax = dat["tmax"] 
	bin_ang = dat["bin_ang"]
	H_ang = dat["H_ang"]
	##plotting
	bin_ang_centers = (bin_ang[:-1] + bin_ang[1:]) / 2.
	H_ang_sum = H_ang.sum(axis=0)
	plt.plot(bin_ang_centers, H_ang_sum, "b-")
	plt.fill_between(bin_ang_centers, 0, H_ang_sum, color = "b", alpha = 0.1)
	plt.title("Angular Density ($\\alpha = %s$, $\\beta = %s$, tmax = %s)" %(alpha, beta,tmax))
	plt.xlabel("Angle")
	plt.savefig("Results/PDFang_2D_alpha%s_beta%s_tmax%s.png" %(alpha,beta,tmax))
	plt.close()
	return

def plotboxwall(dat):
	##loading data
	alpha = dat["alpha"]
	beta = dat["beta"]
	tmax = dat["tmax"] 
	bin_ang = dat["bin_ang"]
	H_ang = dat["H_ang"]
	bin_rad_ang = dat["bin_rad_ang"]
	r = dat["radius"]
	##finding where the wall is... in terms of... I'm trying ok 
	for i in range(bin_rad_ang.size):
		if bin_rad_ang[i] > r: 
			wall = i
			r = r + 100
	bin_ang_centers = (bin_ang[:-1] + bin_ang[1:])/ 2.
	H_box = H_ang[:wall].sum(axis = 0)
	H_wall = H_ang[wall:].sum(axis = 0)
	plt.plot(bin_ang_centers, H_box, "b-", label = "Box")
	plt.plot(bin_ang_centers, H_wall, "k-", label = "Wall")
	plt.legend()
	plt.title("Angular Density at Box and Wall ($\\alpha = %s$, $\\beta = %s$, tmax = %s)" %(alpha, beta,tmax))
	plt.xlabel("Angle")
	plt.savefig("Results/PDFboxwall_2D_alpha%s_beta%s_tmax%s.png" %(alpha,beta,tmax))
	plt.close()
	return 
	
def plotradinout(dat):
	alpha = dat["alpha"]
	beta = dat["beta"]
	tmax = dat["tmax"] 
	bin_ang = dat["bin_ang"]
	H_ang = dat["H_ang"]
	bin_rad_ang = dat["bin_rad_ang"]
	r = dat["radius"]
	pi = np.pi
	##in or out 
	for i in range(bin_ang.size):
		if bin_ang[i] < (-pi/2):
			inwardneg = i
		elif bin_ang[i] > (pi/2):
			inwardpos = i
			pi = pi+10
	bin_rad_centers = (bin_rad_ang[:-1] + bin_rad_ang[1:])/ 2.
	H_in = (np.sum((H_ang[:bin_ang.size,:(inwardneg+1)].sum(axis=1), H_ang[:bin_ang.size,inwardpos:].sum(axis=1)), axis = 0))/(2*np.pi*bin_rad_centers)
	H_out = (H_ang[:bin_ang.size,(inwardneg+1):inwardpos].sum(axis=1))/(2*np.pi*bin_rad_centers)
	plt.plot(bin_rad_centers, H_in, "b-", label = "In")
	plt.plot(bin_rad_centers, H_out, "k-", label = "Out")
	plt.legend()
	plt.title("Radial Density at Radially Inwards and Outwards ($\\alpha = %s$, $\\beta = %s$, tmax = %s)" %(alpha, beta,tmax))
	plt.ylabel("Radius")
	plt.savefig("Results/PDFrinout_2D_alpha%s_beta%s_tmax%s.png" %(alpha,beta,tmax))
	plt.close()
	return 

	
	
			
		
	
## ============================================================================	

if __name__=="__main__":
	main()
