import numpy as np
import matplotlib.pyplot as plt
from sys import argv

## ============================================================================	

def plot2d(histfile):
	##loading data
	dat = np.load(histfile)
	H_rad = dat["H_rad"]
	bins_rad = dat["bin_center_rad"]
	domain_rad = dat["domain_rad"]
	U_landscape = dat["U_landscape"]
	prob_eq = dat["prob_eq"]
	alpha = dat["alpha"]
	beta = dat["beta"]
	tmax = dat["tmax"]
	
	## Plot
	plt.plot(bins_rad, H_rad, "b-")
	plt.fill_between(bins_rad, 0, H_rad, color="b", alpha=0.1)
	plt.plot(domain_rad, U_landscape, label = '$U(r)$', color = 'r')
	plt.plot(domain_rad, prob_eq, label = '$\\rho_{\\rm eq}(r)$', color = 'g')
	plt.ylim([0,round(H_rad.max(),1)+0.1])
	plt.title("Radial Density ($\\alpha = %s$, $\\beta = %s$, tmax = %s)" %(alpha, beta,tmax))
	plt.xlabel("Position")  
	plt.legend(loc = 'upper right')
	plt.savefig("Datafiles/PDFr_2D_tmax%s_alpha%s_beta%s.jpg" %(tmax,alpha,beta))
	
	return
	
## ============================================================================	

if __name__=="__main__":
	plot2d(argv[1])
