import numpy as np
import matplotlib.pyplot as plt

## ============================================================================	
## ============================================================================	

def main():
	"""
	Illustration of splitting the plot function from simulation / histogram function.
	"""
	
	simulate()
	plot()
	
	return

## ============================================================================	

def simulate():
	"""
	Makes some random data, histograms the data, and saves histogram to npz.
	"""

	## Random data (1D)
	npoints=5000
	v = np.random.randn(npoints)

	## Make and save histogram arrays
	H, bins = np.histogram(v, 2*npoints**0.3, normed=True)
	np.savez("data",H=H,bins=bins)
	print "Saved datafile to data.npz"
	
	return

## ============================================================================	
	
def plot():
	"""
	Reads in histogram file, plots, and saves figure.
	"""

	## Load data from file
	dat = np.load("data.npz")
	H = dat["H"]
	bins = dat["bins"]

	## Calculate bin centres
	x = 0.5*(bins[:-1]+bins[1:])

	## Plot
	plt.plot(x, H, "b-")
	plt.fill_between(x, 0, H, color="b", alpha=0.1)

	plt.savefig("fig.png")
	print "Saved plotfile to fig.png"
	plt.show()
	
	return

## ============================================================================	
## ============================================================================	
if __name__=="__main__":
	main()

