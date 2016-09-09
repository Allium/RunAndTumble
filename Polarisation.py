import numpy as np
import scipy
from scipy.signal import resample
import matplotlib.pyplot as plt
import os, glob, optparse
from time import time

from persistant_motion_2D_PE_landscape import radial_PE_landscape

##=============================================================================

def main():
	"""
	NAME
		Polarisation.py
	
	EXECUTION
		python Polarisation.py [path] [flags]
	
	ARGUMENTS
		histfile	path to density histogram
		dirpath 	path to directory containing histfiles
		
	FLAGS
		-v	--verbose	False
		-s	--show		False
			--nosave	False
		-a	--plotall	False
	"""
	
	me = "Polarisation.main: "
	t0 = time()
	
	## Options
	parser = optparse.OptionParser(conflict_handler="resolve")
	parser.add_option("-s","--show",
		dest="showfig", default=False, action="store_true")
	parser.add_option("-v","--verbose",
		dest="verbose", default=False, action="store_true")
	parser.add_option("--nosave",
		dest="nosave", default=False, action="store_true")
	parser.add_option("-a","--plotall",
		dest="plotall", default=False, action="store_true")
	parser.add_option("-h","--help",
		dest="help", default=False, action="store_true")		
	opt, args = parser.parse_args()
	if opt.help: print main.__doc__; return
	path	= args[0]
	showfig = opt.showfig
	verbose = opt.verbose
	nosave	= opt.nosave
	plotall = opt.plotall
	
	if plotall and os.path.isdir(path):
		showfig = False
		allfiles(path, nosave, verbose)
		
	elif os.path.isfile(path):
		pressure_pdf_file(path, nosave, verbose)
	elif os.path.isdir(path):
		pressure_dir(path, nosave, verbose)
	else:
		raise IOError, me+"You gave me rubbish. Abort."
	
	if verbose: print me+"execution time",round(time()-t0,2),"seconds"
	if showfig: plt.show()
	
	return
	
##=============================================================================
def polarisation_pdf_file(histfile, nosave, verbose):
	"""
	Make radial PDF and P(r) plot for a single file
	"""
	me = "Polarisation.polarisation_pdf_file: "
	t0 = time()

	## Filename
	plotfile = os.path.dirname(histfile)+"/POLARr"+os.path.basename(histfile)[6:-4]+".jpg"
	
	## Get pars from filename and convert to R (box size) and tau (memory time)
	a, b = filename_par(histfile, "_alpha"), filename_par(histfile, "_beta")
	R = b
	tau = a*b
	
	## Space (for axes)
	data = np.load(histfile)
	r = data["bin_center_rad"]
	th = data["bin_center_theta"]
	
	## Load histogram, normalise
	H = data["H_rad"]
	rho = H / np.trapz(np.trapz(H*2*np.pi*r,r,axis=0),th)

	## Potential and equilibrium result
	r_eq = np.linspace(r[0],R+4.0,50*int(R+4.0-r[0]))
	th_eq = th
	U = np.array([radial_PE_landscape(ri, R) for ri in r_eq])
	rho_eq = np.exp(-U) / np.trapz(np.exp(-U)*2*np.pi*r_eq, r_eq)
			
	##---------------------------------------------------------------			
	## PLOT SET-UP
	
	figtit = "Density and polarisation; quadratic potential; $\\alpha="+str(a)+"$, $\\beta = "+str(b)+"$"
	fsa, fsl, fst = 16,12,16
	
	fig, axs = plt.subplots(3,1,sharex=True)
		
	##---------------------------------------------------------------	
	## PDF PLOT

	ax = axs[0]
	
	ax.plot(r,rho.sum(axis=0),   "b-", label="RTP")
	ax.fill_between(r, 0, rho, color="b", alpha=0.1)
	ax.plot(r_eq,rho_eq,"r-", label="Eq.")
	ax.fill_between(r_eq, 0, rho_eq, color="r", alpha=0.1)
	
	## Accoutrements
	ax.set_xlim([0.0,R+4.0])
	ax.set_ylim(bottom=0.0, top=ax.get_ylim()[1])
	ax.set_ylabel("$\\rho(r,\\phi)$", fontsize=fsa)
	ax.grid()
	
	## Wall
	ax.plot(r_eq, U*ax.get_ylim()[1]/U[-1], "k--", label="$U(r)$")
	
	##---------------------------------------------------------------
	## POLARISATION
	
	## CALCULATIONS
	m1c, m1s = calc_angular_moment(th,rho,1)
	m2c, m2s = calc_angular_moment(th,rho,2)
	
	m1c_eq, m1s_eq = calc_angular_moment(th_eq,rho_eq,1)
	m2c_eq, m2s_eq = calc_angular_moment(th_eq,rho_eq,2)
	
	## FIRST MOMENT
	ax = axs[1]
	
	## Pressure
	ax.plot(r,m1c, "b-",  label="RTP $\\langle\\cos\\theta\\rangle(r)")
	ax.plot(r,m1s, "b--", label="RTP $\\langle\\sin\\theta\\rangle(r)")
	ax.plot(r_eq,m1c_eq,"r-",  label="Eq. $\\langle\\cos\\theta\\rangle(r)")
	ax.plot(r_eq,m1s_eq,"r--", label="Eq. $\\langle\\sin\\theta\\rangle(r)")
	
	## Accoutrements
	ax.set_ylabel("First moments", fontsize=fsa)
	ax.grid()
	
	## Wall
	ax.plot(r_eq, U*ax.get_ylim()[1]/U[-1], "k--", label="$U(r)$")

	## SECOND MOMENT
	ax = axs[2]
	
	## Pressure
	ax.plot(r,m2c, "b-",  label="RTP $\\langle\\cos^2\\theta\\rangle(r)")
	ax.plot(r,m2s, "b--", label="RTP $\\langle\\sin^2\\theta\\rangle(r)")
	ax.plot(r_eq,m2c_eq,"r-",  label="Eq. $\\langle\\cos^2\\theta\\rangle(r)")
	ax.plot(r_eq,m2s_eq,"r--", label="Eq. $\\langle\\sin^2\\theta\\rangle(r)")
	
	## Accoutrements
	ax.set_ylabel("Second moments", fontsize=fsa)
	ax.grid()
	
	## Wall
	ax.plot(r_eq, U*ax.get_ylim()[1]/U[-1], "k--", label="$U(r)$")
	
	## Bottom figure
	ax.set_xlabel("$r$", fontsize=fsa)
	ax.legend(loc="best",fontsize=fsl)
	
	##---------------------------------------------------------------
	
	## Tidy figure
	fig.suptitle(figtit,fontsize=fst)
	fig.tight_layout();	plt.subplots_adjust(top=0.9)	
	
	if not nosave:
		fig.savefig(plotfile)
		if verbose: print me+"plot saved to",plotfile
	
	return

##=============================================================================
def allfiles(dirpath, nosave, verbose):
	for filepath in np.sort(glob.glob(dirpath+"/*.npz")):
		pressure_pdf_file(filepath, nosave, verbose)
		plt.close()
	return

## ============================================================================

def filename_par(filename, searchstr):
	"""
	Scrape filename for parameters and return a dict.
	"""
	start = filename.find(searchstr) + len(searchstr)
	finish = start + 1
	while unicode(filename[start:].replace(".",""))[:finish-start].isnumeric():
		finish += 1
	return float(filename[start:finish])
	
##=============================================================================

def calc_angular_moment(th,rho,k):
	"""
	Calculate polarisation given density function of coordinate r and body angle th.
	"""
	mkc = 2*np.pi* np.trapz(rho*np.power(np.cos(th),k),th,axis=1) / np.trapz(rho,th,axis=1)
	mks = 2*np.pi* np.trapz(rho*np.power(np.sin(th),k),th,axis=1) / np.trapz(rho,th,axis=1)
	return (mkc, mks)
	
## ============================================================================
## ============================================================================
if __name__=="__main__":
	main()