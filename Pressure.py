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
		Pressure.py
	
	EXECUTION
		python Pressure.py [path] [flags]
	
	ARGUMENTS
		histfile	path to density histogram
		dirpath 	path to directory containing histfiles
		
	FLAGS
		-v	--verbose	False
		-s	--show		False
			--nosave	False
		-a	--plotall	False
	"""
	
	me = "Pressure.main: "
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
def pressure_pdf_file(histfile, nosave, verbose):
	"""
	Make radial PDF and P(r) plot for a single file
	"""
	me = "Pressure.pressure_pdf_file: "
	t0 = time()

	## Filename
	plotfile = os.path.dirname(histfile)+"/PRESSr"+os.path.basename(histfile)[6:-4]+".jpg"
	
	## Get pars from filename and convert to R (box size) and tau (memory time)
	a, b = filename_par(histfile, "_alpha"), filename_par(histfile, "_beta")
	R = b
	tau = a*b
	
	## Space (for axes)
	data = np.load(histfile)
	r = data["bin_center_rad"]
	
	## Load histogram, convert to normalised pdf
	H = data["H_rad"]
	try:	H = H.sum(axis=1)
	except ValueError:	pass
	rho = H / np.trapz(H*2*np.pi*r, r)

	## Potential and equilibrium result
	r_eq = np.linspace(r[0],R+4.0,50*int(R+4.0-r[0]))
	U = np.array([radial_PE_landscape(ri, R) for ri in r_eq])
	rho_eq = np.exp(-U) / np.trapz(np.exp(-U)*2*np.pi*r_eq, r_eq)
			
	##---------------------------------------------------------------			
	## PLOT SET-UP
	
	figtit = "Density and pressure; quadratic potential; $\\alpha="+str(a)+"$, $\\beta = "+str(b)+"$"
	fsa, fsl, fst = 16,12,16
	
	fig, axs = plt.subplots(2,1,sharex=True)
		
	##---------------------------------------------------------------	
	## PDF PLOT

	ax = axs[0]
	
	ax.plot(r,rho,   "b-", label="RTP")
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
	## PRESSURE
	
	## CALCULATIONS
	p		= calc_pressure(r, rho, r_eq, U, spatial=True)
	p_eq	= calc_pressure(r_eq, rho_eq, r_eq, U, spatial=True)
	
	## PLOT
	ax = axs[1]
	
	## Pressure
	ax.plot(r,p,   "b-", label="RTP")
	ax.plot(r_eq,p_eq,"r-", label="Eq.")
	
	## Accoutrements
	ax.set_ylim(bottom=0.0, top=ax.get_ylim()[1])
	ax.set_xlabel("$r$", fontsize=fsa)
	ax.set_ylabel("$P(r)$", fontsize=fsa)
	ax.legend(loc="best",fontsize=fsl)
	ax.grid()
	
	## Wall
	ax.plot(r_eq, U*ax.get_ylim()[1]/U[-1], "k--", label="$U(r)$")

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

##=============================================================================
def pressure_dir(dirpath, nosave, verbose):
	"""
	Plot pressure at "infinity" against alpha for all files in directory.
	"""
	me = "Pressure.pressure_dir: "
	t0 = time()
		
	## File discovery
	histfiles = np.sort(glob.glob(dirpath+"/*.npz"))
	numfiles = histfiles.size
	if verbose: print me+"found",numfiles,"files"

	## Initialise
	A = np.zeros(numfiles) 
	B = np.zeros(numfiles)
	P = np.zeros(numfiles)
	P_eq = np.zeros(numfiles)
		
	## Loop over files
	for i,histfile in enumerate(histfiles):
	
		a, b = filename_par(histfile, "_alpha"), filename_par(histfile, "_beta")
		R = b
		tau = a*b
		
		## Space (for axes)
		data = np.load(histfile)
		r = data["bin_center_rad"]
		
		## Load histogram, convert to normalised pdf
		H = data["H_rad"]
		rho = H / np.trapz(H*2*np.pi*r, r)

		## Potential and equilibrium result
		r_eq = np.linspace(r[0],R+4.0,50*int(R+4.0-r[0]))
		U = np.array([radial_PE_landscape(ri, R) for ri in r_eq])
		rho_eq = np.exp(-U) / np.trapz(np.exp(-U)*2*np.pi*r_eq, r_eq)
		
		## Pressure array
		P[i] 	= calc_pressure(r, rho, r_eq, U)
		P_eq[i] = calc_pressure(r_eq, rho_eq, r_eq, U)
		A[i], B[i] = a, b
		
	## ------------------------------------------------	
	## Normalise pressure
	
	P /= P_eq
	
	fsa, fsl, fst = 16,12,16
	figtit = "Pressure ($P_{\\rm eq}(\\alpha,\\beta)=1.0$)"

	## ------------------------------------------------	
	## Phase diagram plot
	
	plotfile = dirpath+"/PRESSab_scatter.jpg"
	fig = plt.figure(); ax = plt.gca()
	
	## Plot data with variable colour and size
	plt.scatter(A, B, marker="o", c=np.log(P), s=500*P/P.max(), edgecolor="None", label="RT")
	cbar = plt.colorbar()
	
	ax.set_xscale("log");	ax.set_yscale("log")
	# ax.set_xlim(left=0.0);	ax.set_ylim(bottom=0.0)
	ax.set_xlabel("$\\alpha$",fontsize=fsa)
	ax.set_ylabel("$\\beta$",fontsize=fsa)
	cbar.set_label("$\\log P$", fontsize=fsa, rotation=270)
	fig.suptitle(figtit, fontsize=fst)
	ax.grid()
	
	if not nosave:
		fig.savefig(plotfile, dpi=2*fig.dpi)
		if verbose: print me+"Plot saved to",plotfile
	
	## ------------------------------------------------	
	## Lines plot
		
	AA = np.unique(A)
	BB = np.unique(B)
	
	## 2D pressure array: [A,B]
	PP = -np.ones([AA.size,BB.size])
	PP_eq = np.zeros(PP.shape)
	for i in range(AA.size):
		Aidx = (A==AA[i])
		for j in range(BB.size):
			Bidx = (B==BB[j])
			Pidx = Aidx*Bidx
			try:
				PP[i,j] = P[Pidx]
				PP_eq[i,j] = P_eq[Pidx]
			except ValueError:
				# pass
				PP[i,j] = np.nan
				PP_eq[i,j] = np.nan
	
	## Plotting
	
	fig = plt.figure(); ax = plt.gca()
	
	plotfile = dirpath+"/PRESSab_line.jpg"
	for i in range(BB.size):
		ax.plot(AA, PP[:,i], "o-", label="$\\beta="+str(BB[i])+"$")
	ax.set_xlabel("$\\alpha$",fontsize=fsa)
	# plotfile = dirpath+"/PRESSba_line.jpg"
	# for i in range(AA.size):
		# ax.plot(BB, PP[i,:], "o-", label="$\\alpha="+str(AA[i])+"$")
	# ax.set_xlabel("$\\beta$",fontsize=fsa)
	
	ax.set_xlim(left=0.0);	ax.set_ylim(bottom=0.0)
	# ax.set_xlim(right=2.0);	ax.set_ylim(top=2.0); plotfile = plotfile[:-4]+"_zoom.jpg"
	
	ax.set_ylabel("Pressure (normalised)",fontsize=fsa)
	fig.suptitle(figtit, fontsize=fst)
	ax.grid()
	ax.legend(loc="best",fontsize=fsl)

	if not nosave:
		fig.savefig(plotfile, dpi=2*fig.dpi)
		if verbose: print me+"Plot saved to",plotfile
	
	## ------------------------------------------------	

	if verbose: print me+"Plotting",round(time()-t0,2),"seconds."
	
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

def calc_pressure(r,rho,r_eq,U,spatial=False):
	"""
	Calculate pressure given density a a function of coordinate.
	Uses array resampling rather than interpolation multiplication.
	"""
	
	## Resample potential
	if U.size!=r.size:	U = scipy.interpolate.interp1d(r_eq,U)(r)
	
	## Force from potential
	force = -np.gradient(U, r[1]-r[0])
	
	## Pressure
	if spatial == True:
		P = -np.array([np.trapz(force[:i]*rho[:i], r[:i]) for i in range(1,r.size+1)])
	else:
		P = -np.trapz(force*rho, r)
	
	return P
	
## ============================================================================
## ============================================================================
if __name__=="__main__":
	main()