from persistant_motion_2D_PE_landscape import main

def script(): 
	for alpha in [0.5,5.0]:
		for beta in [0.5,5.0]:
				main(alpha, beta)

script()