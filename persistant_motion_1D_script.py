from persistant_motion_with_PE_landscape import main

def script(): 
	for alpha in [1.0]:
		for beta in [1.0]:
			for which_force in ["infinite"]:
				main(alpha, beta, which_force)

script()