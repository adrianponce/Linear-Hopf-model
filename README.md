# Linear-Hopf-model
Estimation of network statistics of the Hopf whole-brain model using a linear noise approximation

The present codes estimate the network statistics of the Hopf whole-brain model (Deco et al. 2017). This model corresponds to a network of nonlinear oscillators described by a normal form of a supercritical Hopf bifurcation. In the case of weak noise and small non-linearities, one can estimate the network statistics using a linear approximation, without the need of simulations. In the linear approximation, the stationary instantaneous and lagged covariance matrices, the cross-spectrum, and the PSDs of the model can be obtained through algebraic operations including the Jacobian matrix. This can be done both in the homogeneous and the heterogeneous cases, and also in the presence of time delays.

Content:
Network without delays:
- StochSim_HopfNet.m : stochastic numerical simulations of the network of N hopf nodes.
- HopfModel_LNA.m :  calculates the network's statistics (covariance, the lagged-covariance, the power spectral density, and the cross-spectrum) using the linear approximation.
- run_simulations_vs_LNA : compares the model statistics obtained using numerical simulations and using the linear approximation.

Network with delays:
- StochSim_DelayedHopfNet.m : stochastic numerical simulations of the network of N hopf nodes with delayed interactions.
- DelayedHopfModel_LNA.m :  calculates the network's statistics using the linear approximation, in the presence of delayed interactions.
- run_simulations_vs_LNA : compares the model statistics obtained using numerical simulations and using the linear approximation, in the presence of delayed interactions.

- Connectome250.mat : 
	- C :  Example coupling matrix (N=250 nodes).
	- D :  Example distances matrix (N=250 nodes).

 - Connectome1000.mat : 
	- C :  Example coupling matrix (N=1000 nodes).
	- D :  Example distances matrix (N=1000 nodes).
 	- FC : Example functional connectivity matrix (N=1000 nodes). 

If you use this toolbox as part of a published academic work, please cite it as:

Ponce-Alvarez A and Deco G (2023) The Hopf whole-brain model and its linear approximation.

_________________
Other references:

Deco, G. et al. (2017) Single or multiple frequency generators in on-going brain activity: A mechanistic whole-brain model of empirical MEG data. NeuroImage 152, 538–550.
