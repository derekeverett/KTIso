# KTIso



KTIso is a code to numerically propagate the stress tensor from asymptotically 
early times to a hydrodynamic matching time in heavy ion collisions.

The equations of motion have been derived from the Boltzmann transport equation with the Isotropization Time
Approximation for the collision kernel [https://arxiv.org/abs/1905.05139](https://arxiv.org/abs/1905.05139). 

This code has only been tested for a specific case of a boost-invariant distribution proportional to a delta function in $p^z$ (particles which can move and scatter in the tranverse plane, but can not have a longitudinal momentum component). 

It's purpose is to act as a pre-hydrodynamic dynamic model with a free parameter that controls the collision rate. In the ITA eqns of motion, collisions isotropize momentum space on a time scale given by $\tau_{\rm iso}=\frac{\eta}{sT}$, where $\eta/s$ is a constant specific shear viscosity, and $T$ the local temperature. 