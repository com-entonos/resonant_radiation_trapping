# resonant_radiation_trapping

<b>What:</b>
The generalized Holstein-Biberman-Payne equation describes the resonant radiation temporal and spatial evolution. This may be important for energy balance, ionization source, the spectrum of escaping radiation, ...

<b>This:</b>
Uses propagator (Green's functions) for resonant radiation trapping. The spatial propagator ("Q") can be expressed as an exact divergence so surface integrals are used to construct a matrix to describe spatial redistribution of radiation. 1D, cylindrical and spherical geometries are supported. Similarly, the frequency propagator/matrix ("J") can be complete or partial frequency redistribution. For partial, there are many options including Jefferies-White approximation (Voigt/Dirac/Doppler lineshape) and coherent and incoherent exact angle averaged scattering.

<b>References:</b>
Code was used in the following papers:

J.E. Lawler, G.J. Parker and W.N.G. Hitchon, "Radiation trapping simulations using the propagator function method",  J. Quant. Spectros.  Radiat. Transfer. v49, p627 (1993).

G.J. Parker, W.N.G. Hitchon and J.E. Lawler, "Radiation trapping simulations using the propagator function method: complete and partial frequency redistribution", J. Phys. B v26, p4643, (1993).

A.F. Molisch, G.J. Parker, B.P. Oehry, W. Schupita and G. Magerl, "Radiation trapping with partial frequency redistribution: comparison of approximations and exact solutions", J. Quant. Spectros.  Radiat. Transfer v53, p269, (1995).

<b>Repository:</b>

	source/			build directory
	source/src		code (legacy f77)
	examples/		example outputs
	examples/calc_decay	find fundamental mode
	examples/calc_ss	find steady-state solution with external production
	examples/jw_figs	find effective lifetime via Jefferies-White approximation
	reference_manual	reference manual
	user_manual		user manual

<b>Build and run:</b>

	cd source && make rad EXE=gnu && cd ../examples/calc_decay/ && ../../source/rad_gnu


Original tar ball can be found here: http://parker9.com/ph_prog.html#res

