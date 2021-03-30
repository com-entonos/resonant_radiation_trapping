# resonant_radiation_trapping
legacy (f77) code for resonant radiation trapping simulations in plasmas

Uses propagator (Green's functions) for resonant radiation trapping. The spatial propagator ("Q") can be expressed as an exact divergence so surface integrals are used to construct a matrix to describe spatial redistribution of radiation. The frequency propagator ("J") can be complete or partial frequency redistribution. For partial, there are many options including Jefferies-White approximation w/ or w/o Voigt lineshape, exact angle averaged coherent or incoherent scattering. Code was used in the following papers:

J.E. Lawler, G.J. Parker and W.N.G. Hitchon, "Radiation trapping simulations using the propagator function method",  J. Quant. Spectros.  Radiat. Transfer. v49, p627 (1993).  

G.J. Parker, W.N.G. Hitchon and J.E. Lawler, "Radiation trapping simulations using the propagator function method: complete and partial frequency redistribution", J. Phys. B v26, p4643, (1993).  

A.F. Molisch, G.J. Parker, B.P. Oehry, W. Schupita and G. Magerl, "Radiation trapping with partial frequency redistribution: comparison of approximations and exact solutions", J. Quant. Spectros.  Radiat. Transfer v53, p269, (1995).

Source (f77) is in "source" directory. Makefile can use either gfortran or Intel's ifort. Examples can be found in the "examples" directory. "reference_manual" is the manual for the code while "user_manual" is the user manual.

Original tar ball can be found here: http://parker9.com/ph_prog.html#res
