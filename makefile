all:
	gfortran thermo.f90 sandiego.f90 euler2d_4.f90 weno5.f90 hllc.f90 paraview_subroutines.f90 compute_primitive.f90 drag_coeff.f90 -o sd.exe
clean:
	rm thermo.mod sd.exe
