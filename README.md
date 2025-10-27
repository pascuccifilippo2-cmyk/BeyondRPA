gap_number_YINS.f90 contains the Fortran90 code to self-consistently solve the gap and number equation at zero temperature for an electron-hole bilayer system with the inclusion of the Hartree-Fock term and without a screening process. 
RPA.nb contains the Mathematica 14 code that evaluates all the first-order corrections to the normal and anomalous polarization function for a fixed density using the data produced by the gap_number_YINS.f90 code in the output "Deltakvkukd[value]YINS.txt"

Run the gap_number_YINS.f90 to generate the "Deltakvkukd[value]YINS.txt" file which contains on the first column the wave-vector k array and on the second the corresponding Deltak superfluid gap value. 
Then, in the Mathematica code, insert the physical parameter of the system used in the Fortran code and the path of "Deltakvkukd[value]YINS.txt" file where required. Compile the Mathematica code and the value of the different polarization functions for a discretised interval of k is generated.  
