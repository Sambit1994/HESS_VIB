# HESS_VIB
Molecular vibrational frequency analysis from Hessian
----
The program utilizes user-provided molecular coordinates (geometry.xyz) and a Hessian matrix (Hessian.dat, provided in a lower triangular matrix form) to calculate 3N-6 vibrational frequencies. The method is based on the procedure discussed in vib.pdf (see for further insights)

* The running script is VIB_SCRIPT.py, which reads the input geometry and Hessian and runs the program using the modules (MOI_module.py and Ortho_module.py) to provide the resulting frequencies. The user has the option to collect results using either QR Decomposition or the Gram-Schmidt orthogonalisation method.
  
* The MOI_module.py module generates the translated coordinates to COM and the moment of inertia.
  
* The Ortho_module.py is used to generate the transformation matrix to generate the 3N-6*3N-6 matrix.
  
* The script should (practically) work for all non-linear molecules with 3N-6 vibrational modes. However, the default code is only tested for water (H2O) and 1-chloro-2-fluoroethene(Z)(C2H2ClF). The atomic masses only for these systems are included in the default script.

* The steps are explained as comments in the script.
