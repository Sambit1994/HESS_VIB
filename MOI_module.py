# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 14:08:57 2020

@author: SAMBIT KUMAR DAS
"""
def print_function( par1, par2 ):
    
    import numpy as np

### Constants used in the calculations

    h = 6.626068960e-27 # in erg s
    c = 29979245800 # in cm s-1
    amu_t_g  = 1.660538782e-24 # amu to gram
    ang_t_cm = 1e-8 # angstorm to cm
    amu_ang2_t_g_cm2 = amu_t_g * (ang_t_cm**2) # amu(ang**2) to g(cm**2)
    cmi_t_GHZ = 29.9792458 # cm-1 to GHz
    rc_c = h/(8*((np.pi)**2)*c) # value of h/[8(pi**2)c] 

### Taking input from the main running Vibration script
### par1: intial X, Y, Z coodinates in Angstrom; saved in 'R_ini'
### par2: Mass of atoms in amu; saved in 'M'
    
    R_ini = np.array(par1)
    M = np.array(par2)
    Natms = len(M)
    
    R = (R_ini.reshape(Natms,3)).transpose()

### Finding the COM of the molecule and translating the coordinates; eqn.3 in vib.pdf
    
    M_sum = 0
    for i in range(Natms):
        M_sum = M_sum + float(M[i])

    X_COM = float(0.0)
    Y_COM = float(0.0)
    Z_COM = float(0.0)
    for i in range(Natms):
        X_COM = X_COM + M[i]*(R[0,i])
        Y_COM = Y_COM + M[i]*(R[1,i])
        Z_COM = Z_COM + M[i]*(R[2,i])

    X_COM = X_COM / M_sum
    Y_COM = Y_COM / M_sum
    Z_COM = Z_COM / M_sum
    
    #COM = np.array([X_COM, Y_COM, Z_COM])
    
    R_Translated = []
    
    for i in range(Natms):
        X_shift = R[0,i] - X_COM
        Y_shift = R[1,i] - Y_COM
        Z_shift = R[2,i] - Z_COM
        R_Translated.append(float(X_shift))
        R_Translated.append(float(Y_shift))
        R_Translated.append(float(Z_shift))
        
    R_Translated = np.array(R_Translated).reshape(Natms,3)      
    RT = R_Translated
 
### Coordinates have been translated to the COM
### R_Translated/RT are the translated coordinates

        
### Creating the Moment of Inertia matrix (I); eqn.4 in vib.pdf
    
    I_0_0 = float(0.0)
    I_1_1 = float(0.0)
    I_2_2 = float(0.0)
    I_0_1 = float(0.0)
    I_0_2 = float(0.0)
    I_1_2 = float(0.0)
        
    for i in range(Natms):
        I_0_0 = I_0_0 + (M[i]*(((RT[i,1])**2) + ((RT[i,2])**2))) 
        I_1_1 = I_1_1 + (M[i]*(((RT[i,0])**2) + ((RT[i,2])**2)))        
        I_2_2 = I_2_2 + (M[i]*(((RT[i,0])**2) + ((RT[i,1])**2)))        
        I_0_1 = I_0_1 + (M[i]*(RT[i,0]*RT[i,1])) 
        I_0_2 = I_0_2 + (M[i]*(RT[i,0]*RT[i,2]))
        I_1_2 = I_1_2 + (M[i]*(RT[i,1]*RT[i,2]))
        
    I_0_1 = -(I_0_1)
    I_0_2 = -(I_0_2)
    I_1_2 = -(I_1_2)
            
    I = np.zeros((3,3)).astype(float)
    I[0,0] = I_0_0
    I[1,1] = I_1_1
    I[2,2] = I_2_2
    I[0,1] = I_0_1
    I[1,0] = I[0,1]
    I[0,2] = I_0_2
    I[2,0] = I[0,2]
    I[1,2] = I_1_2
    I[2,1] = I[1,2]              
    
### Moment of Inertia matrix (I) has been created
### Finding the Eigen values (evals) and eigen vectors (evecs) of I matrix    

    evals, evecs = np.linalg.eig(I)

### Finding the rotationals constants in cm-1(RC_vals_cmi) and GHz(RC_vals_GHZ)
    
    RC_vals_cmi = rc_c/(np.array(evals)*(amu_ang2_t_g_cm2))
    RC_vals_GHZ = RC_vals_cmi*cmi_t_GHZ
    
    return (RC_vals_cmi, RC_vals_GHZ, evecs, RT)