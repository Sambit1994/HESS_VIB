# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 16:42:08 2020

@author: SAMBIT KUMAR DAS
"""

import numpy as np
import MOI_module as MOI
import Ortho_module as Ortho

### Constants used in the calculations

amu_t_kg = 1.660538782e-27
mass_e = 9.10938215e-31

### Reading the geometry.xyz file
### coordinates are saved in 'At_coord'
### Elements are saved in 'Atms' and used to find the atomic-masses in amu

At_coord = []
Atms = []

with open('geometry.xyz') as infile:
    for line in infile:
        At_coord.append(list(map(float,line.split()[1:4])))
        Atms.append(line.split()[0])

Natms = len(At_coord)

At_coord = (np.array(At_coord)).reshape(Natms,3)
Atms = np.array(Atms)

### Atomic masses are collected in 'M'

M = []

for i  in range(Natms):
    #print (Atms[i])
    if Atms[i] == 'H':
        mass = float(1.0079)
    elif Atms[i] == 'D':
        mass = float(2.0104)
    elif Atms[i] == 'O':
        mass = float(15.9994)
    elif Atms[i] == 'C':
        mass = float(12.0107)
    elif Atms[i] == 'Cl':
        mass = float(34.9688527)
    elif Atms[i] == 'F':
        mass = float(18.9984032)
    else:
        mass = float(0.0)
        print ('The atomic masses of only the above elements are present')
    M.append(float(mass))

M = np.array(M)

### XYZ coordinates and Atomic masses are supplied to MOI module to generate the following:
### Rotational constants in cm-1: RC_vals_cmi
### Rotational constants in GHz: RC_vals_GHZ
### Eigen vectors of MOI matrix: evecs
### Translated XYZ coordinates to COM: RT

RC_vals_cmi, RC_vals_GHZ, evecs, RT = MOI.print_function(At_coord, M)

### collecting Rotational constants(in GHz) to 'Rt_constants'

Rt_constants = np.round(RC_vals_GHZ, decimals=5)

### Atomic masses, Translated XYZ and Eigen vectors of MOI are supplied to Orthogonalisaton module:
### Matrix used for transforming original mass weighted hessian matrix is retured: D_matrix; eqn.6 in vib.pdf

D_matrix = Ortho.print_function(M, RT, evecs)

### Reading the Hessian.dat file and creating and 3N*3N matrix
### The procedure first collect the data in a lower triangular matrix 'LTM' of dimension 3N*3N
### Then complete symmetric matrix is created and saved as 'Hess_matrix'

arr = []
with open('Hessian.dat') as infile:
    for line in infile:
        arr.append(list(map(float,line.split()[0:(3*Natms)])))
        

## nHmatrix is same as 3*Natms
nHmatrix = len(arr)

matrix = np.zeros([nHmatrix, nHmatrix]).astype(float)

for i in range(nHmatrix):
    ele = arr[i]
    for j in range(len(ele)):
        matrix[i,j] = ele[j] 

LTM = matrix
UTM = LTM.transpose()

## Using numpy.triu function the upper triangular elements of the matrix are found

UTM_ele = UTM[np.triu_indices(nHmatrix, k=1)]
idx = np.triu_indices(nHmatrix, k=1)
LTM[idx] = UTM_ele

Hess_matrix = LTM

### 3N*3N Hessian matrix has been created: Hess_matrix

### creating 3N*3N mass-wted Hessian matrix; 'Hess_matrix_wted'; sec.2.1 in vib.pdf
### The atomic masses in electron units are collected in 'Me'

Me = []

for i in range(Natms):
    for j in range(3):
        #print (M[i])
        Me.append(M[i])
        
## converting mass units from amu to kg to electron mass

Me = np.array(Me) * amu_t_kg
Me = np.array(Me) / mass_e

Hess_matrix_wted = []
for i in range(nHmatrix):
    for j in range(nHmatrix):
        #print (Hess_matrix[i,j], M[i], M[j])
        wted = Hess_matrix[i,j] / np.sqrt(Me[i]*Me[j])
        Hess_matrix_wted.append(float(wted))


arr = np.array(Hess_matrix_wted).astype(float)

MWH_matrix = arr.reshape(nHmatrix,nHmatrix) 

### mass weighted Hessian matrix has been created: 'MWH_matrix'

### Finding the eigen values of MWH_matrix "ev" and converting them to vibrational frequencied (in cm-1) 'vib_freqs'

e_val = np.linalg.eig(MWH_matrix)
ev = np.array(e_val[0])
ev_abs = np.abs(ev)

au_t_cmi   = 2.194746315e+05

## converting the eigen values to vibrational freqencies in cm-1

vibs = (np.sqrt(ev_abs)) * au_t_cmi
vibs = vibs*(np.sign(ev))
vib_freqs = np.round(vibs, decimals=5)    

### 'vib_freqs' are the eigen values (vibrational frequencies in cm-1) of mass weighted Hessian matrix; sec-2 of vib.pdf

### Tranforming the 3N*3N matrix to 3N-6*3N-6 matrix (Internal Coordinates) using 'D_matrix' obtained from Ortho_module
### 'INT_matrix' is the 3N-6*3N-6 matrix; sec:2.4 of vib.pdf
### The eigen values of the 'INT_matrix' is collected in 'ev_INT'
### The eigen values are converted into vibrational frequencies(in cm-1) and collected in 'vib_freqs_INT'

D_matrix_T = D_matrix.transpose()
Transform1 = np.matmul(MWH_matrix,D_matrix)
Transform2 = np.matmul(D_matrix_T,Transform1)

INT_matrix = Transform2

e_val_INT = np.linalg.eig(INT_matrix)

ev_INT = np.array(e_val_INT[0])
ev_INT_abs = np.abs(ev_INT)

vibs_INT = (np.sqrt(ev_INT_abs)) * au_t_cmi
vibs_INT = vibs_INT*(np.sign(ev_INT))
vib_freqs_INT = np.round(vibs_INT, decimals=5)

### The eigen vectors of the Transformed matrix are collected in 'evec_INT'
### The eigen vectors are used to find the Reduced masses; sec:2.6 of vib.pdf

evec_INT = np.array(e_val_INT[1])
evec_INT_T = evec_INT.transpose()

## To check the obtained eigen values of 'INT_matrix' are same as the diagonal elements

p1 = np.dot(INT_matrix,evec_INT)
p2 = np.matmul(evec_INT_T,p1)
p3 = np.abs(p2)

vibs_INT_dia = (np.sqrt(p3))*au_t_cmi
vibs_INT_dia = vibs_INT_dia*(np.sign(p2))
#print (vibs_INT_dia)

### Finding the Reduced masses of the vibrational modes
### 'L_matrix' is the eigen vector matrix of the  'INT_matrix', i.e., same as 'evec_INT'
### 'DL' matrix is created first; eqn.9 of the vib-file
### 'IL' diagonal matrix is created where the diagonal elements are sqrt(atomic-masses); eqn.10 of the vib.pdf
### 'L_cart' matrix is created; eqn.11 of the vib-file to generate the Redcuced masses
### 'Red_mass' contains the reduced masses of the vibrational modes (in amu)

L_matrix = evec_INT
DL = np.dot(D_matrix,L_matrix)

TNmol = len(Me)

IL = np.zeros([TNmol,TNmol])
for i in range(TNmol):
    IL[i,i] = 1/np.sqrt(Me[i])
    
    
L_cart = np.dot(IL,DL)
L_cart = L_cart**2

RM_vals = []
for i in range((3*Natms)-6):
   RM_val =  sum(np.array(L_cart[:,i]))
   RM_val = (((1/RM_val)*mass_e)/amu_t_kg)
   RM_vals.append(RM_val)

Red_mass = np.round(RM_vals, decimals=5)

print ('\n')
print ('The Rotational constants(in GHz) are:')
print (*Rt_constants, sep = ', ')
print ('\n')
print ('The Vibratioanl Frequencies(in cm-1) obtained from the mass weighted Hessian Matrix are:') 
print (*vib_freqs, sep = ', ')
print ('\n')
print ('The Vibratioanl Frequencies(in cm-1) obtained after transforming the mass weighted Hessian Matrix to internal coordinates are:') 
print (*vib_freqs_INT, sep = ', ')
print ('\n')
print ('The corresponding reduced masses of the vibrational modes are:')
print (*Red_mass, sep = ', ')
