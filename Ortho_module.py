# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 16:42:57 2020

@author: SAMBIT KUMAR DAS
"""

def print_function( par1, par2, par3 ):

    import numpy as np

### Taking input from the running Vibration script
### par1: Mass of atoms in amu; saved in 'M'
### par2: Translated X, Y, Z coordinates; saved in 'RT'
### par3: eigen vectors of the MOI matrix (I); saved in 'X'
    
    M = np.array(par1)
    RT = np.array(par2)
    X = np.array(par3)
    Natms = len(M)
    
### creating Translation Vectors; D1, D2 and D3; sec:2.3 of vib.pdf

    D1 = []
    D2 = []
    D3 = []
    for i in range(Natms):
        D1.append(np.sqrt(M[i]))
        D1.append(float(0.0))
        D1.append(float(0.0))
        D2.append(float(0.0))
        D2.append(np.sqrt(M[i]))
        D2.append(float(0.0))
        D3.append(float(0.0))
        D3.append(float(0.0))
        D3.append(np.sqrt(M[i]))
    
    
    D1 = (np.array(D1).reshape((3*Natms),1))
    D1 = D1.dot(1/(np.linalg.norm(D1)))
    
    D2 = (np.array(D2).reshape((3*Natms),1))
    D2 = D2.dot(1/(np.linalg.norm(D2)))
    
    D3 = (np.array(D3).reshape((3*Natms),1))
    D3 = D3.dot(1/(np.linalg.norm(D3)))
    
### Translation vectors has been created

### Creating Rotational Vectors; D4, D5 and D6; sec_2.3 of vib.pdf 

    D4 = []
    D5 = []
    D6 = []

    for i in range(Natms):
        for j in range(3):
            DR1 = ( (RT[i,1]*X[j,2]) - (RT[i,2]*X[j,1]) )* (np.sqrt(M[i]))
            D4.append(DR1)
            DR2 = ( (RT[i,2]*X[j,0]) - (RT[i,0]*X[j,2]) )* (np.sqrt(M[i]))
            D5.append(DR2)
            DR3 = ( (RT[i,0]*X[j,1]) - (RT[i,1]*X[j,0]) )* (np.sqrt(M[i]))    
            D6.append(DR3)
        

    D4 = (np.array(D4)).reshape((3*Natms),1)
    D4 = D4.dot(1/(np.linalg.norm(D4)))
    
    D5 = (np.array(D5).reshape((3*Natms),1))
    D5 = D5.dot(1/(np.linalg.norm(D5)))
    
    D6 = (np.array(D6).reshape((3*Natms),1))
    D6 = D6.dot(1/(np.linalg.norm(D6)))

### Rotational vectors has been created

### Creating 3N-6 vectors using orthogonalisation; sec:2.3 of vib.pdf
### The procedure first creates a random matrix 'D' of size 3N*3N
### The first 6 columns are substituted with nomalised vectors from D1-D6; keeping the last 3 as random
### The matrix is then orthogonalised with both QR Decomposition(QR) and Gram-Schmidt orthogonalisation(GS)
### After the orthogonalisation the last 3N-6 columns are collected to D_matrix
### The user can select either to choose the D_matrix obtained from QR or GS technique

    D = np.random.rand((3*Natms),(3*Natms))
    D[:,0] = D1[:,0]
    D[:,1] = D2[:,0]
    D[:,2] = D3[:,0]
    D[:,3] = D4[:,0]
    D[:,4] = D5[:,0]
    D[:,5] = D6[:,0]

### QR Decomposition

    q, r1 = np.linalg.qr(D)

#D7_QR, D8_QR, and D9_QR are the last 3N-6 orthogonal vectors; collected from D_matrix (valid only for 3 atomic systems like H2O)

    D7_QR = q[:,6]
    D8_QR = q[:,7]
    D9_QR = q[:,8]

    D_matrix_QR = (q.transpose())[6:(3*Natms)]
    D_matrix_QR = D_matrix_QR.transpose()

### QR orthogonalisation done; D_matrix is D_matrix_QR

### GS orthogonalisation

    V = []
    for i in range(len(D)):
        if i == 0:
            v1 = D[:,0]
            V.append(v1)
        elif i > 0 :
            x = D[:,i]
            for j in range(len(V)):
                v2 = x
                pdt1 = np.dot(V[j],x)
                pdt2 = np.dot(V[j],V[j])
                v2 = v2 - (pdt1/pdt2)*V[j]
                x = v2
            V.append(v2)

    V_matrix = []
    for i in range(len(V)):
        norm = np.linalg.norm(V[i])
        Vn = V[i]*(1/norm)
        V_matrix.append(Vn)

    V_matrix = (np.array(V_matrix).reshape(3*Natms,3*Natms)).transpose()

#D7_ortho, D8_ortho, and D9_ortho are the last 3N-6 orthogonal vectors; collected from D_matrix (valid only for 3 atomic systems like H2O)    
    
    D7_ortho = V_matrix[:,6]    
    D8_ortho = V_matrix[:,7]    
    D9_ortho = V_matrix[:,8]
    
    D_matrix_GS = (V_matrix.transpose())[6:(3*Natms)]
    D_matrix_GS = D_matrix_GS.transpose()
    
### GS orthogoanlisation over; D_matrix is D_matrix_GS 

    print ('Kindly choose either of the following two methods to get the Transformation matrix (D_matrix) and continue')    
    print ('GS: Gram-Schmidt Orthogonalization procedure')
    print ('QR: QR Decomposition in Python')

    ON_method = input('Choose Orthogonalisation method, GS or QR: ')
    
    if ON_method == str('GS'):
        D_matrix = D_matrix_GS
    elif ON_method == str('QR'):
        D_matrix = D_matrix_QR

    return (D_matrix)