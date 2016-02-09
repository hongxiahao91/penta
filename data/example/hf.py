import numpy as np
import math
from scipy.integrate import quad
from overlap import OverLap_Matrix
from kinetic import Kinetic_Matrix
from nuclear_electron import NE_Matrix
from electron_electron import EE_Matrix
from read_xyz import load_xyz 
def Inital_Coefficient_Matrix(OrbCoeff,FCenter,CartAng):
	S=OverLap_Matrix(OrbCoeff,FCenter,CartAng)
	C=np.zeros((len(OrbCoeff),len(OrbCoeff)))
	C[0][0]=(2*(1+S[0][1]))**-0.5
	C[0][1]=(2*(1-S[0][1]))**-0.5
	C[1][0]=(2*(1+S[0][1]))**-0.5
	C[1][1]=-(2*(1-S[0][1]))**-0.5
#	M_Coeff=np.eye(len(OrbCoeff))
	return C
def Density_Matrix():
	pass
def Fock_Matrix():
	pass
def HF_solver(OrbCoeff,FCenter,CartAng):
	S=OverLap_Matrix(OrbCoeff,FCenter,CartAng)
	T=Kinetic_Matrix(OrbCoeff,FCenter,CartAng)
	NE=NE_Matrix(OrbCoeff,FCenter,CartAng)
	H=T+NE
	EE=EE_Matrix(OrbCoeff,FCenter,CartAng)
	C=np.zeros((len(OrbCoeff),len(OrbCoeff)))
	P=np.zeros((len(OrbCoeff),len(OrbCoeff)))
	G=np.zeros((len(OrbCoeff),len(OrbCoeff)))
	F=np.zeros((len(OrbCoeff),len(OrbCoeff)))
####################################################
	C=Inital_Coefficient_Matrix(OrbCoeff,FCenter,CartAng)
	for u in range(len(OrbCoeff)):
		for v in range(len(OrbCoeff)):
			p=0.0
			for a in range(int(len(OrbCoeff)/2.0)):
				p=p+2*C[u][a]*C[v][a]
			P[u][v]=p
	for i in range(len(OrbCoeff)):
		for j in range(len(OrbCoeff)):
			g=0.0
			for k in range(len(OrbCoeff)):
				for l in range(len(OrbCoeff)):
					g=g+P[i][j]*(EE[i][j][k][l]-0.5*EE[i][l][k][j])	
			G[i][j]=g
	F=H+G
	return F

################################################################################################



mol=load_xyz('./h2.xyz')
R=mol['coordinates']
FCenter=np.array([R[0],R[1]])
Z=mol['numbers']
#obasis=load_basis('sto-3g')
PrimCoeff = np.array([[0.444635,0.535328,0.154329],[0.444635,0.535328,0.154329]])
OrbCoeff = np.array([[0.168856,0.623913,3.42525],[0.168856,0.623913,3.42525]])
CartAng=np.array([[0,0,0],[0,0,0]])
################################################################################################

print HF_solver(OrbCoeff,FCenter,CartAng)
