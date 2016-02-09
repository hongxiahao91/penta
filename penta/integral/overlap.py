import numpy as np
import math


def overlap(alpha,beta,A,B,a,b):
#Initial Conditions
        s=np.zeros((3,7,7))
        ki=np.zeros((3,7,7))
        P=np.zeros(3)
        for i in range(3):
                P[i]=(alpha*A[i]+beta*B[i])/(alpha+beta)
                s[i][0][0]=1
                s[i][1][0]=-(A[i]-P[i])
#Recurrence Index
                for j in range(6):
                        if j>1:
                                s[i][j][0]=-(A[i]-P[i])*s[i][j-1][0]+(j-1)/(2.0*(alpha+beta))*s[i][j-2][0]
                for k in range(6):
                        if k>0:
                                for j in range(6):
                                        s[i][j][k]=s[i][j+1][k-1]+(A[i]-B[i])*s[i][j][k-1]

	EAB=math.exp(-alpha*beta/(alpha+beta)*np.dot((A-B),(A-B)))
	return EAB*(math.pi/(alpha+beta))**1.5*s[0][a[0]][b[0]]*s[1][a[1]][b[1]]*s[2][a[2]][b[2]]


def doublefactorial(n):
     if n <= 0:
         return 1
     else:
         return n * doublefactorial(n-2)
def NormCoeff(alpha,a):
	return (2*alpha/math.pi)**0.75*((4*alpha)**((a[0]+a[1]+a[2])/2.0))/(doublefactorial(2*a[0]-1)*doublefactorial(2*a[1]-1)*doublefactorial(2*a[2]-1))**0.5



def OverLap_Matrix(OrbCoeff,FCenter,CartAng):
	M=np.zeros((len(OrbCoeff),len(OrbCoeff)))
	for p in range(len(OrbCoeff)):
		for q in range(len(OrbCoeff)):
			Sum=0
			for i in range(3):
				for j in range(3):
					Sum=Sum+NormCoeff(OrbCoeff[p][i],CartAng[p])*NormCoeff(OrbCoeff[q][j],CartAng[q])*PrimCoeff[p][i]*PrimCoeff[q][j]*overlap(OrbCoeff[p][i],OrbCoeff[q][j],FCenter[p],FCenter[q],CartAng[p],CartAng[q])
			M[p][q]=Sum
	return M
#Input file
################################################################################################
R = np.array([[0.,0.,-0.7315],[0.,0.,0.7316]])
PrimCoeff = np.array([[0.444635,0.535328,0.154329],[0.444635,0.535328,0.154329]])
OrbCoeff = np.array([[0.480844,1.776691,9.753935],[0.168856,0.623913,3.42525]])
FCenter=np.array([R[0],R[1]])
CartAng=np.array([[0,0,0],[0,0,0]])
################################################################################################
