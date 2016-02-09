import numpy as np
import math
from scipy.integrate import quad
def Nea(alpha,beta,A,B,R,a,b):
#Initial Conditions
	eta=np.zeros((3,7,7))
	P=np.zeros(3)
	def integrand(t):
		for i in range(3):
			P[i]=(alpha*A[i]+beta*B[i])/(alpha+beta)
			eta[i][0][0]=1
			eta[i][1][0]=-(A[i]-P[i]+t**2*(P[i]-R[i]))
#Recurrence Index
			for j in range(6):
				if j >1:
					eta[i][j][0]=-(A[i]-P[i]+t**2*(P[i]-R[i]))*eta[i][j-1][0]+(j-1)/(2.0*(alpha+beta))*(1-t**2)*eta[i][j-2][0]
		
#Transfer Equation
			for k in range(6):
				if k>0:
					for j in range(6):


							eta[i][j][k]=eta[i][j+1][k-1]+(A[i]-B[i])*eta[i][j][k-1]
		return 0.5*eta[0][a[0]][b[0]]*eta[1][a[1]][b[1]]*eta[2][a[2]][b[2]]*math.exp(-(alpha+beta)*t**2*np.dot((P-R),(P-R)))
#Exponential Overlap
        EAB=math.exp(-alpha*beta/(alpha+beta)*np.dot((A-B),(A-B)))

	return EAB*2.0*math.pi/(alpha+beta)*quad(integrand,-1,1)[0]
		






def doublefactorial(n):
     if n <= 0:
         return 1
     else:
         return n * doublefactorial(n-2)
def NormCoeff(alpha,a):
        return (2*alpha/math.pi)**0.75*((4*alpha)**((a[0]+a[1]+a[2])/2.0))/(doublefactorial(2*a[0]-1)*doublefactorial(2*a[1]-1)*doublefactorial(2*a[2]-1))**0.5



def NE_Matrix(OrbCoeff,FCenter,CartAng):
        M=np.zeros((len(OrbCoeff),len(OrbCoeff)))
        for p in range(len(OrbCoeff)):
                for q in range(len(OrbCoeff)):
                        Sum=0
                        for Nuc in range(len(Z)):
                                for i in range(3):
                                        for j in range(3):
                                                Sum=Sum+-Z[Nuc]*NormCoeff(OrbCoeff[p][i],CartAng[p])*NormCoeff(OrbCoeff[q][j],CartAng[q])*PrimCoeff[p][i]*PrimCoeff[q][j]*Nea(OrbCoeff[p][i],OrbCoeff[q][j],FCenter[p],FCenter[q],R[Nuc],CartAng[p],CartAng[q])
				
                                M[p][q]=Sum
        return M

################################################################################################
R = np.array([[0.,0.,-0.7315],[0.,0.,0.7316]])
PrimCoeff = np.array([[0.444635,0.535328,0.154329],[0.444635,0.535328,0.154329]])
OrbCoeff = np.array([[0.480844,1.776691,9.753935],[0.168856,0.623913,3.42525]])
FCenter=np.array([R[0],R[1]])
CartAng=np.array([[0,0,0],[0,0,0]])
Z=np.array([2,1])
################################################################################################
