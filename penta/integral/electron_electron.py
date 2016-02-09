import numpy as np
import math
from scipy.integrate import quad
def EE(alpha,beta,gamma,delta,A,B,C,D):
#Initial Conditions
	P=np.zeros(3)
	Q=np.zeros(3)
	def integrand(t):
		for i in range(3):
			P[i]=(alpha*A[i]+beta*B[i])/(alpha+beta)
			Q[i]=(gamma*C[i]+delta*D[i])/(gamma+delta)
		return math.exp(-((alpha+beta)*(gamma+delta)/(alpha+beta+gamma+delta)*np.dot((P-Q),(P-Q)))*t**2)
		
#Exponential Overlap
        EAB=math.exp(-alpha*beta/(alpha+beta)*np.dot((A-B),(A-B)))
        ECD=math.exp(-gamma*delta/(gamma+delta)*np.dot((C-D),(C-D)))

	return 2*math.pi**2.5/((alpha+beta)*(gamma+delta)*(alpha+beta+gamma+delta)**0.5)*EAB*ECD*quad(integrand,0,1)[0]
		






def doublefactorial(n):
     if n <= 0:
         return 1
     else:
         return n * doublefactorial(n-2)
def NormCoeff(alpha,a):
        return (2*alpha/math.pi)**0.75*((4*alpha)**((a[0]+a[1]+a[2])/2.0))/(doublefactorial(2*a[0]-1)*doublefactorial(2*a[1]-1)*doublefactorial(2*a[2]-1))**0.5



def EE_Matrix(OrbCoeff,FCenter,CartAng):
	M=np.zeros((len(OrbCoeff),len(OrbCoeff),len(OrbCoeff),len(OrbCoeff)))
	for a in range(len(OrbCoeff)):
		for b in range(len(OrbCoeff)):
			for c in range(len(OrbCoeff)):
				for d in range(len(OrbCoeff)):
        				Sum=0
        				for i in range(3):
        					for j in range(3):
							for k in range(3):
								for l in range(3):
        		        					Sum=Sum+NormCoeff(OrbCoeff[a][i],CartAng[a])*NormCoeff(OrbCoeff[b][j],CartAng[b])*NormCoeff(OrbCoeff[c][k],CartAng[c])*NormCoeff(OrbCoeff[d][l],CartAng[d])*PrimCoeff[a][i]*PrimCoeff[b][j]*PrimCoeff[c][k]*PrimCoeff[d][l]*EE(OrbCoeff[a][i],OrbCoeff[b][j],OrbCoeff[c][k],OrbCoeff[d][l],FCenter[a],FCenter[b],FCenter[c],FCenter[d])
					M[a][b][c][d]=Sum
        return M

################################################################################################
R = np.array([[0.,0.,-0.7315],[0.,0.,0.7316]])
PrimCoeff = np.array([[0.444635,0.535328,0.154329],[0.444635,0.535328,0.154329]])
OrbCoeff = np.array([[0.480844,1.776691,9.753935],[0.168856,0.623913,3.42525]])
FCenter=np.array([R[0],R[1]])
CartAng=np.array([[0,0,0],[0,0,0]])

Z=np.array([2,1])
################################################################################################

