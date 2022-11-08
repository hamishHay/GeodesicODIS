import numpy as np 
import matplotlib.pyplot as plt 


R = 1560e3
Omega = 2.05e-5
w = -Omega 
g = 1.3
H = 10e3
obliq = np.deg2rad(-1.0)
alpha = 1e-8
m=1
t = 0.01 * 2*np.pi/w

lam = w/(2*Omega)
lambs = 4*Omega**2*R**2 / (g*H)

U21 = 0.5 * Omega**2 * R**2 * obliq

def pn(n):
    return (n+1)*(n+m)/(n*(2*n+1))

def qn(n):
    return n*(n+1-m)/((n+1)*(2*n+1))

def Kn(n):
    return lam + m/(n*(n+1)) - n*(n+1)/(lambs * lam) + 1j*alpha / (2*Omega)
    
def Ln(n):
    return lam + m/(n*(n+1)) + 1j*alpha / (2*Omega)

Psi11 = U21/(2*Omega) / (qn(1) - Kn(2)*Ln(1)/pn(2))
Phi21 = -1j*Ln(1)/pn(2) * Psi11

lats = np.arange(0,180.1,0.1)*np.pi/180.0 
lons = np.arange(0, 360, 0.1)*np.pi/180.0

lons,lats = np.meshgrid(lons,lats)

Y31 = -1.5 * (5*np.cos(lats)**2 - 1)*np.sin(lats)*np.exp(1j * 1 * lons)

Y21 = -3 * np.cos(lats)*np.sin(lats) * np.exp(1j * 1 * lons)

Y11 = -np.sin(lats)*np.exp(1j*1*lons)

U = (-1/3. * (Psi11 * Y21) + 1j*Phi21*Y21)*np.exp(-1j*w*t)
dUdt = -1j*w*U

U += U.real - 1j*U.imag         # Add complex conjugate
U /= 2*R*np.sin(lats)


V = (1j*Psi11 *Y11 + Phi21*(4*Y31 - 9*Y11)/5.)*np.exp(-1j*w*t)
V += V.real - 1j*V.imag         # Add complex conjugate 
V /= 2*R*np.sin(lats)

eta = -1j*6/R**2 * H/w * Phi21 * Y21  * np.exp(-1j*w*t)
eta += eta.real - 1j*eta.imag
eta *= 0.5

# Calc eta a cell centers
# Calc U, V at cell edges

skip = 100

plt.contourf(lons,np.pi*0.5 - lats,eta)
plt.colorbar()
plt.quiver(lons[::skip,::skip],np.pi*0.5 -lats[::skip,::skip],U[::skip,::skip],-V[::skip,::skip])
plt.savefig("init.test.png")


