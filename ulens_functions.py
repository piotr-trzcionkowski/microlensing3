#Functions with microlensing models

import numpy as np
import matplotlib.pyplot as plt



#amplification in Paczynski curve
def ulens_A(t,t0,tE,u0):

	tau2 = ((t-t0)/tE)**2
	u02 = u0**2
	u = np.sqrt(u02 + tau2)
	A = (u**2 + 2)/(u*np.sqrt(u**2 + 4))
	return A	

#magnitude, classical lens, no effects, blending included
def ulens(t, t0, tE, u0, I0, fs):
    ampl= ulens_A(t, t0, tE, u0)
    F = ampl*fs + (1-fs)
    I = I0 - 2.5*np.log10(F)
    return I

#magnitude, classical lens, no effects, fixed blending
def ulens_fixedbl(t, t0, tE, u0, I0):
    fs=1.
    ampl= ulens_A(t, t0, tE, u0)
    F = ampl*fs + (1-fs)
    I = I0 - 2.5*np.log10(F)
    return I

#amplification in parallax model
def ulensparallax_A(t, t0, tE, u0, piEN, piEE, alpha, delta, t0par):
  qe=[]
  qn=[]
  for hjd in t:
    qn_single,qe_single=geta(hjd, alpha, delta, t0par) 
    qn.append(qn_single)
    qe.append(qe_single)
  qn = np.array(qn)
  qe = np.array(qe)

  dtau = piEN*qn + piEE*qe
  dbeta = piEN*qe - qn*piEE

  tau = (t-t0)/tE + dtau
  beta = u0 + dbeta
  u = np.sqrt(tau**2 + beta**2)
  A = (u**2 + 2)/(u*np.sqrt(u**2 + 4))
  return A 

#magnitude in parallax model
def ulensparallax(t, t0, tE, u0, piEN, piEE, I0, fs, alpha, delta, t0par):
  ampl = ulensparallax_A(t, t0, tE, u0, piEN, piEE, alpha, delta, t0par)
  F = ampl*fs + (1-fs)
  I = I0 - 2.5*np.log10(F)

  return I 
