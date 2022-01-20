import requests
import csv,sys
from pylab import *
from scipy import stats
import numpy as np
from scipy.optimize import leastsq

### these imports were added only for github version of the code
import numpy as np
import matplotlib.pyplot as plt
import parsubs
import ulens_functions


### 1st part

t=np.linspace(5000,6000,1000)
t0 = 5500
te = 100
u0 = 0.1
I0 = 15.0
fs = 1.0
mag = ulens(t, t0, te, u0, I0, fs)
piEN = 0.1
piEE = -0.5
alpha = 270
delta = -30
t0par = t0
mag_par = ulensparallax(t, t0, te, u0, piEN, piEE, I0, fs, alpha, delta, t0)
plt.plot(t, mag)
plt.plot(t, mag_par)
plt.axis('tight')
plt.gca().invert_yaxis()
plt.xlabel("Time [days]")
plt.ylabel("magnitude")


#### 2nd Part: Inputs

inputfile="starBLG181.5.I.40760.dat" #you need to input photometric data from the main directory
name=inputfile
t0_0=0 #smart
te_0=60
u0_0=0.15
I0_0=0 #smart
fs_0=1.0
piEN_0 = 0.
piEE_0 = 0.0

### some examples of inputs are below

#BLG343.8.I.24509
# alpha=267.10663 
# delta=-23.38713 
#t0par=3120.33844267
# inputfile="starBLG343.8.I.24509.dat"
# name = "BLG343.8.I.24509"

#from http://gsaweb.ast.cam.ac.uk/alerts/alert/Gaia19dke/
alpha=291.49451   
delta=28.40686

#starBLG131.1.I.74150
#alpha=267.593292
#delta=-34.438667
# t0par=2821.404253
# inputfile="starBLG131.1.I.74150.dat"
# name = "starBLG131.1.I.74150"



### 3rd part

if (inputfile!=""):
  my_data = np.genfromtxt(inputfile, delimiter=' ')
  if my_data.ndim == 1:
    my_data = my_data.reshape(1, -1)
  t = my_data[:,0] ##OGLE data has truncated hjd already (-2450000)
  mag = my_data[:,1]
  err = (my_data[:,2])
else:
  #reading from web:
  name="Gaia19dke"
  url = "http://gsaweb.ast.cam.ac.uk/alerts/alert/"+name+"/lightcurve.csv/"
  req = requests.get(url)
  text = req.text
  lines = text.split('\n')
  gaiahjd = []
  gaiamag = []
  gaiamagerr = []
  #
  for l in lines:
      col = l.split(',')
      if (len(col)>1):
          if (len(col)==3 and (col[1] != 'JD')):
              if (col[2] != 'null' and col[2] != 'untrusted' ):
                  hjd = float(col[1])
                  if (hjd>2450000.0): hjd-=2450000.
                  gaiahjd.append(float(hjd))
                  gaiamag.append(float(col[2]))
                  gaiamagerr.append(0.1)

  t=np.array(gaiahjd)
  mag = np.array(gaiamag)
  err = np.array(gaiamagerr)

  
### plotting datapoints

plt.scatter(t,mag,color='k',alpha=.9,s=3)
plt.errorbar(t,mag,yerr=err, ls='none', ecolor='k')
plt.axis('tight')
plt.gca().invert_yaxis()
plt.show()



### 4th part: fitting and ploting fitted curve over input data

maxmag = mag.min()
index = mag.argmin()
smartt0 = t[index]
t0_0 = smartt0
smartI0 = np.median(mag)
I0_0 = smartI0

t0par = t0_0
v0 = [t0_0, te_0, u0_0, piEN_0, piEE_0, I0_0, fs_0]  #starting values (I0 and t0 are from data directly)
fn = lambda v, x: ulensparallax(x, v[0], v[1], v[2], v[3], v[4], v[5], v[6], alpha, delta, t0par) #values from ulens fitting
en = lambda v, x, y, err: ((fn(v,x)-y)/err) #errors
vn, successn = leastsq(en, v0, args=(t,mag,err), maxfev=1000000) #leastsquares
t0=vn[0]
te=vn[1]
u0=vn[2]
piEN=vn[3]
piEE=vn[4]
I0=vn[5]
fs=vn[6]

#print(vn)
verbose=1

chi2=sum(en(vn,t,mag,err)**2)#/(x.size-vn.size)
ndof = len(t)-7
chi2dof = chi2/ndof
if (verbose>0) :
    print(("t0    = %f ")%(t0))
    print(("tE    = %f ")%(te))
    print(("u0    = %f ")%(u0))
    print(("piEN  = %f ")%(piEN))
    print(("piEE  = %f ")%(piEE))
    print(("I0    = %f ")%(I0))
    print(("fs    = %f ")%(fs))
    print(("alpha = %f ")%(alpha))
    print(("delta = %f ")%(delta))
    print(("t0par = %f\n")%(t0par))
    print(("chi2      = %f\nndof      = %f\nchi2/ndof = %f")%(chi2, ndof, chi2dof))
    
plt.subplot(111)
t1=t.min()
t2=t.max()
npoints = 1000
XX = linspace(t1,t2,npoints)

offsetres = I0+5*err.max()
resy = fn(vn,XX) - fn(vn,XX) + offsetres
resx = mag-fn(vn,t) + offsetres

errorbar(t, mag, yerr=err, fmt='k.', alpha=0.6)
errorbar(t, resx, yerr=err, fmt='k.', alpha=0.6)

plot(XX, fn(vn,XX), 'm-')
plot(XX, resy, 'm-')

ax = gca()
ax.set_xlim([t1,t2])
plt.gca().invert_yaxis()
plt.title(name)
plt.xlabel("Julian Date - 2450000")
plt.ylabel("Magnitude")

plt.show()
