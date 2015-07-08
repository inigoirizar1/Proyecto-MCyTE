# En este programa reescalaremos TAU para que los saltos de EMCEE sean de considerable magnitud para este parametro (mayores que 0.0001)
# Import block - clean if not needed
import os
import emcee
import triangle
import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
record_prueba = open('record_4_variables(3).txt','w')
counter = 0

def read_scattering_angles():
	lines1 = open('file1.txt','r')
	lines2 = open('file2.txt','r')
	lines3 = open('file3.txt','r')
	mu = []
	mu0 = []
	deltaphi = []
	for line in lines1:
		mu.append(float(line))
	for line in lines2:
		mu0.append(float(line))
	for line in lines3:
		deltaphi.append(float(line))
        return mu, mu0, deltaphi

def write_model(mu, mu0, deltaphi, T1, T2, P3, omega3):
	newmodel = open('model.ini','w')
	newmodel.write("'LAMBDA (MICRONS) '     0.890" + "\n")
	newmodel.write("'NUMBER OF LAYERS '     -7" + "\n")
	newmodel.write("'OUTPUT SELECTION '     'points' 45" + "\n")
	newmodel.write("'PLANET           '     'saturn'" + "\n")
	newmodel.write("'P.GRAPH. LATITUDE'     +3.5" + "\n")
	newmodel.write("'MU               '     ")
	for item in mu:
		newmodel.write(str(item) + " ")
	newmodel.write("\n")
	newmodel.write("'MU0              '     ")
	for item in mu0:
		newmodel.write(str(item) + " ")
	newmodel.write("\n")
	newmodel.write("'FI-FI0 (deg)     '     ")
	for item in deltaphi:
		newmodel.write(str(item) + " ")
	newmodel.write("\n")
	newmodel.write("'-----------------'" + "\n")
	newmodel.write("\n")
	newmodel.write("'TYPE OF LAYER    '     'gas'" + "\n")
	newmodel.write("'P_TOP (bars)     '    0.000" + "\n")
	newmodel.write("'P_BOTTOM (bars)  '    0.0200000" + "\n")
	newmodel.write("'K_CH4 (1/Km-amag)'    21.5" + "\n")
	newmodel.write("\n")
	newmodel.write("'TYPE OF LAYER    '     'mixed'" + "\n")
	newmodel.write("'P_TOP (bars)     '    0.0200000" + "\n")
	newmodel.write("'P_BOTTOM (bars)  '    0.0400000" + "\n")
	newmodel.write("'K_CH4 (1/Km-amag)'    21.5" + "\n")
	newmodel.write("'PHASE FUNCTION   '    'do_mie'" + "\n")
	newmodel.write("'TAU_EXT PARTICLES'    "+str(T1)  + "\n")
	newmodel.write("'FOURIER TERMS    '    20" + "\n")
	newmodel.write("'Re(m) & -Im(m)   '    1.430  -0.0010" + "\n")
	newmodel.write("'SIZE DISTRIBUTION'    'hansen' 0.10 0.100" + "\n")
	newmodel.write("'RANGE rmin & rmax'    0.005  3" + "\n")
	newmodel.write("\n")
	newmodel.write("'TYPE OF LAYER    '     'gas'" + "\n")
	newmodel.write("'P_TOP (bars)     '    0.0400000"+ "\n")
	newmodel.write("'P_BOTTOM (bars)  '    0.10000" + "\n")
	newmodel.write("'K_CH4 (1/Km-amag)'    21.5" + "\n")
	newmodel.write("\n")
	newmodel.write("'TYPE OF LAYER    '     'mixed'" + "\n")
	newmodel.write("'P_TOP (bars)     '    0.10000" + "\n")
	newmodel.write("'P_BOTTOM (bars)  '    0.60000" + "\n")
	newmodel.write("'K_CH4 (1/Km-amag)'    21.5" + "\n")
	newmodel.write("'PHASE FUNCTION   '    'do_mie'" + "\n")
	newmodel.write("'TAU_EXT PARTICLES'    "+ str(T2) + "\n")
	newmodel.write("'FOURIER TERMS    '    20" + "\n")
	newmodel.write("'Re(m) & -Im(m)   '    1.430  -0.0100" + "\n")
	newmodel.write("'SIZE DISTRIBUTION'    'hansen' 1.00 0.100" + "\n")
	newmodel.write("'RANGE rmin & rmax'    0.005  5" + "\n")
	newmodel.write("\n")
	newmodel.write("'TYPE OF LAYER    '     'gas'" + "\n")
	newmodel.write("'P_TOP (bars)     '    0.60000" + "\n")
	newmodel.write("'P_BOTTOM (bars)  '    "+ str(P3) + "\n")# P3 = 1.500000
	newmodel.write("'K_CH4 (1/Km-amag)'    21.5" + "\n")
	newmodel.write("\n")
	newmodel.write("'TYPE OF LAYER    '     'mixed'" + "\n")
	newmodel.write("'P_TOP (bars)     '    "+ str(P3)  + "\n")
	newmodel.write("'P_BOTTOM (bars)  '    2.4010000" + "\n")
	newmodel.write("'K_CH4 (1/Km-amag)'    21.5" + "\n")
	newmodel.write("'PHASE FUNCTION   '    'isotropic'" + "\n")
	newmodel.write("'OMEGA PARTICLES  '    "+ str(omega3) + "\n")
	newmodel.write("'TAU_EXT PARTICLES'    10.000" + "\n")
	newmodel.write("\n")
	newmodel.write("'TYPE OF LAYER    '     'gas'" + "\n")
	newmodel.write("'P_TOP (bars)     '    2.4010000" + "\n")
	newmodel.write("'P_BOTTOM (bars)  '    6" + "\n")
	newmodel.write("'K_CH4 (1/Km-amag)'    21.5" + "\n")
	newmodel.write("\n")
	newmodel.write("'END              '" + "\n")
	newmodel.write("\n")
	newmodel.close()	
        return None

def execute_atmos():
	os.system('./atmos')	
        return None

def read_obs():
        lon = []
        iparf = []
        mu0obs = []
        muobs = []
	lines11 = open('MT3-3.5N.dat','r').readlines()[1:]#Para leer longitudes y I/F
	for line in lines11:
		parts = line.split()
		lon.append(float(parts[0]))
		iparf.append(float(parts[1]))
		muobs.append(float(parts[2]))
		mu0obs.append(float(parts[3]))
        return lon, iparf, muobs, mu0obs

def read_pointsdat():
        mod = []
        londism = []
	lines4 = open("points.dat","r").readlines()[80:]#El output del 'atmos'
	i = int()
	j = int()
	k = int()
	while (i<179):
		for line in lines4:
			i += 1
			k = 1 + 4*j
			if (i == k):
				value = '0' + line.lstrip()
				mod.append(float(value))
				j += 1
	l = int()
	m = int()
	n = int()
	if (l<=177):
		for line in lines11:
			l += 1
			n = 1 + 4*m
			if (n == l):
				parts = line.split()
				londism.append(float(parts[0]))
				m += 1
        return mod,londism

def interpolate_model(lon,londism,mod):# Interpolamos el modelo.
	interpolated_y = np.interp(lon,londism,mod)	
        return interpolated_y

def create_delta(interpolated_y):
	iparf = read_obs()	
	return iparf - interpolated_y

def insert_noise(y):
    for item in y:
	item += np.random.uniform(-1.0,1.0) * item
    return y

def func(x, T1, T2, P3, omega3):
	global counter
	values = []
        mu, mu0, deltaphi = read_scattering_angles()
	write_model(mu, mu0, deltaphi, T1, T2, P3, omega3)
	execute_atmos()
        lon, iparf, muobs, mu0obs = read_obs()
        mod,londism = read_pointsdat()
        interpolated_y = interpolate_model(lon,londism,mod)
	interpolated_y_noise = insert_noise(interpolated_y)
	counter +=1 
	return interpolated_y_noise

# Define the probability function as likelihood  prior.
def lnprior(theta):
	T1, T2, P3, omega3 = theta
	if (0.0 < T1 < 10.0 and 5.0 < T2 < 20.0 and 1.0 < P3 < 2.0 and 0.0 < omega3 < 1.0):
		return 0.0
	return -np.inf

def lnlike(theta, x, y, yerr):#inserta la definicon de yerr en algun lado!! He cambiado yerr por sigma!!!!!
        valuesChi2 = []
	sigma = []
	T1, T2, P3, omega3 = theta   
	interpolated_y_noise = func(x,T1,T2,P3,omega3)    
	for k in range(0,np.size(y)):
            item = ((y[k]-interpolated_y_noise[k])/(0.05*y[k]))**2
            valuesChi2.append(item)
	suma = -0.5*(np.sum(valuesChi2))
        record_prueba.write(str(T1) + "\t" + str(T2) + "\t" + str(P3) + "\t" + str(omega3) + "\t" + str(suma) + "\n")
	return suma

def lnprob(theta, x, y, sigma):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)

# Main block
if __name__ == '__main__':
   # Choose the "true" parameters.
   T1_true = 0.1
   T2_true = 12.0
   P3_true = 1.5000
   omega3_true = 0.9
   
   # Generate some synthetic data from the model.
   x = []
   y = []
   yerr = 0.5
   lines11 = open('MT3-3.5N.dat','r').readlines()[1:]
   for line in lines11:
   	parts = line.split()
	x.append(float(parts[0]))
	#y.append(float(parts[1]))
   y = func(x,0.1,12.0,1.3,0.7)
   #y += np.abs(f_true*y) * np.random.randn(NN)
   #y += yerr * np.random.randn(NN)
   print np.size(x)
   print np.size(y)
   #print np.size(interpolated_y_noise)
   
   nll = lambda *args: -lnlike(*args)
   result = op.minimize(nll, [T1_true, T2_true, P3_true, omega3_true], args=(x, y, yerr))
   T1_ml, T2_ml, P3_ml, omega3_ml = result["x"]
   print T1_ml, T2_ml, P3_ml, omega3_ml   

   #initial_position = [0.08,10.0, 1.1, 0.7]
   # Set up the sampler.
   ndim, nwalkers = 4,50
   pos = [result["x"] + 1e-2*np.random.randn(ndim) for i in range(nwalkers)]
   sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))

   # Clear and run the production chain.
   print("Running MCMC...")
   sampler.run_mcmc(pos, 250, rstate0=np.random.get_state())
   print("Done.")

   plt.clf()
   fig, axes = plt.subplots(4, 1, sharex=True, figsize=(8, 9))
   axes[0].plot(sampler.chain[:, :, 0].T, color="k", alpha=0.4)
   axes[0].yaxis.set_major_locator(MaxNLocator(5))
   axes[0].axhline(T1_true, color="#888888", lw=2)
   axes[0].set_ylabel("$T1$")

   axes[1].plot(sampler.chain[:, :, 1].T, color="k", alpha=0.4)
   axes[1].yaxis.set_major_locator(MaxNLocator(5))
   axes[1].axhline(T2_true, color="#888888", lw=2)
   axes[1].set_ylabel("$T2$")

   axes[2].plot(sampler.chain[:, :, 1].T, color="k", alpha=0.4)
   axes[2].yaxis.set_major_locator(MaxNLocator(5))
   axes[2].axhline(P3_true, color="#888888", lw=2)
   axes[2].set_ylabel("$P3$")

   axes[3].plot(sampler.chain[:, :, 1].T, color="k", alpha=0.4)
   axes[3].yaxis.set_major_locator(MaxNLocator(5))
   axes[3].axhline(omega3_true, color="#888888", lw=2)
   axes[3].set_ylabel("$omega3$")

   fig.tight_layout(h_pad=0.0)
   fig.savefig("line-time_4_variables(3).png")

   # Make the triangle plot.
   burnin = 50
   samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))

   fig = triangle.corner(samples, labels=["$T1$", "$T2$","$P3$","$omega3$"],
                      truths=[T1_true, T2_true, P3_true, omega3_true])
   fig.savefig("line-triangle_4_variables(3).png")
