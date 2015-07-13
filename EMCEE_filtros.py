# En este programa reescalaremos TAU para que los saltos de EMCEE sean de considerable magnitud para este parametro (mayores que 0.0001)
# Import block - clean if not needed
import os
import emcee
import triangle
import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator

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

def write_model(mu, mu0, deltaphi, T1, P2, T2, P3, value1, value2):
	newmodel = open('model.ini','w')
	newmodel.write("'LAMBDA (MICRONS) '     "+ str(value1) + "\n")
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
	newmodel.write("'P_BOTTOM (bars)  '    0.00100" + "\n")
	newmodel.write("'K_CH4 (1/Km-amag)'    "+ str(value2) + "\n")
	newmodel.write("\n")
	newmodel.write("'TYPE OF LAYER    '     'mixed'" + "\n")
	newmodel.write("'P_TOP (bars)     '    0.0010000" + "\n")
	newmodel.write("'P_BOTTOM (bars)  '    "+ str(P2) + "\n")# P2 = 0.04
	newmodel.write("'K_CH4 (1/Km-amag)'    "+ str(value2) + "\n")
	newmodel.write("'PHASE FUNCTION   '    'do_mie'" + "\n")
	newmodel.write("'TAU_EXT PARTICLES'    "+str(T1)  + "\n")# T1 = 0.05
	newmodel.write("'FOURIER TERMS    '    20" + "\n")
	newmodel.write("'Re(m) & -Im(m)   '    1.43  -0.0010" + "\n")
	newmodel.write("'SIZE DISTRIBUTION'    'hansen' 0.10 0.100" + "\n")
	newmodel.write("'RANGE rmin & rmax'    0.05  3" + "\n")
	newmodel.write("\n")
	newmodel.write("'TYPE OF LAYER    '     'gas'" + "\n")
	newmodel.write("'P_TOP (bars)     '    0.0400000"+ "\n")
	newmodel.write("'P_BOTTOM (bars)  '    "+ str(P3) + "\n")
	newmodel.write("'K_CH4 (1/Km-amag)'    "+ str(value2) + "\n")
	newmodel.write("\n")
	newmodel.write("'TYPE OF LAYER    '     'mixed'" + "\n")
	newmodel.write("'P_TOP (bars)     '    "+ str(P3) + "\n")# P3 = 0.1
	newmodel.write("'P_BOTTOM (bars)  '    0.50000" + "\n")
	newmodel.write("'K_CH4 (1/Km-amag)'    "+ str(value2) + "\n")
	newmodel.write("'PHASE FUNCTION   '    '2hg'" + "\n")
        newmodel.write("'OMEGA PARTICLES  '    0.990" + "\n")# omega = 0.95
	newmodel.write("'TAU_EXT PARTICLES'    "+ str(T2) + "\n")# tau2 = 15.0
	newmodel.write("'FOURIER TERMS    '    20" + "\n")
	newmodel.write("'f & g1 & g2      '    0.7 0.7 -0.3" + "\n")
	newmodel.write("\n")
	newmodel.write("'TYPE OF LAYER    '     'gas'" + "\n")
	newmodel.write("'P_TOP (bars)     '    0.50000" + "\n")
	newmodel.write("'P_BOTTOM (bars)  '    1.0000" + "\n")# P3 = 1.500000
	newmodel.write("'K_CH4 (1/Km-amag)'    "+ str(value2) + "\n")
	newmodel.write("\n")
	newmodel.write("'TYPE OF LAYER    '     'mixed'" + "\n")
	newmodel.write("'P_TOP (bars)     '    1.0000"  + "\n")
	newmodel.write("'P_BOTTOM (bars)  '    1.4000" + "\n")
	newmodel.write("'K_CH4 (1/Km-amag)'    "+ str(value2) + "\n")
	newmodel.write("'PHASE FUNCTION   '    'isotropic'" + "\n")
	newmodel.write("'OMEGA PARTICLES  '    0.995"+ "\n")
	newmodel.write("'TAU_EXT PARTICLES'    999.000" + "\n")
	newmodel.write("\n")
	newmodel.write("'TYPE OF LAYER    '     'gas'" + "\n")
	newmodel.write("'P_TOP (bars)     '    1.400000" + "\n")
	newmodel.write("'P_BOTTOM (bars)  '    6" + "\n")
	newmodel.write("'K_CH4 (1/Km-amag)'    "+ str(value2) + "\n")
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
        noisy_y = []
        for item in y:
            value = np.random.uniform(-0.1,0.1) * item
            value22 = item + value
            noisy_y.append(value22)        
        return noisy_y

def func(x, T1, P2, T2, P3, value1, value2):
	values = []
        mu, mu0, deltaphi = read_scattering_angles()
	write_model(mu, mu0, deltaphi, T1, P2, T2, P3, value1, value2)
	execute_atmos()
        lon, iparf, muobs, mu0obs = read_obs()
        mod,londism = read_pointsdat()
        interpolated_y = interpolate_model(lon,londism,mod)
	interpolated_y_noise = insert_noise(interpolated_y) 
	return interpolated_y_noise

def optimize():
        nll = lambda *args: -lnlike(*args)
        result = op.minimize(nll, [T1_true, T2_true, P3_true, omega_true], args=(x, y, yerr))
        T1_ml, T2_ml, P3_ml, omega_ml = result["x"]
        return result["x"]

def calculate_diference(x, T1_true, P2_true, T2_true, P3_true, T1_ml, P2_ml, T2_ml, P3_ml, value1, value2, yerr):
        values = []
        y = func(x, T1_true, P2_true, T2_true, P3_true, value1, value2)
        model = func(x, T1_ml, P2_ml, T2_ml, P3_ml, value1, value2)
	for z in range(0,np.size(y)):
            item = ((y[z]-model[z])/(yerr*y[z]))**2
            values.append(item)
	error = -0.5*(np.sum(values))
	return error

def create_comparison_plot(x, T1_true, P2_true, T2_true, P3_true, T1_ml, P2_ml, T2_ml, P3_ml, value1, value2):
        Chi2_values = []
        errorbar = []
        y = func(x, T1_true, P2_true, T2_true, P3_true, value1, value2)
        model = func(x, T1_ml, P2_ml, T2_ml, P3_ml, value1, value2)
        for item in model:
            value = 0.05*item
            errorbar.append(value)
        for k in range(0,np.size(y)):
            item = ((y[k]-model[k])/(0.05*y[k]))**2
            Chi2_values.append(item)
            suma1 = np.sum(Chi2_values)
        #Plot the data:
        f = plt.figure()
        gs = gridspec.GridSpec(2,1 ,height_ratios=[10,4])
        ax1 = plt.subplot(gs[0])
        ax1.plot(x,y,'-',label = 'observations')
        #ax1.errorbar(x,model,yerr=errorbar)
        ax1.plot(x,model,'-',label = 'model')
        ax1.legend()
        #ax1.title('Error:'+ str(suma1))
        plt.ylabel('I/F')
        ax2 = plt.subplot(gs[1])
        ax2.plot(x,Chi2_values,'-',label = '$\chi^2$')
        ax2.legend()
        plt.xlabel('Longitude (SIII)')
        plt.ylabel('$\chi^2$')
        plt.show()
        f.savefig('output.png', format='png', dpi=1000)


# Define the probability function as likelihood  prior.
def lnprior(theta):
	T1, P2, T2, P3 = theta
	if (0.0 < T1 < 1.0 and 0.0015 < P2 < 0.45 and 5.0 < T2 < 20.0 and 0.015 < P3 < 0.45 and P2 < P3):
		return 0.0
	return -np.inf

def lnlike(theta, x, y, yerr):#inserta la definicon de yerr en algun lado!! He cambiado yerr por sigma!!!!!
        valuesChi2 = []
	T1, P2, T2, P3 = theta   
	interpolated_y_noise = func(x, T1, P2, T2, P3, value1, value2)    
	for k in range(0,np.size(y)):
            item = ((y[k]-interpolated_y_noise[k])/(0.05*y[k]))**2
            valuesChi2.append(item)
	suma = -0.5*(np.sum(valuesChi2))
        #record.write(str(T1) + "\t" + str(P2)  + "\t" + str(T2)  + "\t" + str(P3) + "\t" + str(value1) + "\t" + str(value2) + "\t" + str(suma) + "\n")
	return suma

def lnprob(theta, x, y, sigma):# Sigma
        lp = lnprior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + lnlike(theta, x, y, yerr)

def run_EMCEE(initial_position, dimensions, walkers, step_number, burn):
        ndim, nwalkers = dimensions, walkers
        pos = [initial_position + 1e-2*np.random.randn(ndim) for i in range(nwalkers)]
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))
        # Clear and run the production chain.
        print("Running MCMC...")
        sampler.run_mcmc(pos, step_number, rstate0=np.random.get_state())
        print("Done.")
        plt.clf()
        fig, axes = plt.subplots(4, 1, sharex=True, figsize=(8, 9))
        axes[0].plot(sampler.chain[:, :, 0].T, color="k", alpha=0.4)
        axes[0].yaxis.set_major_locator(MaxNLocator(5))
        axes[0].axhline(T1_ml, color="#888888", lw=2)
        axes[0].set_ylabel("$T1$")
        axes[1].plot(sampler.chain[:, :, 1].T, color="k", alpha=0.4)
        axes[1].yaxis.set_major_locator(MaxNLocator(5))
        axes[1].axhline(P2_ml, color="#888888", lw=2)
        axes[1].set_ylabel("$P2$")
        axes[2].plot(sampler.chain[:, :, 2].T, color="k", alpha=0.4)
        axes[2].yaxis.set_major_locator(MaxNLocator(5))
        axes[2].axhline(T2_ml, color="#888888", lw=2)
        axes[2].set_ylabel("$T2$")
        axes[3].plot(sampler.chain[:, :, 3].T, color="k", alpha=0.4)
        axes[3].yaxis.set_major_locator(MaxNLocator(5))
        axes[3].axhline(P3_ml, color="#888888", lw=2)
        axes[3].set_ylabel("$P3$")
        fig.tight_layout(h_pad=0.0)
        fig.savefig("line-time_4_variables_"+str(counter)+".png")
        # Make the triangle plot.
        burnin = burn
        samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))
        fig = triangle.corner(samples, labels=["$T1$","$P2$", "$T2$","$P3$"], truths=[T1_ml, P2_ml, T2_ml, P3_ml])
        fig.savefig("line-triangle_4_variables_"+str(counter)+".png")
######################################################################################################
# Hay que meter todo este bloque en un loop por el cual pasen valores de lambda y KCH4.
######################################################################################################
# Main block
if __name__ == '__main__':
   # Choose the "true" parameters.
   T1_true = 0.05
   P2_true = 0.03
   T2_true = 15.0
   P3_true = 0.1
   omega_true = 0.95
   
   lmbd = [0.451, 0.555, 0.727, 0.752, 0.889, 0.938]
   kch4 = [0.000, 0.000, 7.200, 0.027, 27.3, 0.050]
   
   for k in xrange(len(lmbd)):
        value1 = lmbd[k]
        value2 = kch4[k]
        print value1
        print value2
        counter += 1
        #record = open('record_'+str(counter)+'.txt','w')
        mu = []
        mu0 = []
        deltaphi = []
        x = []
        y = []
        yerr = 0.05
        lines11 = open('MT3-3.5N.dat','r').readlines()[1:]
        for line in lines11:
   	     parts = line.split()
	     x.append(float(parts[0]))
	     #y.append(float(parts[1]))
        y = func(x,0.05,15.0,0.1,0.95,value1,value2)
        nll = lambda *args: -lnlike(*args)
        result = op.minimize(nll, [T1_true, P2_true, T2_true, P3_true], args=(x, y, yerr))
        T1_ml, P2_ml, T2_ml, P3_ml = result["x"]
        calculate_diference(x, T1_true, P2_true, T2_true, P3_true, T1_ml, P2_ml, T2_ml, P3_ml, value1, value2, yerr)
        create_comparison_plot(x, T1_true, P2_true, T2_true, P3_true, T1_ml, P2_ml, T2_ml, P3_ml, value1, value2)
        # Problema con optimize
        #initial_position = [T1_true, T2_true, P3_true, omega_true]
        #run_EMCEE(initial_position, 4, 100, 700, 50)
        #record.close()
   
   
   
