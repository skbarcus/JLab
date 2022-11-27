import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import stats

start_amp = 5           #He3 = 3.140. H3 = 3.876.Starting guess for Gaus histo fit amplitude. 
start_avg = 1.91        #He3 = 1.911. H3 = 1.927. Starting guess for Gaus histo fit average. 
start_std_dev = 0.05    #He3 = 0.01903. H3 = 0.1418. Starting guess for Gaus histo fit standard deviation. 

nbins = 16              #He3 = 16. H3 = 10.Number of bins used in the fitted histogram.

x_data = np.loadtxt('bootstrapped_charge_radii_he3.txt',skiprows=0)
print(x_data)

#def gaus_fit(x, C, mu, sigma):
#    return ( C * np.exp(-1.0 * (x - mu)**2 / (2 * sigma**2)) )

# Creating histogram
#fig, ax = plt.subplots(figsize =(10, 7))
#nbins = 20
#ax.hist(x_data, bins = nbins)

#popt, pcov = curve_fit(gaus_fit, xdata=binscenters, ydata=data_entries_1, p0=[start_amp, start_avg, start_std_dev])
#print('3He popt =',popt)
 
# Show plot
#plt.show()

hist, bin_edges = np.histogram(x_data,bins=nbins)
print('hist =',hist)
print('bin_edges =',bin_edges)
#hist=hist/sum(hist)  #Use if normalizing Gaussian to unity.

n = len(hist)
x_hist=np.zeros((n),dtype=float) 
for ii in range(n):
    x_hist[ii]=(bin_edges[ii+1]+bin_edges[ii])/2
    
y_hist=hist
       
print('x_hist =',x_hist)
print('y_hist =',y_hist)

#Calculating the Gaussian PDF values given Gaussian parameters and random variable X
def gaus(X,C,X_mean,sigma):
    return C*np.exp(-(X-X_mean)**2/(2*sigma**2))

mean = sum(x_hist*y_hist)/sum(y_hist)                  
sigma = sum(y_hist*(x_hist-mean)**2)/sum(y_hist) 

#Gaussian least-square fitting process
param_optimised,param_covariance_matrix = curve_fit(gaus,x_hist,y_hist,p0=[max(y_hist),mean,sigma],maxfev=5000)

#print fit Gaussian parameters
print("Fit parameters: ")
print("=====================================================")
print("C = ", param_optimised[0], "+-",np.sqrt(param_covariance_matrix[0,0]))
print("X_mean =", param_optimised[1], "+-",np.sqrt(param_covariance_matrix[1,1]))
print("sigma = ", param_optimised[2], "+-",np.sqrt(param_covariance_matrix[2,2]))
print("\n")

#STEP 4: PLOTTING THE GAUSSIAN CURVE -----------------------------------------
fig = plt.figure()
x_hist_2=np.linspace(np.min(x_hist),np.max(x_hist),500)
plt.plot(x_hist_2,gaus(x_hist_2,*param_optimised),'r.:',label='Gaussian fit')
plt.legend(fontsize=20)

#Normalise the histogram values
#weights = np.ones_like(x_data) / len(x_data) #Use if normalizing Gaussian to unity.
#plt.hist(x_data, weights=weights) #Use if normalizing Gaussian to unity.

plt.hist(x_data,bins=nbins)
#setting the label,title and grid of the plot
plt.title('Bootstrapped World Data Charge Radii for $^3$He',fontsize=24)
plt.xlabel("Charge Radius (fm)",fontsize=20)
plt.ylabel("Occurrences",fontsize=20)
plt.grid("on")
plt.show()
