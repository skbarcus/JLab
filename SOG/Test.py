
import numpy as np
#matplotlib inline
import matplotlib.pyplot as plt

f = open('Test.txt', 'r')
data = np.genfromtxt(f, delimiter='   ')
#delete(data,0,0) # Erases the first row (i.e. the header)

X = data[:,0]
Y = data[:,1]

plt.plot(data[:,0],data[:,1],'o')

print data[:,0] 
print data[:,1]
print X
print X[0]
print Y

f.close()
