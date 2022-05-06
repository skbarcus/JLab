import numpy as np
import matplotlib.pyplot as plt

pi = 3.141592654
deg2rad = pi/180.0
rad2deg = 180.0/pi

x = (1,2,3,4,5)

y = np.power(x,2)

fig, ax = plt.subplots(figsize=(12,6))

plt.scatter(x,y,edgecolors='r',facecolors='none',s=80)

#Now plot same points inside first points but filled and with variable size.
#def size(x):
#    val = 80/x
#    return val

size = [80/(2*n) for n in x]

plt.scatter(x,y,color='r',s=size)

plt.show()

