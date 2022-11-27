import numpy as np
from sklearn.preprocessing import normalize
import random

x = np.arange(16.0).reshape(4, 4)

print('x', x)

x1 = np.hsplit(x,np.array([3,6]))

print('x1',x1)

x2 = np.hsplit(x,np.array([1,2]))

print('x2',x2)

x3 = np.hsplit(x,np.array([1,3]))

print('x3',x3)


x4 = np.hsplit(x,np.array([1,4]))

print('x4',x4)

x5 = np.hsplit(x,np.array([1]))

#x5 = np.array(x5)
#print('x5.shape()',x5.shape)
print('x5',x5)
x5_0 = np.array(x5[0]) 
x5_1 = np.array(x5[1]) 
print('x5_0.shape',x5_0.shape)
print('x5_0',x5_0)
print('x5_1.shape',x5_1.shape)
print('x5_1',x5_1)

x6 = np.delete(x,2,0)
print('x6.shape',x6.shape)
print('x6',x6)

x7 = normalize(x, axis=0)
print('x7',x7)

test = [0,0,0,0,0,0]
print(len(test))

print(random.uniform(0, 1))
print(random.uniform(0, 1))
print(random.uniform(0, 1))
print(random.uniform(0, 1))
print(random.uniform(0, 1))
