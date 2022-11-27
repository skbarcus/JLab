from matplotlib import pyplot as plt
import numpy as np
import math
import time
from matplotlib.transforms import Transform
from matplotlib.ticker import(AutoLocator, AutoMinorLocator)
from collections import OrderedDict

start = time.time() #Start a timer.

#X direction (vertical) uniformity by row.
#SBS 4
xuni_x4 = np.linspace(1,24,24)
xuni_y4 = (0.0476,0.0566,0.0487,0.0503,0.0521,0.051,0.0498,0.0499,0.0376,0.045,0.0471,0.0454,0.0493,0.0477,0.0424,0.0382,0.0459,0.0442,0.0486,0.0478,0.0481,0.0409,0.0371,0.0283)
xuni_y4_uncer = (0.0382,0.03,0.0326,0.0275,0.0262,0.0232,0.0299,0.0269,0.0311,0.025,0.0252,0.0297,0.0257,0.0298,0.0287,0.0326,0.0241,0.0293,0.0261,0.0283,0.0237,0.0318,0.0314,0.0348)
#SBS 8
xuni_x8 = (4,6,10,11,13,15,17,19,21,22,23)
xuni_y8 = (0.0297,0.0254,0.0195,0.0226,0.0325,0.0375,0.0335,0.0393,0.0491,0.0197,0.0385)
xuni_y8_uncer = (0.0269,0.0307,0.0331,0.0373,0.0316,0.0265,0.0332,0.0238,0.0265,0.0442,0.0294)
#xuni_y8 = (0,0,0,0.0297,0,0.0254,0,0,0,0.0195,0.0226,0,0.0325,0,0.0375,0,0.0335,0,0.0393,0,0.0491,0.0197,0.0385,0)
#xuni_y8_uncer = (0,0,0,0.0269,0,0.0307,0,0,0,0.0331,0.0373,0,0.0316,0,0.0265,0,0.0332,0,0.0238,0,0.0265,0.0442,0.0294,0)

fig, ax = plt.subplots(figsize=(10,10))
xsbs4 = plt.errorbar(xuni_x4,xuni_y4,yerr=xuni_y4_uncer,color='blue',label="SBS-4",marker='o',linestyle='',capsize=6.0,markersize=8.0)

xsbs8 = plt.errorbar(xuni_x8,xuni_y8,yerr=xuni_y8_uncer,color='red',label="SBS-8",marker='o',linestyle='',capsize=6.0,markersize=8.0)

#print(xsbs8)
#xsbs8.remove()
#xuni_y8[0].remove()

plt.title("HCal X-Direction (Vertical) Uniformity",fontsize=36)
plt.xlabel("Row Number",fontsize=28)
plt.xticks(fontsize=24,ticks=np.linspace(1,24,24))
plt.ylabel("Detection Efficieny",fontsize=28)
plt.yticks(fontsize=24)
plt.legend(loc="upper left",fontsize=15)
#ax.xaxis.set_major_locator(plt.MaxNLocator(24))

plt.show()

#Y direction (vertical) uniformity by row.
#SBS 4
yuni_x4 = (3,4,5,6,7,8,9,10,11,12)
yuni_y4 = (0.0727,0.0724,0.0699,0.0664,0.0643,0.0678,0.0635,0.0589,0.0572,0.0577)
yuni_y4_uncer = (0.0437,0.0351,0.038,0.0294,0.0256,0.0277,0.0279,0.0356,0.0296,0.0284)

#SBS 8
yuni_x8 = (8,9,11)
yuni_y8 = (0.0416,0.0374,0.0256)
yuni_y8_uncer = (0.0785,0.0270,0.0444)

fig, ax = plt.subplots(figsize=(10,10))
plt.errorbar(yuni_x4,yuni_y4,yerr=yuni_y4_uncer,color='blue',label="SBS-4",marker='o',linestyle='',capsize=6.0,markersize=8.0)
plt.errorbar(yuni_x8,yuni_y8,yerr=yuni_y8_uncer,color='red',label="SBS-8",marker='o',linestyle='',capsize=6.0,markersize=8.0)

plt.title("HCal Y-Direction (Horizontal) Uniformity",fontsize=36)
plt.xlabel("Column Number",fontsize=28)
plt.xticks(fontsize=24,ticks=np.linspace(1,12,12))
plt.ylabel("Detection Efficieny",fontsize=28)
plt.yticks(fontsize=24)
#ax.xaxis.set_major_locator(plt.MaxNLocator(24))

plt.legend(loc="upper left",fontsize=15)

plt.show()

print("The script took %.2f seconds (%.2f minutes or %.2f hours) to run." % (time.time() - start, (time.time() - start)/60.,(time.time() - start)/60./60.)) #Print time to run.
