# calculate tritium elastic cross section
# ref: Nuclear Physics A579 (1994) 596-626, note that there is an typo in eq 1 : -1/2->-1/4

#import pandas as pd
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys
import os
import glob
#import f1f209
from numpy.polynomial.polynomial import polyfit
from matplotlib.backends.backend_pdf import PdfPages
#import trit_elastic as doug
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
plt.rcParams['font.size'] = 15


def fq(qfm,R,Q):	#qfm=Q in fm-1

	gamma=0.8/np.sqrt(1.5) #fm

	temp=np.cos(qfm*R)+2*R**2/gamma**2*np.sin(qfm*R)/qfm/R
	temp=temp*Q/(1+2*R**2/gamma**2)
	temp=temp*np.exp(-0.25*qfm**2*gamma**2)

	return(temp.sum())


def fm2deg(eb, qfm):
	alpha  = 1/137.036
	hbarc2 = 0.38938 # GeV2*mbarn
	hbarc2 = hbarc2/1000.0*100 # GeV2*fm2: barn=100 fm2
	hbarc2 = 0.197327**2 # GeV2 fm2
	mass   = 2.80894327 #GeV, tritium
	
	q2 = qfm**2*hbarc2
	deg= np.arccos(1-1./(2*eb**2/q2-eb/mass) )
	deg=deg*180./np.pi
	return(deg)	

def h3_elastic(eb,theta):

#-----------parameters-------------------
	alpha  = 1/137.036
	hbarc2 = 0.38938 # GeV2*mbarn
	hbarc2 = hbarc2/1000.0*100 # GeV2*fm2: barn=100 fm2
	hbarc2 = 0.197327**2 # GeV2 fm2
	mu     = 2.9789248  #3H Nuclear Magnetic Moment 
	mass   = 2.80894327 #GeV, tritium
	a      = 3
	z      = 1
	
	
	# ref p611 table 1 SOG parameters for tritium
	Ri=np.array([0.1,0.5,0.9,1.3,1.6,2.0,2.4,2.9,3.4,4.0,4.6,5.2])
	Q_ch=np.array([0.054706,0.172505,0.313852,0.072056,0.225333,0.020849,0.097374,0.022273,0.011933,0.009121,0.0,0.0])
	Q_m=np.array([0.075234,0.164700,0.273033,0.037591,0.252089,0.027036,0.098445,0.040160,0.016696,0.015077,0.0,0.0])
	#-------------kinematics-----------------
	
	#eb =float(sys.argv[1])
	#theta  = float(sys.argv[2])
	
	deg = theta /180.*np.pi
	ep = eb/(1+eb/mass*(1-np.cos(deg)))
	q2 = 2*eb*ep*(1-np.cos(deg))
	tau = q2/4./mass**2
	
	qfm = np.sqrt(q2/hbarc2) # q in fm-1
	
	#-----------------------from Doug's MCEEP code---------------------------------------------
	#     Simple Calculation Of The Effective Q  
	#     Hesenberg and Blok, Ann. Rev. Nucl. Part. Sci. 33 (1983) 574.
	#----------------------------------------------------------------------
	#
	Qeff=qfm* ( 1. + 1.5*z*alpha*np.sqrt(hbarc2)/eb/(1.12*a**0.33333)  )
	
	mu=mu*a/z  # ???
	
	FM  = fq(Qeff,Ri,Q_m)
	FC  = fq(Qeff,Ri,Q_ch)
	
	
	mott = 4*z**2*alpha**2*hbarc2*ep**2/q2**2
	mott = mott *(np.cos(deg/2))**2 
	#mott2 = z**2*(hbarc2*alpha**2/4/eb**2)*np.cos(deg/2)**2/np.sin(deg/2)**4
	#print mott, mott2
	mott = mott * ep/eb
	C=FC**2/(1+tau)
	D=FM**2*tau*mu**2*(1/(1+tau)+2*np.tan(deg/2)**2)
	sigma_c=10*mott*C
	sigma_m=10*mott*D
	sigma=sigma_c+sigma_m#fm-2->mbarn
	#sigma=10*mott*(C+D)#fm-2->mbarn
	#print "me  :%e mbarn/sr" %sigma
#	sigma= doug.trit_elastic(eb*1000,deg)
#	print "doug: %e mbarn/sr" %sigma*10
	
	return (FC,FM,sigma_c,sigma_m,sigma,qfm)





def he3_elastic(eb,theta,newpar=0):
	#-----------parameters-------------------
	alpha  = 1/137.036
	hbarc2 = 0.38938 # GeV2*mbarn
	hbarc2 = hbarc2/1000.0*100 # GeV2*fm2: barn=100 fm2
	hbarc2 = 0.197327**2 # GeV2 fm2
	mu     = -2.12762485   #3He Nuclear Magnetic Moment 
	mass   = 2.808413  #GeV, He3
	a      = 3
	z      = 2


	# ref p611 table 1 SOG parameters for he3
	Ri   = np.array([0.1,0.5,0.9,1.3,1.6,2.0,2.4,2.9,3.4,4.0,4.6,5.2])
	Q_ch = np.array([0.027614,0.170847,0.219805,0.170486,0.134453,0.100953,0.074310,0.053970,0.023689,0.017502,0.002034,0.004338])
	Q_m  = np.array([0.059785,0.138368,0.281326,0.000037,0.289808,0.019056,0.114825,0.042296,0.028345,0.018312,0.007843,0.0])

	# new fitting from Scott Barcus
	if newpar==1:
		Q_m  = np.array([0.159649,0.0316168,0.277843,0.0364955,0.0329718,0.233469,0.117059,0.0581085,0.0485212,1.78E-12 ,0.0240927,8.95E-12])
		Ri   = np.array([0.3,0.7,0.9,1.1,1.5,1.6,2.2,2.7,3.3,4.2,4.3,4.8])
		Q_ch = np.array([0.0996392,0.214304,0.0199385,0.195676,0.0785533,0.167223,0.126926,0.0549379,0.0401401,0.0100803,0.0007217,4.99E-12])

	#-------------kinematics-----------------
	#	eb =float(sys.argv[1])
	#	theta  = float(sys.argv[2])

	deg = theta /180.*np.pi
	ep  = eb/(1+eb/mass*(1-np.cos(deg)))
	q2  = 2*eb*ep*(1-np.cos(deg))
	tau = q2/4./mass**2

	qfm = np.sqrt(q2/hbarc2) # q in fm-1

	#-----------------------from Doug's MCEEP code---------------------------------------------
	#     Simple Calculation Of The Effective Q  
	#     Hesenberg and Blok, Ann. Rev. Nucl. Part. Sci. 33 (1983) 574.
	#----------------------------------------------------------------------
	#
	Qeff=qfm * ( 1. + 1.5*z*alpha*np.sqrt(hbarc2)/eb/(1.12*a**0.33333)  )



	FM  = fq(Qeff,Ri,Q_m )
	FC  = fq(Qeff,Ri,Q_ch)


	mott = 4*z**2*alpha**2*hbarc2*ep**2/q2**2
	mott = mott *(np.cos(deg/2))**2 
	#mott2 = z**2*(hbarc2*alpha**2/4/eb**2)*np.cos(deg/2)**2/np.sin(deg/2)**4
	#print mott, mott2
	mott = mott * ep/eb
	C=FC**2/(1+tau)
	mu=mu*a/z  
	D=FM**2*tau*mu**2*(1/(1+tau)+2*np.tan(deg/2)**2)
	sigma_c=10*mott*C
	sigma_m=10*mott*D
	sigma=sigma_c+sigma_m#fm-2->mbarn
	#print "me :%e mbarn/sr" %sigma
#	sigma= doug.he3_elastic(eb*1000,deg)
#	print "doug: %e mbarn/sr" %(sigma*10)

	# print newpar, theta, FC, FM
	return (FC,FM,sigma_c,sigma_m,sigma,qfm)


#----------------------
#
#      main
#
#
#----------------------

plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['font.size'] = 12

# mass   = 2.808413  #GeV, He3

name  = sys.argv[1]
eb    = 2.216

if len(sys.argv)>2:
	eb = float(sys.argv[2])
if len(sys.argv)>3:
	newpar = float(sys.argv[3])

xmin  = 5
xmax  = 30
angle = np.linspace(xmin,xmax,400)

lfc    =  []
lfm    =  []
lsigm  =  []
lsigc  =  []
lsig   =  []
lqfm   =  []

nlfc   =  []
nlfm   =  []
nlsigm =  []
nlsigc =  []
nlsig  =  []
nlqfm  =  []


#-------get form factor and cross section from model--------

for deg in angle:
	if name=="he3":
		(fc,fm,sigc,sigm,sig,qfm)       = he3_elastic(eb,deg,0)
		(fcn,fmn,sigcn,sigmn,sign,qfmn) = he3_elastic(eb,deg,1)
		tt = "Helium3"
	elif name=="h3":
		(fc,fm,sigc,sigm,sig,qfm)=h3_elastic(eb,deg)
		tt="Tritium"

	lfc.append(fc)
	lfm.append(fm)
	lsigm.append(sigm)
	lsigc.append(sigc)
	lsig.append(sigm+sigc)
	lqfm.append(qfm**2)

	nlfc.append(fcn)
	nlfm.append(fmn)
	nlsigm.append(sigmn)
	nlsigc.append(sigcn)
	nlsig.append(sigmn+sigcn)
	nlqfm.append(qfmn**2)



lsig   = np.array(lsig)
lsigc  = np.array(lsigc)
lsigm  = np.array(lsigm)

nlsig  = np.array(nlsig)
nlsigc = np.array(nlsigc)
nlsigm = np.array(nlsigm)


#-------plot form factor (two fits) --------

ll=3 # linewidth
f1=plt.figure()
ax1 = f1.add_subplot(111)
ax1.axis([angle[0],angle[-1],1e-6,1])
# ax1.plot(angle,np.abs(lfc),'r--',linewidth=ll,label="$F_{ch}$")
# ax1.plot(angle,np.abs(lfm),'b:',linewidth=ll,label="$F_M$")
ax1.plot(angle,np.abs(nlfc), 'r-' , linewidth=ll,label="$F_{ch}$ Barcus")
ax1.plot(angle,np.abs(nlfm), 'b-' , linewidth=ll,label="$F_M$ Barcus"   )
ax1.plot(angle,np.abs(lfc),  'm--', linewidth=ll,label="$F_{ch}$ Amroun et al.")
ax1.plot(angle,np.abs(lfm),  'c--', linewidth=ll,label="$F_M$ Amroun et al."    )

ax1.legend(loc="upper right")
ax1.set_yscale('log')
plt.suptitle(tt+" Form Factor")
# ax1.set_xlabel('$fm^{-2}$')
ax1.set_xlabel('degree')
ax1.grid()

#-----second axis in fm-2------

ax2 = ax1.twiny()

ax1Ticks = ax1.get_xticks() 
temp=0
ax2value=[]

for deg in ax1Ticks:
 	(a,a,a,a,a,axq)=h3_elastic(eb,deg)
 	ax2value.append(axq**2)
ax2label=["%.1f" %z for z in ax2value]
# print ax2label
temp=['$fm^{-2}$']

ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(ax1Ticks)
ax2.set_xticklabels(ax2label)#,ha='right',minor=True)
ax2.set_xlabel(temp[0])
ax2.xaxis.set_label_coords(1.07,1.02)
ax1.xaxis.set_label_coords(1.07,-0.02)

# for la, x in zip(ax2label, ax1Ticks):    
#     ax1.annotate(la, xy=(x, 0.01), xycoords=('data', 'axes fraction'),
#         xytext=(0, -25), textcoords='offset points', va='top', ha='center')

# ax1.annotate(temp, xy=(80, 0.01), xycoords=('data', 'axes fraction'),
# 	xytext=(-22, -25), textcoords='offset points', va='top', ha='center')

# ax1.xaxis.set_label_coords(1.0,-0.0)


f1.set_size_inches(8,6)
fn = name+"-ff-%g.png" %eb
f1.savefig(fn)

print(fn+" created")
#print fn+" created"
plt.show()

del(f1)
#-------plot cross section (Amroun only) --------

f1=plt.figure()
ax1 = f1.add_subplot(111)
# ax1 = f1.add_axes((0.1,0.16,0.78,0.75))
# pos1 = ax1.get_position()
# print pos1
ax1.axis([angle[0],angle[-1],1e-16,0.01])

# ax1.axis([lqfm[1],30,1e-13,10])
#ax1.scatter(angle,lsigc,edgecolors='r',facecolors='r',marker='o',label="$\sigma(F_M=0)$")
#ax1.scatter(angle,lsigm,edgecolors='b',facecolors='b',marker='^',label="$\sigma(F_{ch}=0)$")
ax1.plot(lqfm,nlsig,'k-',  linewidth=ll,label="$\sigma_{tot}$")
ax1.plot(lqfm,nlsigc,'r--',linewidth=ll,label="$\sigma(F_M=0)$")
ax1.plot(lqfm,nlsigm,'b:', linewidth=ll,label="$\sigma(F_{ch}=0)$")
ax1.set_ylabel('mb/sr')#, color='b')
# ax1.tick_params('y', colors='b')
#ax1.annotate('Something', (0,0), (0, -20), xycoords='axes fraction', textcoords='offset points', va='top')
ax1.legend(loc="upper right")#,frameon=False)
ax1.set_yscale('log')
plt.suptitle(tt+" Elastic Cross Section with %g GeV Beam (Barcus's fits)" %eb)
ax1.set_xlabel("degree")
# ax1.set_xlabel("$fm^{-2}$")

ax1.grid()
# #=========second x axis=====================

ax2 = ax1.twiny()

ax1Ticks = ax1.get_xticks() 
temp=0
ax2value=[]

for deg in ax1Ticks:
 	(a,a,a,a,a,axq)=h3_elastic(eb,deg)
 	ax2value.append(axq**2)
ax2label=["%.1f" %z for z in ax2value]
# print ax2label
temp=['$fm^{-2}$']

ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(ax1Ticks)
ax2.set_xticklabels(ax2label)#,ha='right',minor=True)
ax2.set_xlabel(temp[0])
ax2.xaxis.set_label_coords(1.07,1.02)
ax1.xaxis.set_label_coords(1.07,-0.02)


f1.set_size_inches(8,6)
fn = name+"-sig-%g.png" %eb
f1.savefig(fn)
print(fn+" created")
#print fn+" created"
plt.show()

del(f1)

exit()


# --------------- save old code-------------

# # create second Axes. Note the 0.0 height

# #ax2 = f1.add_axes((0.1,0.1,0.8,0.0))
# ax2 = f1.add_axes((pos1.x0,pos1.y0-0.1,pos1.x1-pos1.x0,0.0))
# # print ax2.get_position()

ax1Ticks = ax1.get_xticks() 
#print ax1Ticks 
temp=0
ax2value=[]

# for deg in ax1Ticks:
#  	(a,a,a,a,a,axq)=h3_elastic(eb,deg)
#  	qfm.append(axq)
# ax2label=["%.2f" %z for z in qfm]
# # print ax2label
# temp=['$fm^{-1}$']

for qq in ax1Ticks[:-1]:
        dd = fm2deg(eb,np.sqrt(qq))
        ax2value.append(dd)
 	#print eb, eb/(1+eb/mass*(1-np.cos(dd/180.0*np.pi))), dd, qq 
        ax2label=["%.2f" %z for z in ax2value]
        temp='$degree$'

# print ax2label
# ax2.set_xbound(ax1.get_xbound())
# ax2.set_xticks(ax1Ticks)
# ax2Ticks = ax2.get_xticks() 
# print ax2Ticks
# ax2.yaxis.set_visible(False) # hide the yaxis
# # ax2.set_xbound(ax1.get_xbound())
# #ticks = x_min, x_max = ax1.get_xlim()
# #ticks = [(tick - x_min)/(x_max - x_min) for tick in ax1.get_yticks()]
# ax2.set_xticklabels(ax2label)
# ax2.set_xlabel(r"fm2")
# temp=['$fm^{-2}$']
#for la, x in zip(ax2label+temp, np.append(ax1Ticks, 42.5)):
for la, x in zip(ax2label, ax1Ticks):    
    ax1.annotate(la, xy=(x, 0.01), xycoords=('data', 'axes fraction'),
        xytext=(0, -25), textcoords='offset points', va='top', ha='center')

ax1.annotate(temp, xy=(80, 0.01), xycoords=('data', 'axes fraction'),
	xytext=(-22, -25), textcoords='offset points', va='top', ha='center')

ax1.xaxis.set_label_coords(1.08,-0.0)



#===============second y axis=========
# ax3 = ax1.twinx()
# ax3.plot(lqfm,lsigc/lsig*100,'g-.',linewidth=ll-0.5,label="$\sigma(F_M=0)/\sigma_{tot}}$")
# ax3.set_ylabel("percent", color='g')
# ax3.tick_params('y', colors='g')
# ax3.axis([lqfm[1],lqfm[-1],0,120])
# ax3.legend(loc="upper right")

# f1.tight_layout()
