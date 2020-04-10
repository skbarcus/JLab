import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.widgets import Slider, Button, RadioButtons
import numpy as np

fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.30)
plt.ylim(-75,10)
plt.xlim(-150,30)
#plt.subplots_adjust(left=0.25, bottom=0.25)
#fig.ylim(0,10)

#Length units are cm.
f0 = 5
delta_f = 0.1
x1 = 1
trgt_len = 40
trgt_wid = 2.1
dend = 1
dSHMS = 116.77
dAcc = 16.14
theta_init = 13*np.pi/180

#Starting values for collimator positions (x1,y1).
phi_init = 0*np.pi/180
Dx1_init = -50
Dy1_init = -20
Ux1_init = -10
Uy1_init = -20

#Define the length and width of the two collimators.
dx = 3.23*2.54
dy = 0.875*2.54
ux = 3.23*2.54
uy = 1.86*2.54

#Draw the target
trgt_top, = plt.plot([-trgt_len/2,trgt_len/2],[trgt_wid/2,trgt_wid/2],'r')
trgt_top, = plt.plot([-trgt_len/2,trgt_len/2],[-trgt_wid/2,-trgt_wid/2],'r')

#Draw the SHMS acceptance.
SHMS, = plt.plot([-dSHMS*np.cos(theta_init)-0.5*dAcc*np.sin(theta_init),-dSHMS*np.cos(theta_init)+0.5*dAcc*np.sin(theta_init)],[-dSHMS*np.sin(theta_init)+0.5*dAcc*np.cos(theta_init),-dSHMS*np.sin(theta_init)-0.5*dAcc*np.cos(theta_init)])

#Define line from beam right of the downstream endcap to the lower angle side of the SHMS.
line1, = plt.plot([-trgt_len/2+dend,-dSHMS*np.cos(theta_init)-0.5*dAcc*np.sin(theta_init)],[trgt_wid/2,-dSHMS*np.sin(theta_init)+0.5*dAcc*np.cos(theta_init)],'b--')

#Define line from beam left of the downstream endcap to the higher angle side of the SHMS.
line2, = plt.plot([-trgt_len/2+dend,-dSHMS*np.cos(theta_init)+0.5*dAcc*np.sin(theta_init)],[-trgt_wid/2,-dSHMS*np.sin(theta_init)-0.5*dAcc*np.cos(theta_init)],'b--')

#Define line from beam right of the upstream endcap to the lower angle side of the SHMS.
line3, = plt.plot([trgt_len/2-dend,-dSHMS*np.cos(theta_init)-0.5*dAcc*np.sin(theta_init)],[trgt_wid/2,-dSHMS*np.sin(theta_init)+0.5*dAcc*np.cos(theta_init)],'g--')

#Define line from beam left of the upstream endcap to the higher angle side of the SHMS.
line4, = plt.plot([trgt_len/2-dend,-dSHMS*np.cos(theta_init)+0.5*dAcc*np.sin(theta_init)],[-trgt_wid/2,-dSHMS*np.sin(theta_init)-0.5*dAcc*np.cos(theta_init)],'g--')

#Draw the four lines making up the downstream collimator. 
#Top left is (x1,y1), top right is (x2,y2), bottom left is (x3,y3), bottom right is (x4,y4).
#Line from (x1,y1) to (x2,y2)
ld1, = plt.plot([Dx1_init,Dx1_init + dx*np.cos(phi_init)],[Dy1_init,Dy1_init + dx*np.sin(phi_init)],'r')
#Line from (x1,y1) to (x3,y3)
ld2, = plt.plot([Dx1_init,Dx1_init + dy*np.sin(phi_init)],[Dy1_init,Dy1_init - dy*np.cos(phi_init)],'r')
#Line from (x2,y2) to (x4,y4)
ld3, = plt.plot([Dx1_init + dx*np.cos(phi_init),Dx1_init + dy*np.sin(phi_init) + dx*np.cos(phi_init)],[Dy1_init + dx*np.sin(phi_init),Dy1_init - dy*np.cos(phi_init) + dx*np.sin(phi_init)],'r')
#Line from (x3,y3) to (x4,y4)
ld4, = plt.plot([Dx1_init + dy*np.sin(phi_init),Dx1_init + dy*np.sin(phi_init) + dx*np.cos(phi_init)],[Dy1_init - dy*np.cos(phi_init),Dy1_init - dy*np.cos(phi_init) + dx*np.sin(phi_init)],'r')

#Draw the four lines making up the upstream collimator.
#Top left is (x1,y1), top right is (x2,y2), bottom left is (x3,y3), bottom right is (x4,y4).
#Line from (x1,y1) to (x2,y2)
ld5, = plt.plot([Ux1_init,Ux1_init + ux*np.cos(phi_init)],[Uy1_init,Uy1_init + ux*np.sin(phi_init)],'r')
#Line from (x1,y1) to (x3,y3)
ld6, = plt.plot([Ux1_init,Ux1_init + uy*np.sin(phi_init)],[Uy1_init,Uy1_init - uy*np.cos(phi_init)],'r')
#Line from (x2,y2) to (x4,y4)
ld7, = plt.plot([Ux1_init + ux*np.cos(phi_init),Ux1_init + uy*np.sin(phi_init) + ux*np.cos(phi_init)],[Uy1_init + ux*np.sin(phi_init),Uy1_init - uy*np.cos(phi_init) + ux*np.sin(phi_init)],'r')
#Line from (x3,y3) to (x4,y4)
ld8, = plt.plot([Ux1_init + uy*np.sin(phi_init),Ux1_init + uy*np.sin(phi_init) + ux*np.cos(phi_init)],[Uy1_init - uy*np.cos(phi_init),Uy1_init - uy*np.cos(phi_init) + ux*np.sin(phi_init)],'r')

axcolor = 'lightgoldenrodyellow'
#Downstream Collimator phi slider.
axDphi = plt.axes([0.1, 0.05, 0.3, 0.03], facecolor=axcolor)
sDphi = Slider(axDphi, 'Dphi', 0., 180.0, valinit=phi_init, valstep=delta_f)

#Upstream Collimator phi Slider.
axUphi = plt.axes([0.6, 0.05, 0.3, 0.03], facecolor=axcolor)
sUphi = Slider(axUphi, 'Uphi', 0., 180.0, valinit=phi_init, valstep=delta_f)

#Downstream Collimator X slider.
axDx = plt.axes([0.1, 0.15, 0.3, 0.03], facecolor=axcolor)
sDx = Slider(axDx, 'Dx', -200., 30.0, valinit=Dx1_init, valstep=delta_f)

#Downstream Collimator Y slider.
axDy = plt.axes([0.1, 0.1, 0.3, 0.03], facecolor=axcolor)
sDy = Slider(axDy, 'Dy', -75., 10.0, valinit=Dy1_init, valstep=delta_f)

#Upstream Collimator X slider.
axUx = plt.axes([0.6, 0.15, 0.3, 0.03], facecolor=axcolor)
sUx = Slider(axUx, 'Ux', -200., 30.0, valinit=Ux1_init, valstep=delta_f)

#Upstream Collimator Y slider.
axUy = plt.axes([0.6, 0.1, 0.3, 0.03], facecolor=axcolor)
sUy = Slider(axUy, 'Uy', -75., 10.0, valinit=Uy1_init, valstep=delta_f)

#Theta slider
axtheta = plt.axes([0.1, 0.2, 0.8, 0.03], facecolor=axcolor)
stheta = Slider(axtheta, 'Theta', 11.0, 21.0, valinit=13.0, valstep=delta_f)

def update(val):
    
    theta = stheta.val*np.pi/180
    Dphi = sDphi.val*np.pi/180
    Uphi = sUphi.val*np.pi/180
    Dx1 = sDx.val
    Dy1 = sDy.val
    Ux1 = sUx.val
    Uy1 = sUy.val

    #Define the new coordinates for the SHMS face based on the slider theta value.
    SHMS.set_xdata([-dSHMS*np.cos(theta)-0.5*dAcc*np.sin(theta),-dSHMS*np.cos(theta)+0.5*dAcc*np.sin(theta)])
    SHMS.set_ydata([-dSHMS*np.sin(theta)+0.5*dAcc*np.cos(theta),-dSHMS*np.sin(theta)-0.5*dAcc*np.cos(theta)])

    #Keep line1 attached to the lower angle side of the SHMS.
    line1.set_xdata([-trgt_len/2+dend,-dSHMS*np.cos(theta)-0.5*dAcc*np.sin(theta)])
    line1.set_ydata([trgt_wid/2,-dSHMS*np.sin(theta)+0.5*dAcc*np.cos(theta)])

    #Keep line2 attached to the higher angle side of the SHMS.
    line2.set_xdata([-trgt_len/2+dend,-dSHMS*np.cos(theta)+0.5*dAcc*np.sin(theta)])
    line2.set_ydata([-trgt_wid/2,-dSHMS*np.sin(theta)-0.5*dAcc*np.cos(theta)])

    #Keep line3 attached to the lower angle side of the SHMS.
    line3.set_xdata([trgt_len/2-dend,-dSHMS*np.cos(theta)-0.5*dAcc*np.sin(theta)])
    line3.set_ydata([trgt_wid/2,-dSHMS*np.sin(theta)+0.5*dAcc*np.cos(theta)])

    #Keep line4 attached to the higher angle side of the SHMS.
    line4.set_xdata([trgt_len/2-dend,-dSHMS*np.cos(theta)+0.5*dAcc*np.sin(theta)])
    line4.set_ydata([-trgt_wid/2,-dSHMS*np.sin(theta)-0.5*dAcc*np.cos(theta)])

    #Update the movement of the downstream collimator with the sliders.
    ld1.set_xdata([Dx1,Dx1 + dx*np.cos(Dphi)])
    ld1.set_ydata([Dy1,Dy1 + dx*np.sin(Dphi)])
    ld2.set_xdata([Dx1,Dx1 + dy*np.sin(Dphi)])
    ld2.set_ydata([Dy1,Dy1 - dy*np.cos(Dphi)])
    ld3.set_xdata([Dx1 + dx*np.cos(Dphi),Dx1 + dy*np.sin(Dphi) + dx*np.cos(Dphi)])
    ld3.set_ydata([Dy1 + dx*np.sin(Dphi),Dy1 - dy*np.cos(Dphi) + dx*np.sin(Dphi)])
    ld4.set_xdata([Dx1 + dy*np.sin(Dphi),Dx1 + dy*np.sin(Dphi) + dx*np.cos(Dphi)])
    ld4.set_ydata([Dy1 - dy*np.cos(Dphi),Dy1 - dy*np.cos(Dphi) + dx*np.sin(Dphi)])

    #Update the movement of the upstream collimator with the sliders.
    ld5.set_xdata([Ux1,Ux1 + ux*np.cos(Uphi)])
    ld5.set_ydata([Uy1,Uy1 + ux*np.sin(Uphi)])
    ld6.set_xdata([Ux1,Ux1 + uy*np.sin(Uphi)])
    ld6.set_ydata([Uy1,Uy1 - uy*np.cos(Uphi)])
    ld7.set_xdata([Ux1 + ux*np.cos(Uphi),Ux1 + uy*np.sin(Uphi) + ux*np.cos(Uphi)])
    ld7.set_ydata([Uy1 + ux*np.sin(Uphi),Uy1 - uy*np.cos(Uphi) + ux*np.sin(Uphi)])
    ld8.set_xdata([Ux1 + uy*np.sin(Uphi),Ux1 + uy*np.sin(Uphi) + ux*np.cos(Uphi)])
    ld8.set_ydata([Uy1 - uy*np.cos(Uphi),Uy1 - uy*np.cos(Uphi) + ux*np.sin(Uphi)])
    
    fig.canvas.draw_idle()

sDphi.on_changed(update)
sUphi.on_changed(update)
sDx.on_changed(update)
sDy.on_changed(update)
sUx.on_changed(update)
sUy.on_changed(update)
stheta.on_changed(update)

plt.show()
