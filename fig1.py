import numpy as np
import math
from pylab import *
import matplotlib.patches as patches
from scipy.optimize import minimize_scalar

import nspy as nspy


## Plot
fig = figure(figsize=(6,4), dpi=80)
rc('font', family='serif')
rc('xtick', labelsize='xx-small')
rc('ytick', labelsize='xx-small')

gs = GridSpec(1, 1)
#gs.update(wspace = 0.34)
#gs.update(hspace = 0.4)

ax = subplot(gs[0,0])
ax.axis('off')
ax.set_aspect('equal')

#image dimensions
mlim= 1.3
xmin = -1.5
xmax = mlim
ymin = -0.5
ymax = mlim

xlen = xmax-xmin
ylen = ymax-ymin
aspr=ylen/xlen

ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)


#variables
#incl = pi/8
nspy.incl = pi/2.5
nspy.req = 1
nspy.u = 0.0
nspy.rg = nspy.u * nspy.req
nspy.muc = -nspy.rg / (nspy.req - nspy.rg)
nspy.o21 = 0.3


#coordinates of spot and observer
spot_theta=pi/5
spot_phi=0.67

#obs_theta=incl
obs_theta=pi/2.7
obs_phi=-pi/1.5

#default front and back line styles
fmt={'color':'k','linestyle':'solid',}
fmtb = {'color':'k','linestyle':'dotted',}



#-----------------------------------------------------
#draw xyz axis
ax = nspy.draw_axis(ax, obs_phi,      pi/2, rfac=1.4)
ax = nspy.draw_axis(ax, obs_phi+pi/2, pi/2, rfac=1.4)
ax = nspy.draw_axis(ax, 0.0,          0.0 , rfac=1.8)

#rotation direction, i.e. curly arrow
phistart=-pi/2-0.8
phistop=pi/2+0.5
wsize=0.15
wheight=1.2

ax=nspy.draw_latitude(ax, wsize, fmt=fmt, start=phistart, stop=phistop, rfac=wheight)
yy1=nspy.y(phistop,wsize)*wheight
zz1=nspy.z(phistop,wsize)*wheight

yy2 = -0.01
zz2 = 0.002

fmta= {'color':'black', 'linestyle':'solid', 'lw':1.0, 'head_width': 0.02, 'head_length': 0.04,}
ax.add_patch(patches.FancyArrow(
            yy1, zz1,
            yy2, zz2,
            **fmta
            ))



#-----------------------------------------------------
#borders
fmtO={'color':'red','linestyle':'solid',}
#ax = draw_outline(ax, fmtO)

ax=nspy.draw_longitude(ax, pi/2, fmt=fmt)
ax=nspy.draw_longitude(ax, -pi/2, fmt=fmt)

#ax=nspy.draw_longitude(ax, -0.1, fmt=fmt)
#ax=nspy.draw_longitude(ax, -pi-0.1, fmt=fmt)
#ax=nspy.draw_longitude(ax, -pi-0.1, fmt=fmtb, backside=True)

ax=nspy.draw_latitude(ax, pi/2, fmt=fmt)
ax=nspy.draw_latitude(ax, pi/2, backside=True, fmt=fmtb)


#fmt2 = {'color':'k','linestyle':'solid','lw':0.8,}
#fmt2b = {'color':'k','linestyle':'dotted','lw':0.8,}
#for thetaa in np.linspace(0.0, pi/2, 10):
#    ax = nspy.draw_latitude(ax, thetaa, fmt=fmt2)
#    ax = nspy.draw_latitude(ax, thetaa, fmt=fmt2b, backside=True)
#                       
#for phii in np.linspace(0.0, 2*pi, 30):
#    ax = nspy.draw_longitude(ax, phii, fmt=fmt2,
#                        start=0.0,
#                        stop=pi/2
#                        )


#-----------------------------------------------------
#spot position
ax=nspy.draw_longitude(ax, spot_phi, fmt=fmt)
ax=nspy.draw_radial(ax, spot_phi, spot_theta, fmt)
#ax = draw_axis(ax, spot_phi, spot_theta, rfac=1.2)

fmta= {'color':'black', 'linestyle':'solid', 'lw':1.4, 'head_width': 0.02, 'head_length': 0.05,}
alen = 0.6
ax.add_patch(patches.FancyArrow(
            nspy.y(spot_phi, spot_theta)*1.0,  nspy.z(spot_phi, spot_theta)*1.0,
            nspy.y(spot_phi, spot_theta)*alen, nspy.z(spot_phi, spot_theta)*alen,
            **fmta
            ))


fmta= {'color':'black', 'linestyle':'solid', 'lw':1.4, 'head_width': 0.02, 'head_length': 0.05,}
ax = nspy.draw_normal(ax, spot_phi, spot_theta, fmt=fmta, alen=0.4)

#for spot_theta in np.linspace(0.1, pi/2, 10):
#    ax = draw_normal(ax, spot_phi, spot_theta, alen=0.1)


ax=nspy.draw_longitude(ax, spot_phi, fmt=fmt, start=0, stop=spot_theta, rfac=0.17)
ax=nspy.draw_longitude(ax, spot_phi, fmt=fmt, start=0, stop=spot_theta, rfac=0.17, backside=True) #theta angle
#add(p, PlotLabel(.53, .58, "<i>\\theta</i>"))

#help axis
ax = nspy.draw_radial(ax, spot_phi, pi/2, fmt)


#circular spot
rho = 0.1
ax = nspy.draw_spot(ax, spot_phi, spot_theta, rho)

#ax=nspy.draw_longitude(ax, spot_phi-rho, fmt=fmt)
#ax=nspy.draw_longitude(ax, spot_phi+rho, fmt=fmt)
#ax=nspy.draw_radial(ax, spot_phi, spot_theta+rho, fmt)
#ax=nspy.draw_radial(ax, spot_phi, spot_theta-rho, fmt)


#-----------------------------------------------------
#observer
fmti = {'color':'k','linestyle':'dashed',}

ax=nspy.draw_radial(ax, obs_phi, obs_theta, fmt=fmti, rfac=1.65)
ax=nspy.draw_radial(ax, obs_phi, pi/2, fmt=fmti)

ax=nspy.draw_longitude(ax, obs_phi, fmt=fmt, start=0, stop=obs_theta, rfac=0.1) #incl angle
ax=nspy.draw_longitude(ax, obs_phi, fmt=fmt, start=0, stop=obs_theta, rfac=0.1, backside=True) #incl angle


ax=nspy.draw_latitude(ax, pi/2, fmt=fmt, start=spot_phi, stop=obs_phi, rfac=0.15) #phi angle
ax=nspy.draw_latitude(ax, pi/2, fmt=fmt, start=spot_phi, stop=obs_phi, rfac=0.15, backside=True) #phi angle

obspx = nspy.y(obs_phi, obs_theta)*1.8
obspy = nspy.z(obs_phi, obs_theta)*1.8
#ax = nspy.draw_observer(ax, -1.5, 0.8, angle=-0.45, sangle=0.4, alen=0.15)
ax = nspy.draw_observer(ax, obspx, obspy, angle=obs_theta-pi/2-0.16, sangle=0.4, alen=0.15)


#observer axis
fmt_obs= {'color':'k', 
        'linestyle':'solid', 
        'lw':1.5, 
        'head_width': 0.02, 
        'head_length': 0.04,}

obspx = nspy.y(obs_phi, obs_theta)*1.25
obspy = nspy.z(obs_phi, obs_theta)*1.25
ax = nspy.draw_axis(ax, obs_phi,      obs_theta,      rfac=0.2, startp=(obspx, obspy),fmt=fmt_obs)
ax = nspy.draw_axis(ax, obs_phi+pi/2, pi/2,           rfac=0.2, startp=(obspx, obspy),fmt=fmt_obs)
ax = nspy.draw_axis(ax, obs_phi,      obs_theta-pi/2, rfac=0.2, startp=(obspx, obspy),fmt=fmt_obs)


#texts
#-----------------------------------------------------
lsize = 15.0
ax.text(0.0, -0.15, "$\\phi$", va='center', ha='center', size=lsize)
ax.text(0.06, 0.2, "$\\theta$", va='center', ha='center', size=lsize)
ax.text(-0.05, 0.15, "$i$", va='center', ha='center', size=lsize)

ax.text(-0.6, -0.4, "$x$", va='center', ha='center', size=lsize)
ax.text(0.06, 1.1, "$y$", va='center', ha='center', size=lsize)
ax.text(-1.3, 0.1, "$z$", va='center', ha='center', size=lsize)
#ax.text(-0.1, 0.9, "$\\hat{\\Omega}$", va='center', ha='center', size=lsize)

ax.text(0.55, 0.96, "$\\hat{r}$", va='center', ha='center', size=lsize)
ax.text(0.38, 0.92, "$\\hat{n}$", va='center', ha='center', size=lsize)

ax.text(-1.15, 0.5, "$\\hat{x}$", va='center', ha='center', size=lsize)
ax.text(-0.83, 0.88,  "$\\hat{y}$", va='center', ha='center', size=lsize)
#ax.text(-1.2, 0.8,  "$\\hat{z}$", va='center', ha='center', size=lsize)

savefig('fig1.png', bbox_inches='tight')
