import numpy as np
import math
from pylab import *
import matplotlib.patches as patches
from scipy.optimize import minimize_scalar

import nspy as nspy

from matplotlib import cm

import palettable as pal



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
#mlim= 1.3
xmin = -3.6
xmax =  3.6
ymin = -1.8
ymax =  2.2

xlen = xmax-xmin
ylen = ymax-ymin
aspr=ylen/xlen

ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)


#variables
nspy.incl = pi/3.0 #default just touch E bottom
#nspy.incl = pi/2.5
#nspy.incl = pi/2.6
nspy.req = 1
nspy.u = 0.01
nspy.rg = nspy.u * nspy.req
nspy.muc = -nspy.rg / (nspy.req - nspy.rg)
nspy.o21 = 0.70


#coordinates of spot and observer
spot_theta=pi/5
spot_phi=0.67

#obs_theta=incl
obs_theta=pi/2.7
obs_phi=-pi/1.5

#default front and back line styles
fmt={'color':'k','linestyle':'solid',}
#fmtb = {'color':'k','linestyle':'dotted',}
fmtb={'color':'k','linestyle':'solid',} #same linestyle for backside



#-----------------------------------------------------
#borders
fmtO={'color':'red','linestyle':'solid',}
#ax=nspy.draw_longitude(ax, pi/2,  0.0, pi, fmt=fmt)
#ax=nspy.draw_longitude(ax, -pi/2, 0.0, pi, fmt=fmt)

Nlines = 20
for sweeps in np.linspace(0.0, 2.0*pi, Nlines):
    ax=nspy.draw_longitude(ax, sweeps,  0.0, pi/2.0, fmt=fmt)
    ax=nspy.draw_longitude(ax, sweeps,  pi/2.0, pi, fmt=fmt) #bottom

    ax=nspy.draw_longitude(ax, sweeps,  0.0, pi/2.0, fmt=fmtb, backside=True)
    ax=nspy.draw_longitude(ax, sweeps,  pi/2.0, pi,  fmt=fmtb, backside=True) #bottom


for sweeps in np.linspace(0.0, pi, Nlines):
    ax=nspy.draw_latitude(ax, sweeps, fmt=fmt)
    ax=nspy.draw_latitude(ax, sweeps, backside=True, fmt=fmtb)



#-------------------------------------------------- 
# Second star
#xyshift = (0.9, 1.2)
#nspy.req = 0.5
#nspy.o21 = 0.3
#
#for sweeps in np.linspace(0.0, 2.0*pi, 20):
#    ax=nspy.draw_longitude(ax, sweeps,  0.0, pi/2.0, fmt=fmt,xyshift=xyshift)
#    ax=nspy.draw_longitude(ax, sweeps,  pi/2.0, pi, fmt=fmt, xyshift=xyshift) #bottom
#
#    ax=nspy.draw_longitude(ax, sweeps,  0.0, pi/2.0, fmt=fmtb, backside=True, xyshift=xyshift)
#    ax=nspy.draw_longitude(ax, sweeps,  pi/2.0, pi,  fmt=fmtb, backside=True, xyshift=xyshift) #bottom
#
#for sweeps in np.linspace(0.0, pi, 20):
#    ax=nspy.draw_latitude(ax, sweeps, fmt=fmt, xyshift=xyshift)
#    ax=nspy.draw_latitude(ax, sweeps, backside=True, fmt=fmtb, xyshift=xyshift)
#
#
#
#nspy.req = 1.0
#nspy.o21 = 0.7


#-----------------------------------------------------
#circular spot
#rho = 0.2
#ax = nspy.draw_spot(ax, spot_phi, spot_theta, rho)

cmap = cm.get_cmap('inferno')
#cmap = cm.get_cmap('plasma')

#spot_theta=pi/3.0
spot_theta=0.0
spot_phi=0.4
fmtS={'facecolor':'red', 'alpha':0.4,}
Ns = 10

for i, rho in enumerate( np.linspace(1.8, 0.3, Ns) ):
    fmtS['facecolor'] = cmap(i / float(Ns-1) )

    ax = nspy.draw_spot(ax, spot_phi, spot_theta, rho, fmt=fmtS)


#spot_theta=pi
#spot_phi=0.4
#fmtS={'facecolor':'blue', 'alpha':0.4,}
#for rho in np.linspace(0.05, 1.2, 5):
#    ax = nspy.draw_spot(ax, spot_phi, spot_theta, rho, fmt=fmtS)

#
#spot_theta=pi/3.1
#spot_phi=-0.4
#fmtS={'facecolor':'blue', 'alpha':0.4,}
#for rho in np.linspace(0.05, 0.8, 5):
#    ax = nspy.draw_spot(ax, spot_phi, spot_theta, rho, fmt=fmtS)
#
#spot_theta=pi/3.2
#spot_phi=0.4 - 0.8*2
#fmtS={'facecolor':'darkorange', 'alpha':0.4,}
#for rho in np.linspace(0.05, 0.8, 5):
#    ax = nspy.draw_spot(ax, spot_phi, spot_theta, rho, fmt=fmtS)
#
#spot_theta=pi/3.3
#spot_phi=0.4 - 0.8*3
#fmtS={'facecolor':'gray', 'alpha':0.4,}
#for rho in np.linspace(0.05, 0.8, 5):
#    ax = nspy.draw_spot(ax, spot_phi, spot_theta, rho, fmt=fmtS)
#
#spot_theta=pi/3.4
#spot_phi=0.4 - 0.8*4
#fmtS={'facecolor':'green', 'alpha':0.4,}
#for rho in np.linspace(0.05, 0.8, 5):
#    ax = nspy.draw_spot(ax, spot_phi, spot_theta, rho, fmt=fmtS)
#
#spot_theta=pi/3.5
#spot_phi=0.4 - 0.8*6
#fmtS={'facecolor':'cyan', 'alpha':0.4,}
#for rho in np.linspace(0.05, 0.8, 5):
#    ax = nspy.draw_spot(ax, spot_phi, spot_theta, rho, fmt=fmtS)

#-----------------------------------------------------
# equatorial plane wave
#wave_start= -pi/2 -0.6 + 0.4,
#wave_stop = pi/2 -0.6 - 0.4,
wave_start= 0.0
wave_stop = 2.0*pi

#wave_start=  pi/3
#wave_stop = -pi/3


#ax = nspy.draw_wave(ax,
#                    phi_start= wave_start,
#                    phi_stop = wave_stop,
#                    Nphi = 50,
#                    r_start = 1.0,
#                    r_stop = 3.5,
#                    theta = pi/2
#                    )

#ax = nspy.draw_wave(ax,
#                    phi_start= -pi/2 -0.6 + 0.4 -pi,
#                    phi_stop =  pi/2 -0.6 - 0.4 -pi,
#                    Nphi = 25,
#                    r_start = 1.2,
#                    r_stop = 8.0,
#                    theta = pi/2
#                    )


#csfont = {'fontname':'Helvetica'}
#ax.text(-1.5, 0.5, "G", size=40, **csfont)



path_T=[[(0.0, 1.0), (1.0, 1.0),
      (1.0, 0.8), (0.7, 0.8),
      (0.7, 0.0), (0.3, 0.0),
      (0.3, 0.8), (0.0, 0.8),
      (0.0, 1.0) ]]

path_E=[[(0.0, 1.0), (1.1, 1.0),
        (1.1, 0.8), (0.4, 0.8),
        (0.4, 0.6), (1.1, 0.6),
        (1.1, 0.4), (0.4, 0.4),
        (0.4, 0.2), (1.1, 0.2),
        (1.1, 0.0), (0.0, 0.0),
        (0.0, 1.0)]]

path_A=[[#(0.0, 1.0), (1.0, 1.0),
        (0.0, 0.8),
        (0.2, 1.0), (0.8, 1.0),
        (1.0, 0.8),
        (1.0, 0.0), (0.8, 0.0),
        (0.8, 0.4), (0.2, 0.4),
        (0.2, 0.0), (0.0, 0.0),
        (0.0, 0.8)],
#path_A2=
        [(0.3, 0.7),(0.7, 0.7),
         (0.7, 0.6), (0.3,0.6),
         (0.3, 0.7)]
        ]

path_Ab=[[#(0.0, 1.0), (1.0, 1.0),
        (0.0, 0.8),
        (0.2, 1.0), (0.8, 1.0),
        (1.0, 0.8),
        (0.7,0.7),(0.3,0.7),(0.3,0.6),(0.7,0.6),(0.7,0.7),(1.0,0.8),
        (1.0, 0.0), (0.8, 0.0),
        (0.8, 0.4), (0.2, 0.4),
        (0.2, 0.0), (0.0, 0.0),
        (0.0, 0.8)]]


path_R=[[(0.0, 1.0), (0.6, 1.0),
        (1.0, 0.8),
        (1.0, 0.5), (0.9, 0.4),
        (1.0, 0.3), (1.0, 0.0),
        (0.7, 0.0), (0.7, 0.3),
        (0.4, 0.3), (0.4, 0.0),
        (0.0, 0.0), (0.0, 1.0)],
#path_R2=
        [(0.4, 0.8),(0.7, 0.8),
         (0.7, 0.5), (0.4,0.5),
         (0.4, 0.8)]
        ]

path_Rb=[[(0.0, 1.0), (0.6, 1.0),
        (0.7,0.8),(0.4,0.8),(0.4,0.5),(0.7,0.5),(0.7,0.8),(0.6,1.0),
        (1.0, 0.8),
        (1.0, 0.5), (0.9, 0.4),
        (1.0, 0.3), (1.0, 0.0),
        (0.7, 0.0), (0.7, 0.3),
        (0.4, 0.3), (0.4, 0.0),
        (0.0, 0.0), (0.0, 1.0)]]


path_G=[[#(0.4, 1.0), (0.6, 1.0),
        (0.0, 1.0), (0.6, 1.0),
        (1.0, 0.8), (1.0, 0.6),
        (0.6, 0.6), (0.6, 0.8),
        (0.4, 0.8), (0.4, 0.2),
        (0.7, 0.2), (0.7, 0.3),
        (0.5, 0.3), (0.5, 0.4),
        (1.0, 0.4), (1.0, 0.2),
        (0.6, 0.0), (0.4, 0.0),
        (0.0, 0.0), (0.0, 1.0)]]


fmt ={'color':'black','linestyle':'solid', 'alpha':1.0, 'linewidth':2.0}
fmtb={'color':'black','linestyle':'solid', 'alpha':0.9, 'linewidth':2.0}
fmtf={'facecolor':'blue', 'alpha':0.10, 'edgecolor':None}


#logo_start = -2.6
#rdist = 2.6

#playign around
logo_start = -2.4
#logo_start = -2.2
rdist = 2.6


#ax = nspy.draw_letters(ax, phi=logo_start -0.30*0, theta=1.6, rs=rdist, fmt=fmt)
#cmap = cm.get_cmap('viridis')
cmap = cm.get_cmap('plasma')
#cmap = cm.get_cmap('Reds')
#cmap = pal.wesanderson.Zissou_5.mpl_colormap
#cmap = pal.wesanderson.Aquatic1_5.mpl_colormap





letters = [path_G, path_Rb, path_E, path_Ab, path_T]
for i, letter in enumerate(letters):
    
    if i == 1:
        skips = [1,6]
    elif i == 3:
        skips = [3,8]

    else:
        skips = []

    #fmtf['facecolor'] = cmap(1.0/5.0 +  i/5.0)

    #ax = nspy.draw_letters(ax, phi=logo_start -0.30*i, theta=1.6, rs=rdist, paths=letter, fmt=fmt,  fmtb=fmtb, fmtf=fmtf, cmap=cmap)
    ax = nspy.draw_letters(ax, 
                           phi=logo_start -0.40*i, 
                           #phi=logo_start -0.52*i, 
                           theta=1.65, 
                           rs=rdist, 
                           paths=letter, 
                           fmt=fmt,  
                           fmtb=fmtb, 
                           fmtf=fmtf, 
                           cmap=cmap, skips=skips
                           )
    

#ax = nspy.draw_letters(ax, phi=logo_start -0.30*0, theta=1.6, rs=rdist, paths=path_G, fmt=fmt,  fmtb=fmtb, fmtf=fmtf)
#ax = nspy.draw_letters(ax, phi=logo_start -0.30*1, theta=1.6, rs=rdist, paths=path_Rb, fmt=fmt, fmtb=fmtb, fmtf=fmtf)
#ax = nspy.draw_letters(ax, phi=logo_start -0.30*2, theta=1.6, rs=rdist, paths=path_E, fmt=fmt,  fmtb=fmtb, fmtf=fmtf)
#ax = nspy.draw_letters(ax, phi=logo_start -0.30*3, theta=1.6, rs=rdist, paths=path_Ab, fmt=fmt, fmtb=fmtb, fmtf=fmtf)
#ax = nspy.draw_letters(ax, phi=logo_start -0.30*4, theta=1.6, rs=rdist, paths=path_T, fmt=fmt,  fmtb=fmtb, fmtf=fmtf)






savefig('fig2.png', bbox_inches='tight')
