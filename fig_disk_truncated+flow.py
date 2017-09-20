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

ax.set_title('Hard/low state', size=18)

#image dimensions
#mlim= 1.3
xmin = -5.6
xmax =  5.6
ymin = -2.8
ymax =  3.2

xlen = xmax-xmin
ylen = ymax-ymin
aspr=ylen/xlen

ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)


#variables
nspy.incl = pi/2.5 
nspy.req = 1
nspy.u = 0.00
nspy.rg = nspy.u * nspy.req
nspy.muc = -nspy.rg / (nspy.req - nspy.rg)
nspy.o21 = 0.10


#coordinates of spot and observer
spot_theta=pi/5
spot_phi=0.67

#obs_theta=incl
obs_theta=pi/2.7
obs_phi=-pi/1.5

#default front and back line styles
fmt={'color':'k','linestyle':'solid',}
fmtb = {'color':'k','linestyle':'dotted',}
#fmtb={'color':'k','linestyle':'solid',} #same linestyle for backside



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



#-----------------------------------------------------
# equatorial disk
slice_thick = 0.8
#slice_thick = 0.0
disk_start= 0.0    + slice_thick
disk_stop = 2.0*pi - slice_thick


diskH=0.010

inner_disk = 4.0
outer_disk = 8.0


ax = nspy.draw_disk(ax,
                    phi_start= disk_start,
                    phi_stop = disk_stop,
                    Nphi = 400,
                    r_start = inner_disk,
                    r_stop = outer_disk,
                    theta = pi/2-diskH,
                    )

#disk rim
ax = nspy.draw_disk(ax,
                    phi_start= disk_start,
                    phi_stop = disk_stop,
                    Nphi = 80,
                    r_start = inner_disk,
                    r_stop = inner_disk+0.1,
                    theta = pi/2.0+diskH
                    )

#right bottom
ax = nspy.draw_disk(ax,
                    phi_start= disk_start,
                    phi_stop = disk_start+0.02,
                    Nphi = 10,
                    r_start = inner_disk,
                    r_stop = outer_disk,
                    theta = pi/2.0+diskH
                    )

#left bottom
ax = nspy.draw_disk(ax,
                    phi_start= disk_stop,
                    phi_stop = disk_stop+0.02,
                    Nphi = 10,
                    r_start = inner_disk,
                    r_stop = outer_disk,
                    theta = pi/2.0+diskH
                    )


#ax = nspy.draw_flow_line(ax,
#                r_start = 1.0,
#                r_stop  = inner_disk,
#                phi = 0.8,
#                )


ax = nspy.draw_turbulent_flow(ax,
                r_start = 1.0,
                r_stop  = inner_disk,
                phi_start= disk_start,
                phi_stop = disk_stop,
                )


savefig('fig_disk_hard.png', bbox_inches='tight')
