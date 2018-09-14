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
#ax.axis('off')
ax.set_aspect('equal')

#ax.set_title('Pulsar', size=18)

#image dimensions
xmin = -9.6
xmax =  9.6
ymin = -5.8
ymax =  5.2

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


#coordinates of observer
#obs_theta=incl
#obs_theta=pi/2.7
#obs_phi=-pi/1.5


#default front and back line styles
fmt={'color':'k','linestyle':'solid',}
fmtb = {'color':'k','linestyle':'dotted',}
#fmtb={'color':'k','linestyle':'solid',} #same linestyle for backside

#borders
fmtO={'color':'red','linestyle':'solid',}



#-----------------------------------------------------
# STAR

if True:
    Nlines = 5
    for sweeps in np.linspace(0.0, 2.0*pi, Nlines):
        ax=nspy.draw_longitude(ax, sweeps,  0.0, pi/2.0, fmt=fmt)
        ax=nspy.draw_longitude(ax, sweeps,  pi/2.0, pi, fmt=fmt) #bottom
    
        ax=nspy.draw_longitude(ax, sweeps,  0.0, pi/2.0, fmt=fmtb, backside=True)
        ax=nspy.draw_longitude(ax, sweeps,  pi/2.0, pi,  fmt=fmtb, backside=True) #bottom
    
    
    for sweeps in np.linspace(0.0, pi, Nlines):
        ax=nspy.draw_latitude(ax, sweeps, fmt=fmt)
        ax=nspy.draw_latitude(ax, sweeps, backside=True, fmt=fmtb)
    


#-----------------------------------------------------
# DISK

if False:
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
    
    
#-----------------------------------------------------
# Light cylinder

if True:
    fmtLC={'color':'k','linestyle':'solid',}

    H   = 4.0 #cylinder height
    RLC = 3.5 #location of light cylinder   


    #bottom
    ax = nspy.draw_latitude(ax,
                            pi/2.,
                            rfac=RLC,
                            backside=True,
                            fmt = fmtLC,
                            xyshift=(0.0, -H)
                            )
    ax = nspy.draw_latitude(ax,
                            pi/2.,
                            rfac=RLC,
                            backside=False,
                            fmt = fmtLC,
                            xyshift=(0.0, -H)
                            )
                             
    #top
    ax = nspy.draw_latitude(ax,
                            pi/2.,
                            rfac=RLC,
                            backside=True,
                            fmt = fmtLC,
                            xyshift=(0.0, +H)
                            )
    ax = nspy.draw_latitude(ax,
                            pi/2.,
                            rfac=RLC,
                            backside=False,
                            fmt = fmtLC,
                            xyshift=(0.0, +H)
                            )



#-----------------------------------------------------
# B field

if True:
    fmtB={'color':'r','linestyle':'solid',}

    phiB = pi/2.
    theinc = pi/8.0 #pulsar inclination

    #ax = nspy.draw_flow_line(ax,
    #                         r_start=1.0,
    #                         r_stop=RLC + 0.5,
    #                         phi = phiB,
    #                         theta= 0.01,
    #                         fmt=fmtB
    #                         )

    req1 = 1.5 #first field line (where it crosses equator)
    req2 = RLC #last field line 
    for r in np.linspace(req1, req2, 5):
        ax = nspy.draw_dipole_field_line(ax,r,phi=phiB, tinc=theinc, fmt=fmtB)


    #reconnecting field lines
    req1 = RLC + 0.1 #first field line (where it crosses equator)
    req2 = 8.0       #last field line 
    for r in np.linspace(req1, req2, 1):
        ax = nspy.draw_open_field_line(ax,r,RLC, phi=phiB, tinc=theinc, fmt=fmtB)


savefig('fig_pulsar.png', bbox_inches='tight')
