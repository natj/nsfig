import numpy as np
import math
from pylab import *
import matplotlib.patches as patches
from scipy.optimize import minimize_scalar

import nspy as nspy

from matplotlib import cm

import palettable as pal



## Plot
#fig = plt.figure(figsize=(3.54, 2.3)) #single column fig
fig = plt.figure(figsize=(7.48, 3.5))  #two column figure

plt.rc('font', family='serif')
plt.rc('xtick', labelsize=7)
plt.rc('ytick', labelsize=7)
plt.rc('axes', labelsize=7)

gs = GridSpec(1, 1)
#gs.update(wspace = 0.34)
#gs.update(hspace = 0.4)

ax = subplot(gs[0,0])
ax.axis('off')
ax.set_aspect('equal')

#ax.set_title('Pulsar', size=18)

#image dimensions
xmin = -15.0
xmax =  15.0
ymin = -7.2
ymax =  7.4

xlen = xmax-xmin
ylen = ymax-ymin
aspr=ylen/xlen

ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)


#variables
#nspy.incl = pi/2.1 
nspy.req = 1
nspy.u = 0.0
nspy.rg = nspy.u * nspy.req
nspy.muc = -nspy.rg / (nspy.req - nspy.rg)
nspy.o21 = 0.20


#coordinates of observer
#obs_theta=incl
#obs_theta=pi/2.7
#obs_phi=-pi/1.5


#default front and back line styles
fmt={'color':'k','linestyle':'solid', 'lw': 0.8}
fmtb = {'color':'k','linestyle':'dotted', 'lw': 0.8}
#fmtb={'color':'k','linestyle':'solid',} #same linestyle for backside

#borders
fmtO={'color':'red','linestyle':'solid',}



#-----------------------------------------------------
# STAR

if True:

    #coloring
    if True:
        cmap = cm.get_cmap('inferno')
        #cmap = cm.get_cmap('plasma')
        
        #spot_theta=pi/3.0
        spot_theta=0.0
        spot_phi=0.4
        fmtS={'facecolor':'red', 'alpha':0.4,}
        Ns = 30
        
        for i, rho in enumerate( np.linspace(2.8, 0.3, Ns) ):
            fmtS['facecolor'] = cmap(i / float(Ns-1) )
            ax = nspy.draw_spot(ax, spot_phi, spot_theta, rho, fmt=fmtS)

    #draw outline
    if False:
        Nlines = 10
        for sweeps in np.linspace(0.0, 2.0*pi, Nlines):
            ax=nspy.draw_longitude(ax, sweeps,  0.0, pi/2.0, fmt=fmt)
            ax=nspy.draw_longitude(ax, sweeps,  pi/2.0, pi, fmt=fmt) #bottom
        
            ax=nspy.draw_longitude(ax, sweeps,  0.0, pi/2.0, fmt=fmtb, backside=True)
            ax=nspy.draw_longitude(ax, sweeps,  pi/2.0, pi,  fmt=fmtb, backside=True) #bottom
        
        
        Nlines = 10
        for sweeps in np.linspace(0.0, pi, Nlines):
            ax=nspy.draw_latitude(ax, sweeps, fmt=fmt)
            ax=nspy.draw_latitude(ax, sweeps, backside=True, fmt=fmtb)

#-----------------------------------------------------
# DISK

if True:
    fmtD={'color':'k','linestyle':'solid', 'lw': 0.7, 'alpha':0.6}


    # equatorial disk
    #slice_thick = 0.8
    slice_thick = pi/2. - 0.2
    disk_start= 0.0    + slice_thick
    disk_stop = 2.0*pi - slice_thick
    
    diskH=0.300
    inner_disk = 10.0
    outer_disk = 16.0
    
    
    ax = nspy.draw_disk(ax,
                        phi_start= disk_start,
                        phi_stop = disk_stop,
                        Nphi = 200,
                        r_start = inner_disk,
                        r_stop = outer_disk,
                        theta = pi/2-diskH,
                        fmt=fmtD
                        )
    
    #disk rim
    ax = nspy.draw_disk(ax,
                        phi_start= disk_start,
                        phi_stop = disk_stop,
                        Nphi = 80,
                        r_start = inner_disk,
                        r_stop = inner_disk+0.1,
                        theta = pi/2.0+diskH,
                        fmt=fmtD
                        )

    #interior pattern
    ax = nspy.draw_disk_interior(
                        ax,
                        phi_start= disk_start,
                        phi_stop = disk_stop,
                        Nphi = 150,
                        Nthe = 5,
                        r_start = inner_disk,
                        theta_up = pi/2.0-diskH,
                        theta_dw = pi/2.0+diskH,
                        fmt=fmtD
                        )
    
    #right bottom
    ax = nspy.draw_disk(ax,
                        phi_start= disk_start,
                        phi_stop = disk_start+0.02,
                        Nphi = 10,
                        r_start = inner_disk,
                        r_stop = outer_disk,
                        theta = pi/2.0+diskH,
                        fmt=fmtD
                        )
    
    #left bottom
    ax = nspy.draw_disk(ax,
                        phi_start= disk_stop,
                        phi_stop = disk_stop+0.02,
                        Nphi = 10,
                        r_start = inner_disk,
                        r_stop = outer_disk,
                        theta = pi/2.0+diskH,
                        fmt=fmtD
                        )

# DISK WEDGES
if False:

    # equatorial disk
    #RIGHT:

    for i in [1,2]:
        if i == 1:
            slice_thick = 0.6
            disk_start= pi/2.+0.2  + slice_thick
            disk_stop = pi/2.+0.1 
        else:
            slice_thick = 0.6
            disk_start= -pi/2.+0.2  + slice_thick
            disk_stop = -pi/2.+0.1 
    
        diskH=0.250
        inner_disk = 10.0
        outer_disk = 16.0
        
        ax = nspy.draw_disk(ax,
                            phi_start= disk_start,
                            phi_stop = disk_stop,
                            Nphi = 50,
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

        #interior pattern
        ax = nspy.draw_disk_interior(
                            ax,
                            phi_start= disk_start,
                            phi_stop = disk_stop,
                            Nphi = 40,
                            Nthe = 20,
                            r_start = inner_disk,
                            theta_up = pi/2.0-diskH,
                            theta_dw = pi/2.0+diskH
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
    fmtLC={'color':'b','linestyle':'dashed', 'lw': 1.5}

    H   = 5.5 #cylinder height
    RLC = 4.5 #location of light cylinder   


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

    #sides
    x1 = nspy.y(-pi/2., pi/2.)*RLC 
    y1 = nspy.z(-pi/2., pi/2.)*RLC
    ax.plot([x1, x1], [y1-H, y1+H], **fmtLC)

    x2 = nspy.y(+pi/2., pi/2.)*RLC 
    y2 = nspy.z(+pi/2., pi/2.)*RLC
    ax.plot([x2, x2], [y2-H, y2+H], **fmtLC)


#-----------------------------------------------------
# Rotation arrow
if True:
    obs_phi=-pi/1.5

    fmta= {'color':'b', 'linestyle':'solid', 'lw':1.5, 'head_width': 0.08, 'head_length': 0.16,}

    #draw xyz axis
    #ax = nspy.draw_axis(ax, obs_phi,      pi/2, rfac=1.4)
    #ax = nspy.draw_axis(ax, obs_phi+pi/2, pi/2, rfac=1.4)
    ax = nspy.draw_axis(ax, 0.0, 0.0 , startp=(0.0, -5.5), rfac=16.0, fmt=fmta)
    
    #rotation direction, i.e. curly arrow
    phistart=-pi/2-0.8
    phistop=pi/2+0.5
    wsize=0.15
    wheight=-8.0
    
    ax=nspy.draw_latitude(ax, wsize, fmt={'color':'b'}, start=phistart, stop=phistop, rfac=wheight)
    yy1=nspy.y(phistop,wsize)*wheight
    zz1=nspy.z(phistop,wsize)*wheight
    
    yy2 = -0.01
    zz2 = 0.002
    
    ax.add_patch(patches.FancyArrow(
                yy1, zz1,
                yy2, zz2,
                **fmta
                ))
    

#-----------------------------------------------------
# B field

if True:
    fmtB={'color':'k','linestyle':'solid', 'lw': 1.1}

    phiB = pi/2. + 0.0
    theinc = pi/15.0 #pulsar inclination

    #closed field lines

    req1 = 1.5 #first field line (where it crosses equator)
    req2 = 0.9*RLC #last field line 
    for r in np.linspace(req1, req2, 3):
        ax = nspy.draw_dipole_field_line(ax,r,phi=phiB, tinc=theinc, fmt=fmtB)


    # open field lines
    req1 = 1.2*RLC
    req2 = 5.0*RLC
    #for r in np.linspace(req1, req2, 3):
    for r in [1.8*RLC, 3.0*RLC, 30.0*RLC]:
        ax = nspy.draw_open_field_line(
                ax,
                r,
                RLC, 
                phi=phiB, 
                tinc=theinc, 
                fmt=fmtB)


    #separatrix field line
    fmtRC={'color':'r','linestyle':'dashed', 'lw': 1.1}
    ax = nspy.draw_dipole_field_line(ax, RLC,phi=phiB, tinc=theinc, fmt=fmtRC)

    #reconnecting field lines
    req1 = RLC + 0.1 #first field line (where it crosses equator)
    req2 = 8.0       #last field line 
    for r in np.linspace(req1, req2, 1):
        ax = nspy.draw_reconnecting_field_line(
                ax,
                r,
                RLC, 
                phi=phiB, 
                tinc=theinc, 
                sheet_len = 2.2,
                fmt=fmtRC, )



plt.subplots_adjust(left=0.0, bottom=0.0, right=0.99, top=0.99, wspace=0.0, hspace=0.0)
savefig('fig_pulsar.png')
savefig('fig_pulsar.pdf')
