import numpy as np
import math
from pylab import *
import matplotlib.patches as patches
from scipy.optimize import minimize_scalar

from matplotlib.path import Path
from matplotlib.patches import PathPatch

from numba import jit
from numba import float64

#variables
#incl = pi/8
#incl = pi/2.25
incl = pi/2.6
req = 1
u = 0.0

rg = u*req
muc = -rg/(req-rg)
o21 = 0.3

#@jit
#def R(theta):
#    return 1.0

@jit(float64(float64))
def R(theta):
    return 1.0*(1-o21*cos(theta)**2)

@jit(float64(float64))
def dR(theta):
    return 2.0*o21*sin(theta)*cos(theta)

@jit(float64(float64))
def fa(theta):
    return dR(theta)/R(theta)

@jit(float64(float64))
def cosg(theta):
    return (1/sqrt(1-u))/sqrt(1+fa(theta)**2)

@jit(float64(float64))
def sing(theta):
    return fa(theta)*cosg(theta)

@jit(float64(float64, float64))
def mu(phi, theta):
    return cos(incl)*cos(theta)+cos(phi)*sin(incl)*sin(theta)

@jit(float64(float64, float64))
def bimp(phi, theta):
    return sqrt(req*(req+rg+req*mu(phi,theta)-rg*mu(phi,theta))/(1+mu(phi,theta)))

@jit(float64(float64, float64))
def x(phi, theta):
    return bimp(phi,theta)*R(theta)*(cos(incl)*cos(theta)+cos(phi)*sin(incl)*sin(theta))

@jit(float64(float64, float64))
def y(phi, theta):
    return bimp(phi,theta)*R(theta)*(sin(theta)*sin(phi))

@jit(float64(float64, float64))
def z(phi, theta):
    return  bimp(phi,theta)*R(theta)*(cos(theta)*sin(incl)-cos(incl)*cos(phi)*sin(theta))

@jit(float64(float64, float64, float64, float64, float64))
def xfull(phi,theta, a,b,c):
    return bimp(phi, theta)*R(theta)* (cos(incl)*(a*cos(theta)-c*sin(theta)) + sin(incl)*(c*cos(theta)*cos(phi) + a*cos(phi)*sin(theta) - b*sin(phi)))

@jit(float64(float64, float64, float64, float64, float64))
def yfull(phi, theta, a,b,c):
    return bimp(phi, theta)*R(theta)* (b*cos(phi) + (c*cos(theta) + a*sin(theta))*sin(phi))

@jit(float64(float64, float64, float64, float64, float64))
def zfull(phi, theta, a,b,c):
    return bimp(phi, theta)*R(theta)* (sin(incl)*(a*cos(theta)-c*sin(theta)) + cos(incl)*(-cos(phi)*(c*cos(theta) + a*sin(theta)) + b*sin(phi)))



def draw_longitude(ax,
                   phi,
                   start=0,
                   stop=pi/2,
                   rfac=1.0,
                   backside=False,
                   fmt={'color':'k','linestyle':'solid',},
                   xyshift=(0.0,0.0)
                   ):

    xx = []
    yy = []

    for theta in np.linspace(start, stop, 200):
        if not(backside):
            if mu(phi,theta) >= muc:
                xx.append(y(phi,theta)*rfac + xyshift[0])
                yy.append(z(phi,theta)*rfac + xyshift[1])
            else:
                lines, = ax.plot(xx, yy)
                lines.set(**fmt)
                xx = []
                yy = []
        else:
            if mu(phi,theta) <= muc:
                xx.append(y(phi,theta)*rfac + xyshift[0])
                yy.append(z(phi,theta)*rfac + xyshift[1])
            else:
                lines, = ax.plot(xx, yy)
                lines.set(**fmt)
                xx = []
                yy = []

    lines, = ax.plot(xx,yy)
    lines.set(**fmt)

    return ax


def draw_latitude(ax,
                   theta,
                   start=-pi,
                   stop=pi,
                   rfac=1.0,
                   backside=False,
                   fmt={'color':'k','linestyle':'solid',},
                   xyshift=(0.0,0.0)
                   ):

    xx = []
    yy = []

    for phi in np.linspace(start, stop, 200):
        if not(backside):
            if mu(phi,theta) >= muc:
                xx.append(y(phi,theta)*rfac + xyshift[0])
                yy.append(z(phi,theta)*rfac + xyshift[1])
            else:
                lines, = ax.plot(xx, yy)
                lines.set(**fmt)
                xx = []
                yy = []
        else:
            if mu(phi,theta) <= muc:
                xx.append(y(phi,theta)*rfac + xyshift[0])
                yy.append(z(phi,theta)*rfac + xyshift[1])
            else:
                lines, = ax.plot(xx, yy)
                lines.set(**fmt)
                xx = []
                yy = []

    lines, = ax.plot(xx, yy)
    lines.set(**fmt)

    return ax

def draw_outline(ax,
                 fmt={'color':'k','linestyle':'solid',},
                 ):

    def rootf(theta1, phi1):
        return abs(mu(phi1, theta1))-muc

    xx = []
    yy = []

    #print "outline:"
    for phi in np.linspace(pi/2, 3*pi/2, 100):
        #print "phi:", phi
        thetamin = pi/2
        thetamax = 0.0

        res = minimize_scalar(rootf, bounds=(thetamin, thetamax), args=(phi,))
        #print "minim:", res.x, " ", res.fun

        theta2 = res.x

        #ax.plot([y(phi, theta2)], [z(phi, theta2)], "b.")
        xx.append( y(phi, theta2))
        yy.append( z(phi, theta2))

    ax.plot(xx, yy, "r-")

    return ax


def draw_radial(ax, phi, theta,
                fmt={'color':'k','linestyle':'solid',},
                rfac=1.0):
    xx = [0.0]
    yy = [0.0]

    xx.append(y(phi, theta)*rfac)
    yy.append(z(phi, theta)*rfac)

    lines, = ax.plot(xx, yy)
    lines.set(**fmt)

    return ax


def draw_axis(ax, phi, theta,
              startp=(0.0, 0.0),
              fmt= {'color':'k', 'linestyle':'solid', 'lw':1.5, 'head_width': 0.04, 'head_length': 0.08,},
              rfac=1.0):

    ax.add_patch(patches.FancyArrow(
              startp[0],startp[1],
              y(phi, theta)*rfac, z(phi, theta)*rfac,
              zorder=1,
              **fmt
              ))

    return ax


#normal vector
def draw_normal(ax, phi, theta,
                fmt={'color':'black', 'linestyle':'solid', 'lw':1.8, 'head_width': 0.02, 'head_length': 0.05,},
                alen=0.3
                ):

    spot_cosg = cosg(theta)
    spot_sing = sing(theta)

    normal_y1 = y(phi, theta)*1.0 
    normal_z1 = z(phi, theta)*1.0,
    normal_y2 = yfull(phi, theta, spot_cosg, 0.0, -spot_sing)*alen
    normal_z2 = zfull(phi, theta, spot_cosg, 0.0, -spot_sing)*alen

    ax.add_patch(patches.FancyArrow(
                normal_y1,      normal_z1,
                normal_y2,      normal_z2,
                **fmt
                ))
    
    return ax


def draw_observer(ax, xx, yy, angle, sangle=0.3, alen=0.2, 
                  fmt={'color':'k','linestyle':'solid','lw':1.0}):

    x1 = xx + alen*cos(angle-sangle/2)
    y1 = yy + alen*sin(angle-sangle/2)

    x2 = xx + alen*cos(angle+sangle/2)
    y2 = yy + alen*sin(angle+sangle/2)

    ax.plot([xx,x1],[yy,y1],**fmt)
    ax.plot([xx,x2],[yy,y2],**fmt)

    xsegm = []
    ysegm = []
    for theta in np.linspace(angle-sangle/2-0.1, angle+sangle/2+0.1, 10):
        xsegm.append(xx + alen*0.7*cos(theta))
        ysegm.append(yy + alen*0.7*sin(theta))
    ax.plot(xsegm, ysegm, **fmt)

    return ax

def draw_spot(ax, sphi, stheta, rho,
              fmt={'facecolor':'k', 'alpha':0.3,}):

    #TODO check visibility using mu

    xx = []
    yy = []
    for chi in np.linspace(0, 2*pi, 30):

        #transform from spot centric-coordinate system into the standard spherical coord
        xs = cos(rho)*cos(sphi)*sin(stheta)+sin(rho)*(cos(stheta)*cos(sphi)*cos(chi)-sin(sphi)*sin(chi))
        ys = cos(rho)*sin(stheta)*sin(sphi)+sin(rho)*(cos(stheta)*cos(chi)*sin(sphi)+cos(sphi)*sin(chi))
        zs = cos(stheta)*cos(rho)-cos(chi)*sin(stheta)*sin(rho)

        #transform to phi, theta
        phi = arctan2(ys,xs)
        theta = arccos(zs/sqrt(xs**2 + ys**2 + zs**2))

        #ax.plot([y(phi,theta)], [z(phi,theta)], "b.")

        xx.append(y(phi, theta))
        yy.append(z(phi, theta))

    #ax.plot(xx, yy, "b-")

    #hack to make zipped iterator to work with python3
    xxyy = list(zip(xx,yy))
    zipped_list = xxyy[:]

    ax.add_patch(Polygon(zipped_list,
                      closed=True, 
                      zorder=10,
                      **fmt))

    return ax


def draw_wave(ax,
              phi_start = 0.0,
              phi_stop = 2*pi,
              Nphi = 50,
              r_start = 1.2,
              r_stop = 8.0,
              theta = pi/2,
              fmt={'color':'k','linestyle':'solid', 'alpha':0.1}
              ):


    for phi in np.linspace(phi_start, phi_stop, Nphi):
        phi += 0.02
        xx = []
        yy = []
        for rfac in np.linspace(r_start, r_stop, 200):
            thetaI = theta + 0.03*sin(8.0*rfac)
            phiI = phi #+ 0.05*sin(8.0*rfac)
            xx.append(y(phiI,thetaI)*rfac)
            yy.append(z(phiI,thetaI)*rfac)
        ax.plot(xx,yy,**fmt)


    for rfac in np.linspace(r_start, r_stop, 400):
        thetaI = theta + 0.03*sin(8.0*rfac)
        xx = []
        yy = []
        for phi in np.linspace(phi_start, phi_stop, 200):
            phi += 0.02
            phiI = phi #+ 0.05*sin(8.0*rfac)
            xx.append(y(phiI,thetaI)*rfac)
            yy.append(z(phiI,thetaI)*rfac)
        ax.plot(xx,yy,**fmt)


    return ax


# draw normal flat accretion disk
def draw_disk(ax,
              phi_start = 0.0,
              phi_stop = 2*pi,
              Nphi = 50,
              r_start = 1.2,
              r_stop = 8.0,
              theta = pi/2,
              fmt={'color':'k','linestyle':'solid', 'alpha':0.5}
              ):


    for phi in np.linspace(phi_start, phi_stop, Nphi):
        phi += 0.02
        xx = []
        yy = []
        for rfac in np.linspace(r_start, r_stop, 200):
            thetaI = theta #+ 0.03*sin(8.0*rfac)
            phiI = phi #+ 0.05*sin(8.0*rfac)

            xx.append(y(phiI,thetaI)*rfac)
            yy.append(z(phiI,thetaI)*rfac)
        ax.plot(xx,yy,**fmt)


    for rfac in np.linspace(r_start, r_stop, 5):
        thetaI = theta #+ 0.03*sin(8.0*rfac)
        xx = []
        yy = []
        for phi in np.linspace(phi_start, phi_stop, 200):
            phi += 0.02
            phiI = phi #+ 0.05*sin(8.0*rfac)
            xx.append(y(phiI,thetaI)*rfac)
            yy.append(z(phiI,thetaI)*rfac)
        ax.plot(xx,yy,**fmt)

    return ax


# draw accreting pulsar disk with accretion fingers
def draw_acc_disk(ax,
              phi_start = 0.0,
              phi_stop = 2*pi,
              phi_finger1 = pi/2., #location of first acc finger
              Nphi = 50,
              r_start = 1.2,
              r_stop = 8.0,
              r_finger = 0.8, #radial extent of first acc finger
              theta = pi/2,
              fmt={'color':'k','linestyle':'solid', 'alpha':0.5}
              ):


    for phi in np.linspace(phi_start, phi_stop, Nphi):
        phi += 0.02
        xx = []
        yy = []

        modfac = finger_shape(phi, phi_finger1)

        r_startf = modfac*r_finger + (1.0 - modfac)*r_start
        for rfac in np.linspace(r_startf, r_stop, 200):
            thetaI = theta #+ 0.03*sin(8.0*rfac)
            phiI = phi #+ 0.05*sin(8.0*rfac)

            xx.append(y(phiI,thetaI)*rfac)
            yy.append(z(phiI,thetaI)*rfac)
        ax.plot(xx,yy,**fmt)


    for rfac in np.linspace(r_start, r_stop, 5):
        thetaI = theta #+ 0.03*sin(8.0*rfac)
        xx = []
        yy = []
        for phi in np.linspace(phi_start, phi_stop, 200):
            phi += 0.02
            phiI = phi #+ 0.05*sin(8.0*rfac)
            xx.append(y(phiI,thetaI)*rfac)
            yy.append(z(phiI,thetaI)*rfac)
        ax.plot(xx,yy,**fmt)

    return ax


def flare_shape_the(phi, phi_flare):
    modfac = 0.8*(np.abs(phi - phi_flare)/0.4)**2.0
    if modfac < 0.0:
        modfac = 0.0
    if modfac > 1.0:
        modfac = 1.0
    return modfac

def flare_shape_rad(rfac, r_start, r_stop):
    modfac_r = (rfac-r_start)/(r_stop-r_start)**3.0
    return modfac_r


# draw normal flat accretion disk
def draw_flaring_disk(ax,
              phi_start = 0.0,
              phi_stop = 2*pi,
              phi_flare = 0.0,
              Nphi = 50,
              r_start = 1.2,
              r_stop = 8.0,
              theta = pi/2,
              fmt={'color':'k','linestyle':'solid', 'alpha':0.5}
              ):

    theta_flare = theta - 0.4

    dphi = np.abs(phi_stop - phi_start)
    for phi in np.linspace(phi_start, phi_stop, Nphi):
        phi += 0.02
        xx = []
        yy = []

        modfac = flare_shape_the(phi, phi_flare)
        thetaI = theta_flare*(1.0-modfac) + modfac*theta

        for rfac in np.linspace(r_start, r_stop, 200):

            modfac_r = flare_shape_rad(rfac, r_start, r_stop)
            thetaI = thetaI*(1.0-modfac_r) + modfac_r*theta

            xx.append(y(phi,thetaI)*rfac)
            yy.append(z(phi,thetaI)*rfac)
        ax.plot(xx,yy,**fmt)


    for rfac in np.linspace(r_start, r_stop, 5):
        xx = []
        yy = []

        thetaI = theta #+ 0.03*sin(8.0*rfac)
        modfac_r = flare_shape_rad(rfac, r_start, r_stop)
        thetaI = thetaI*(1.0-modfac_r) + modfac_r*theta
        for phi in np.linspace(phi_start, phi_stop, 200):
            phi += 0.02
            phiI = phi #+ 0.05*sin(8.0*rfac)

            #modfac = flare_shape_the(phi, phi_flare)
            #thetaI = theta_flare*(1.0-modfac) + modfac*thetaI

            xx.append(y(phiI,thetaI)*rfac)
            yy.append(z(phiI,thetaI)*rfac)
        ax.plot(xx,yy,**fmt)

    return ax





def draw_disk_interior(
              ax,
              phi_start = 0.0,
              phi_stop = 2*pi,
              Nphi = 10,
              Nthe = 10,
              r_start = 1.2,
              theta_up = pi/2-0.1,
              theta_dw = pi/2+0.1,
              fmt={'color':'k','linestyle':'solid', 'alpha':0.5}
              ):

    rfac = r_start

    
    for phi in np.linspace(phi_start, phi_stop, Nphi):
        xx = []
        yy = []
        for the in np.linspace(theta_dw, theta_up, 20):
            xx.append(y(phi,the)*rfac)
            yy.append(z(phi,the)*rfac)
        ax.plot(xx,yy,**fmt)

    for the in np.linspace(theta_dw, theta_up, Nthe):
        xx = []
        yy = []
        for phi in np.linspace(phi_start, phi_stop, 100):
            xx.append(y(phi,the)*rfac)
            yy.append(z(phi,the)*rfac)
        ax.plot(xx,yy,**fmt)

    return ax


def finger_shape(phi, phi_finger):
    modfac = 1.0 - 0.2*((phi - phi_finger)/0.3)**2.0
    if modfac < 0.0:
        modfac = 0.0

    return modfac


# draw interior of accreting pulsar disk with fingers
def draw_acc_disk_interior(
              ax,
              phi_start = 0.0,
              phi_stop = 2*pi,
              phi_finger1 = pi/2.,     #location of first acc finger
              phi_finger2 = 3.0*pi/2., #location of second acc finger
              Nphi = 10,
              r_start = 1.2,
              r_finger = 0.8,
              theta_up = pi/2-0.1,
              theta_dw = pi/2+0.1,
              fmt={'color':'k','linestyle':'solid', 'alpha':0.5}
              ):

    dthe = np.abs(theta_up - theta_dw)
    for phi in np.linspace(phi_start, phi_stop, Nphi):
        xx = []
        yy = []

        if phi < pi:
            modfac = finger_shape(phi, phi_finger1)
        else:
            modfac = finger_shape(phi, phi_finger2)

        rfac = modfac*r_finger + (1.0 - modfac)*r_start

        for the in np.linspace(theta_dw, theta_up, 20):

            if phi < pi:
                modfac_t = ((the - theta_up)/dthe)**2.0
            else:
                modfac_t = ((the - theta_dw)/dthe)**2.0

            rfacm = rfac*modfac_t + (1.0-modfac_t)*r_start
            xx.append(y(phi,the)*rfacm)
            yy.append(z(phi,the)*rfacm)
        ax.plot(xx,yy,**fmt)

    return ax


def draw_flow_line(ax,
              r_start = 1.2,
              r_stop = 8.0,
              phi = pi/2,
              theta = pi/3.5,
              fmt={'color':'r','linestyle':'solid', 'alpha':0.2}
              ):

    #disk plane
    theta_ref = pi/2.0

    xx = []
    yy = []

    Nr = 20
    transf = np.linspace(0.0, 1.0, Nr)
    for i, rfac in enumerate(np.linspace(r_start, r_stop, Nr)):
        phiI = phi
        thetaI = (1.0 - transf[i])*theta + transf[i]*theta_ref

        xx.append(y(phiI,thetaI)*rfac)
        yy.append(z(phiI,thetaI)*rfac)

    ax.plot(xx,yy,**fmt)

    return ax


def randab(a,b):
    return a + (b-a)*np.random.rand()


def draw_turbulent_flow(ax,
              r_start = 1.2,
              r_stop = 8.0,
              phi_start = 0.0,
              phi_stop  = 2.0*pi,
              fmt={'color':'darkred','linestyle':'solid', 'alpha':0.15, 'linewidth':1.0}
              ):

    #Nlines = 500
    #for i in range(Nlines):
    #    theta = randab(0.0, pi)
    #    #theta = randab(0.0, pi/2.0)
    #    phi   = randab(phi_start, phi_stop)
    #    draw_flow_line(ax, r1, r2, phi, theta, fmt)

    r1 = r_start
    r2 = r_stop #+ randab(0.0, 0.5)

    phis = np.linspace(phi_start, phi_stop, 22)
    dphi = np.diff(phis)[0]/22.0
    dp = 0.0

    for i, theta in enumerate(np.linspace(0.0, pi, 22)):
        for phi in phis:
            if phi+dp < phi_stop:
                draw_flow_line(ax, r1, r2, phi+dp, theta, fmt)
        dp += dphi


    return ax 


yoff = 1.0

def draw_phi(ax, phi1, phi2, theta, r1, r2, fmt, fmtb):

    #front
    xx = []
    yy = []
    for phi in np.linspace(phi1, phi2, 100):
        xx.append(y(phi,theta)*r1)
        yy.append(z(phi,theta)*r1 + yoff )
    ax.plot(xx, yy, **fmt)

    #back
    xxB = []
    yyB = []
    for phi in np.linspace(phi1, phi2, 100):
        xxB.append(y(phi,theta)*r2)
        yyB.append(z(phi,theta)*r2 +yoff )
    ax.plot(xxB, yyB, **fmtb)

    return ax, xx, yy

def draw_theta(ax, the1, the2, phi, r1, r2, fmt, fmtb):

    #front
    xx = []
    yy = []
    for theta in np.linspace(the1, the2, 100):
        xx.append(y(phi,theta)*r1)
        yy.append(z(phi,theta)*r1 + yoff )
    ax.plot(xx, yy, **fmt)

    x1 = xx[1]
    y1 = yy[1]
    x2 = xx[-1]
    y2 = yy[-1]

    #back
    xxB = []
    yyB = []
    for theta in np.linspace(the1, the2, 100):
        xxB.append(y(phi,theta)*r2)
        yyB.append(z(phi,theta)*r2 +yoff )
    ax.plot(xxB, yyB, **fmtb)

    x3 = xxB[1]
    y3 = yyB[1]
    x4 = xxB[-1]
    y4 = yyB[-1]

    #connecting lines
    ax.plot([x1, x3], [y1, y3], **fmtb)
    ax.plot([x2, x4], [y2, y4], **fmtb)

    return ax, xx, yy

def draw_geo(ax, phi1, phi2, the1, the2, r1, r2, fmt, fmtb, plot=True):

    #front
    xx = []
    yy = []
    for (phi, theta) in zip(np.linspace(phi1, phi2, 100), np.linspace(the1, the2, 100)):
        xx.append(y(phi,theta)*r1)
        yy.append(z(phi,theta)*r1 + yoff)
        
    if plot:
        ax.plot(xx, yy, **fmt)

    #back
    xxB = []
    yyB = []
    for (phi, theta) in zip(np.linspace(phi1, phi2, 100), np.linspace(the1, the2, 100)):
        xxB.append(y(phi,theta)*r2)
        yyB.append(z(phi,theta)*r2 +yoff )

    if plot:
        ax.plot(xxB, yyB, **fmtb)

    return ax, xx, yy



def clip_xy(xx, yy, ylim):
    
    #yc = np.empty_like(yy)
    yc = np.array([])
    xc = np.array([])

    for j, y in enumerate(yy):
        #if y > ylim:
        #    yc[j] = ylim
        #else:
        #    yc[j] = yy[j]

        if y > ylim:
            yc = np.append(yc, yy[j])
            xc = np.append(xc, xx[j])

    return xc, yc


def draw_letters(ax,
              phi=0.0,
              theta=pi/2,
              rs=1.2,
              paths=[ [(0,0), (0,1),(1,1),(1,0),(0,0)] ],
              fmt ={'color':'k','linestyle':'solid', 'alpha':0.9},
              fmtb={'color':'k','linestyle':'solid', 'alpha':0.1},
              fmtf={'facecolor':'blue','alpha':0.9},
              cmap=cm.get_cmap('inferno'),
              skips=[]
              ):

    #dphi = 0.12
    #dthe = 0.30
    #dr   = 0.05

    #playing around
    dphi = 0.16
    #dthe = 0.29
    dthe = 0.29
    dr   = 0.05

    phi1 = phi - dphi
    phi2 = phi + dphi

    the1 = theta - dthe
    the2 = theta + dthe

    rs1 = rs - dr
    rs2 = rs + dr

    #xx = np.array([])
    #yy = np.array([])
    polys = []


    ##################################################
    #clipping edge for color shading
    thetaI = np.pi/2
    xx_edge = []
    yy_edge = []
    rfac = 2.0
    for phi in np.linspace(0.0, 2*np.pi, 200):
        phi += 0.02
        phiI = phi 
        xx_edge.append(y(phiI,thetaI)*rfac)
        yy_edge.append(z(phiI,thetaI)*rfac)
    ##################################################

   
    jj = 0
    for path in paths:
        xx = np.array([])
        yy = np.array([])

        for ii in range(len(path)-1):
            phia = phi1 + (1.0-path[ii][0]  )*2.0*dphi
            thea = the1 + (1.0-path[ii][1]  )*2.0*dthe
            phib = phi1 + (1.0-path[ii+1][0])*2.0*dphi
            theb = the1 + (1.0-path[ii+1][1])*2.0*dthe

            xxt = np.array([])
            yyt = np.array([])

            if phia == phib:
                ax, xxt, yyt = draw_theta(ax, thea, theb, phia, rs1, rs2, fmt, fmtb)
            elif thea == theb:
                ax, xxt, yyt = draw_phi(ax, phia, phib, thea, rs1, rs2, fmt, fmtb)
            else:
                if ii in skips:
                    ax, xxt, yyt = draw_geo(ax, phia, phib, thea, theb, rs1, rs2, fmt, fmtb, plot=False)
                else:
                    ax, xxt, yyt = draw_geo(ax, phia, phib, thea, theb, rs1, rs2, fmt, fmtb, plot=True)

            xx = np.append(xx, xxt)
            yy = np.append(yy, yyt)

        top = max(yy)
        bot = min(yy)

        for kk, ylim in enumerate(np.linspace(1.01*bot, 0.99*top, 20)): 
            #fmtf['facecolor'] = cmap( 1.0 - kk / 19.0)
            fmtf['facecolor'] = cmap( kk / 19.0)

            xxc,yyc = clip_xy(xx,yy, ylim)
            polys.append( Polygon(zip(xxc, yyc), closed=True, edgecolor='None', **fmtf) )



    #p = ax.add_patch(polys[0])
    #p = polys[0]
    for poly in polys:
        p = ax.add_patch(poly)
    

    #clipping edge for color shading
    thetaI = np.pi/2
    xx_edge = []
    yy_edge = []
    rfac = 2.0
    for phi in np.linspace(0.0, 2*np.pi, 200):
        phi += 0.02
        phiI = phi 
        xx_edge.append(y(phiI,thetaI)*rfac)
        yy_edge.append(z(phiI,thetaI)*rfac)
    #edge_path =  PathPatch(  Path( zip(xx_edge, yy_edge) ) ) #, facecolor='none')

    #ax.add_patch( Polygon( zip(xx_edge, yy_edge), closed=True, facecolor='red') )


    #p.set_clip_on(True)
    #p.set_clip_path( edge_path )
    #polys[0].set_clip_path( edge_path )



    return ax


def draw_dipole_field_line(
        ax,
        L,
        phi = pi/2.,
        tinc = 0.0,
        rmin = 1.0,
        fmt={'color':'b','linestyle':'solid',},
        plot=True,
        ):

    the1 =  0.0
    the2 =  pi

    for (the1, the2) in zip([0.0, 0.0], [pi, -pi]):

        #front
        xx = []
        yy = []
        for theta in np.linspace(the1, the2, 100):
            r1 = L * np.sin(theta + tinc)**2.
            if r1 > rmin:
                xx.append(y(phi,theta)*r1)
                yy.append(z(phi,theta)*r1)

        if plot:
            ax.plot(xx, yy, **fmt)


    return ax


def draw_reconnecting_field_line(
        ax,
        L,
        RLC,
        phi = pi/2.,
        tinc = 0.0,
        rmin = 1.0,
        fmt={'color':'b','linestyle':'solid',},
        plot=True,
        sheet_len = 1.0,
        ):

    sheet_len -= 1.0 #correct length of reconnecting sheet


    #this is the angle where dipole field is parallel to equator
    magic_angle1 = np.deg2rad(90.0-35.3) #top part
    magic_angle2 = np.deg2rad(35.3)+pi/2. #bottom part

    for (the1, the2, dire) in zip(
            [0.0,     0.0,  pi,    -pi/2.],
            [pi/2., -pi/2., pi/2., -pi],
            [+1,       +1,  -1,     -1]
            ):

        #front
        xx = []
        yy = []
        trans = 0.0
        for theta in np.linspace(the1, the2, 100):
            r1 = L * np.sin(theta + tinc)**2.

            #top lines
            if (dire > 0):
                trans = 0.0
                if (abs(theta) > magic_angle1 ):
                    trans = (abs(theta) - magic_angle1)/(pi/2.0 - magic_angle1)
                    #trans = trans**3.0 #smoothen transition

            #bottom lines
            if (dire < 0):
                trans = 1.0
                if (abs(theta) < magic_angle2):
                    trans = -(abs(theta) - magic_angle2)/(magic_angle2-pi/2.)
                    #trans = trans**3.0 #smoothen transition
                else:
                    trans = 0.0
            trans = trans**3.0 #sharper transition


            if r1 > rmin:
                xx1 = y(phi,theta)*r1
                yy1 = z(phi,theta)*r1

                #ref_theta = np.sign(theta)*pi/2. #equator
                ref_theta = np.sign(theta)*pi/2. - tinc #oblique
                xx2 = y(phi, ref_theta)*RLC*(1.0 + sheet_len*trans)
                yy2 = z(phi, ref_theta)*RLC*(1.0 + sheet_len*trans)

                xxi = (1.0-trans)*xx1 + trans*xx2
                yyi = (1.0-trans)*yy1 + trans*yy2

                xx.append(xxi)
                yy.append(yyi)

        if plot:
            ax.plot(xx, yy, **fmt)

    return ax



def draw_open_field_line(
        ax,
        L,
        RLC,
        phi = pi/2.,
        tinc = 0.0,
        rmin = 1.0,
        fmt={'color':'b','linestyle':'solid',},
        plot=True,
        ):

    #this is the angle where dipole field is parallel to equator
    magic_angle1 = np.deg2rad(90.0-35.3) #top part
    magic_angle2 = np.deg2rad(35.3)+pi/2. #bottom part

    #for (the1, the2, dire) in zip(
    #        [0.0,     0.0,  pi,    -pi/2.],
    #        [pi/2., -pi/2., pi/2., -pi],
    #        [+1,       +1,  -1,     -1]
    #        ):
    for (the1, the2, dire) in zip(
            [0.0-tinc, -pi/2.      ],
            [pi/2.   , -pi-tinc    ],
            [+1      , -1          ]
            ):

        #front
        xx = []
        yy = []
        trans = 0.0
        for theta in np.linspace(the1, the2, 100):
            r1 = L * np.sin(theta + tinc)**2.

            #top lines
            if (dire > 0):
                trans = 0.0
                if (abs(theta+tinc) > magic_angle1 ):
                    trans = (abs(theta+tinc) - magic_angle1)/(pi/2.0 - magic_angle1)

            #bottom lines
            if (dire < 0):
                trans = 1.0
                if (abs(theta+tinc) < magic_angle2):
                    trans = -(abs(theta+tinc) - magic_angle2)/(magic_angle2-pi/2.)
                else:
                    trans = 0.0

            #trans = np.sqrt(trans) #sharper transition
            trans = trans**1.

            if r1 > rmin:
                xx1 = y(phi,theta)*r1
                yy1 = z(phi,theta)*r1

                #if trans > 0.0:
                #    trans = 1.0

                if (dire > 0):
                    xx2 = y(phi, np.sign(theta)*magic_angle1)*L*(1.0 + 2.*trans)
                    yy2 = z(phi, np.sign(theta)*magic_angle1)*L*(1.0 + 2.*trans)
                if (dire < 0):
                    xx2 = y(phi, np.sign(theta)*magic_angle2)*L*(1.0 + 2.*trans)
                    yy2 = z(phi, np.sign(theta)*magic_angle2)*L*(1.0 + 2.*trans)

                xxi = (1.0-trans)*xx1 + trans*xx2
                yyi = (1.0-trans)*yy1 + trans*yy2

                xx.append(xxi)
                yy.append(yyi)

        if plot:
            ax.plot(xx, yy, **fmt)

    return ax

# https://stackoverflow.com/questions/6802577/rotation-of-3d-vector
def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


def draw_tilted_plane(
        ax,
        phi_start = 0.0,
        phi_stop = 2*pi,
        rmin = 1.0,
        rmax = 1.0,
        tinc = 0.0,
        phic = 0.0,
        fmt={'color':'b','linestyle':'solid',},
        ):

    xxi = y(phic, pi/2.)*rmin
    yyi = z(phic, pi/2.)*rmin
    zzi = x(phic, pi/2.)*rmin
    rotaxis = [xxi, yyi, zzi]
    rotmat = rotation_matrix(rotaxis, tinc)

    #wheels
    theta = pi/2.
    for rfac in np.linspace(rmin, rmax, 5):
        xx = []
        yy = []
        for phi in np.linspace(phi_start,phi_stop, 200):
            xxi = y(phi,theta)*rfac
            yyi = z(phi,theta)*rfac
            zzi = x(phi,theta)*rfac

            rotp = np.dot(rotmat, [xxi,yyi,zzi])
            xx.append(rotp[0])
            yy.append(rotp[1])
        if plot:
            ax.plot(xx, yy, **fmt)


    #spokes
    for phi in np.linspace(phi_start,phi_stop, 40):
        xx = []
        yy = []
        for rfac in np.linspace(rmin, rmax, 20):
            xxi = y(phi,theta)*rfac
            yyi = z(phi,theta)*rfac
            zzi = x(phi,theta)*rfac

            rotp = np.dot(rotmat, [xxi,yyi,zzi])
            xx.append(rotp[0])
            yy.append(rotp[1])
        if plot:
            ax.plot(xx, yy, **fmt)



    return ax


