import numpy as np
import math
from pylab import *
import matplotlib.patches as patches
from scipy.optimize import minimize_scalar


#variables
#incl = pi/8
incl = pi/2.5
req = 1
u = 0.0

rg = u*req
muc = -rg/(req-rg)
o21 = 0.3


def R(theta):
    return 1.0

def R(theta):
    return 1.0*(1-o21*cos(theta)**2)

def dR(theta):
    return 2.0*o21*sin(theta)*cos(theta)

def fa(theta):
    return dR(theta)/R(theta)

def cosg(theta):
    return (1/sqrt(1-u))/sqrt(1+fa(theta)**2)

def sing(theta):
    return fa(theta)*cosg(theta)

def mu(phi, theta):
    return cos(incl)*cos(theta)+cos(phi)*sin(incl)*sin(theta)

def bimp(phi, theta):
    return sqrt(req*(req+rg+req*mu(phi,theta)-rg*mu(phi,theta))/(1+mu(phi,theta)))

def x(phi, theta):
    return bimp(phi,theta)*R(theta)*(cos(incl)*cos(theta)+cos(phi)*sin(incl)*sin(theta))

def y(phi, theta):
    return bimp(phi,theta)*R(theta)*(sin(theta)*sin(phi))

def z(phi, theta):
    return  bimp(phi,theta)*R(theta)*(cos(theta)*sin(incl)-cos(incl)*cos(phi)*sin(theta))


def xfull(phi,theta, a,b,c):
    return bimp(phi, theta)*R(theta)* (cos(incl)*(a*cos(theta)-c*sin(theta)) + sin(incl)*(c*cos(theta)*cos(phi) + a*cos(phi)*sin(theta) - b*sin(phi)))

def yfull(phi, theta, a,b,c):
    return bimp(phi, theta)*R(theta)* (b*cos(phi) + (c*cos(theta) + a*sin(theta))*sin(phi))

def zfull(phi, theta, a,b,c):
    return bimp(phi, theta)*R(theta)* (sin(incl)*(a*cos(theta)-c*sin(theta)) + cos(incl)*(-cos(phi)*(c*cos(theta) + a*sin(theta)) + b*sin(phi)))



def draw_longitude(ax,
                   phi,
                   start=0,
                   stop=pi/2,
                   rfac=1.0,
                   backside=False,
                   fmt={'color':'k','linestyle':'solid',},
                   ):

    xx = []
    yy = []

    for theta in np.linspace(start, stop, 200):
        if not(backside):
            if mu(phi,theta) >= muc:
                xx.append(y(phi,theta)*rfac)
                yy.append(z(phi,theta)*rfac)
            else:
                lines, = ax.plot(xx, yy)
                lines.set(**fmt)
                xx = []
                yy = []
        else:
            if mu(phi,theta) <= muc:
                xx.append(y(phi,theta)*rfac)
                yy.append(z(phi,theta)*rfac)
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
                   ):

    xx = []
    yy = []

    for phi in np.linspace(start, stop, 200):
        if not(backside):
            if mu(phi,theta) >= muc:
                xx.append(y(phi,theta)*rfac)
                yy.append(z(phi,theta)*rfac)
            else:
                lines, = ax.plot(xx, yy)
                lines.set(**fmt)
                xx = []
                yy = []
        else:
            if mu(phi,theta) <= muc:
                xx.append(y(phi,theta)*rfac)
                yy.append(z(phi,theta)*rfac)
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
              fmt= {'color':'k', 
              'linestyle':'solid', 
              'lw':1.5, 
              'head_width': 0.04, 
              'head_length': 0.08,},
              rfac=1.0):

    ax.add_patch(patches.FancyArrow(
              startp[0],startp[1],
              y(phi, theta)*rfac, z(phi, theta)*rfac,
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
    ax.add_patch(Polygon(zip(xx,yy),
                      closed=True, 
                      #facecolor='black',alpha=0.5))
                      **fmt))

    return ax



