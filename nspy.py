import numpy as np
import math
from pylab import *
import matplotlib.patches as patches
from scipy.optimize import minimize_scalar

from matplotlib.path import Path
from matplotlib.patches import PathPatch

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


    for rfac in np.linspace(r_start, r_stop, 100):
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


def draw_phi(ax, phi1, phi2, theta, r1, r2, fmt, fmtb):

    #front
    xx = []
    yy = []
    for phi in np.linspace(phi1, phi2, 100):
        xx.append(y(phi,theta)*r1)
        yy.append(z(phi,theta)*r1)
    ax.plot(xx, yy, **fmt)

    #back
    xxB = []
    yyB = []
    for phi in np.linspace(phi1, phi2, 100):
        xxB.append(y(phi,theta)*r2)
        yyB.append(z(phi,theta)*r2)
    ax.plot(xxB, yyB, **fmtb)

    return ax, xx, yy

def draw_theta(ax, the1, the2, phi, r1, r2, fmt, fmtb):

    #front
    xx = []
    yy = []
    for theta in np.linspace(the1, the2, 100):
        xx.append(y(phi,theta)*r1)
        yy.append(z(phi,theta)*r1)
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
        yyB.append(z(phi,theta)*r2)
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
        yy.append(z(phi,theta)*r1)
        
    if plot:
        ax.plot(xx, yy, **fmt)

    #back
    xxB = []
    yyB = []
    for (phi, theta) in zip(np.linspace(phi1, phi2, 100), np.linspace(the1, the2, 100)):
        xxB.append(y(phi,theta)*r2)
        yyB.append(z(phi,theta)*r2)

    if plot:
        ax.plot(xxB, yyB, **fmtb)

    return ax, xx, yy



def clip_xy(xx, yy, ylim):
    
    yc = np.empty_like(yy)
    for j, y in enumerate(yy):
        #if y > ylim:
        #    yc[j] = ylim
        #else:
        #    yc[j] = yy[j]

        if y < ylim:
            yc[j] = ylim
        else:
            yc[j] = yy[j]

    return xx, yc


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
    dthe = 0.30
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

        #for kk, ylim in enumerate(np.linspace(0.5, 1.8, 10)): #dark bottom - light top
        #for kk, ylim in enumerate(np.linspace(1.5, 0.5, 10)): 
        for kk, ylim in enumerate(np.linspace(top, bot, 20)): 
            #fmtf['facecolor'] = cmap(1.0 - kk / 19.0)
            fmtf['facecolor'] = cmap( kk / 19.0)

            xxc,yyc = clip_xy(xx,yy, ylim)
            polys.append( Polygon(zip(xxc, yyc), closed=True, **fmtf) )



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




