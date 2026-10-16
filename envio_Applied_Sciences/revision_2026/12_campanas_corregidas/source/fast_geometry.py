"""Compiled double-precision geometric kernels; no fastmath or parallel reduction.

All branch rules match historical_engine. corrected_slider switches on the
documented rigid-rod repair; geared ratio projection is done by revision_engine.
"""
import math
import numpy as np
from numba import njit

@njit(cache=True)
def intersection(cx,cy,r, gx,gy,s, px,py,first,positive):
    dx=gx-cx;dy=gy-cy;d=math.hypot(dx,dy)
    if d<1e-12 or d>r+s+1e-9 or d<abs(r-s)-1e-9:
        return False,0.,0.
    a=(r*r-s*s+d*d)/(2*d);h2=r*r-a*a
    if h2 < -1e-9:return False,0.,0.
    h=math.sqrt(max(h2,0.));ex=dx/d;ey=dy/d
    ax=cx+a*ex-h*ey;ay=cy+a*ey+h*ex
    bx=cx+a*ex+h*ey;by=cy+a*ey-h*ex
    if first:
        pick=(ay>=by) if positive else True
    else:pick=(ax-px)**2+(ay-py)**2 <= (bx-px)**2+(by-py)**2
    return (True,ax,ay) if pick else (True,bx,by)

@njit(cache=True)
def simulate(fam,x,theta,corrected_slider=True):
    # A,B,C,D,P; fixed anchors are reconstructed from parameters by validator.
    joints=np.zeros((len(theta),5,2));margin=1e100
    pbx=0.;pby=0.;pdx=0.;pdy=0.
    if fam==1 and x[1]<=x[0]:return False,joints,0.
    for i in range(len(theta)):
        th=theta[i]
        if fam==1:
            r,l,q,h,x0=x
            ax=r*math.cos(th)+(x0 if corrected_slider else 0.)
            ay=r*math.sin(th)
            rad=l*l-ay*ay
            if rad<0:return False,joints,0.
            bx=ax+math.sqrt(max(rad,0.))+(0. if corrected_slider else x0);by=0.
            ux=bx-ax;uy=-ay;ln=math.hypot(ux,uy)
            ux/=ln;uy/=ln
            cx=ax;cy=ay;ddx=bx;ddy=by
            xx=ax+q*l*ux-h*uy;yy=ay+q*l*uy+h*ux
            delta=l-r
        elif fam==7:
            L0,L2,L4,L3,L5,ratio,phase,rp,phip=x
            ax=L2*math.cos(th);ay=L2*math.sin(th)
            bx=L0+L4*math.cos(ratio*th+phase);by=L4*math.sin(ratio*th+phase)
            ok,cx,cy=intersection(ax,ay,L3,bx,by,L5,pdx,pdy,i==0,True)
            if not ok:return False,joints,0.
            ux=cx-ax;uy=cy-ay;ln=math.hypot(ux,uy);ux/=ln;uy/=ln
            xx=cx+rp*(math.cos(phip)*ux-math.sin(phip)*uy)
            yy=cy+rp*(math.cos(phip)*uy+math.sin(phip)*ux)
            ddx=cx;ddy=cy
            delta=min(abs((cx-ax)*(by-cy)-(cy-ay)*(bx-cx)),abs(ax*by-ay*(bx-L0)))
        else:
            L0,L2,L3,L4=x[:4]
            ax=L2*math.cos(th);ay=L2*math.sin(th)
            ok,bx,by=intersection(ax,ay,L3,L0,0.,L4,pbx,pby,i==0,True)
            if not ok:return False,joints,0.
            delta=abs((bx-ax)*(-by)-(by-ay)*(L0-bx))
            ux=bx-ax;uy=by-ay;ln=math.hypot(ux,uy);ux/=ln;uy/=ln
            if fam==0:
                rp,phi=x[4:6]
                xx=ax+rp*(math.cos(phi)*ux-math.sin(phi)*uy)
                yy=ay+rp*(math.cos(phi)*uy+math.sin(phi)*ux)
                cx=ax;cy=ay;ddx=bx;ddy=by
            else:
                if fam in (2,3,6):
                    rc,phi=x[4:6]
                    if fam==3:
                        ux=L0-bx;uy=-by;ln=math.hypot(ux,uy);ux/=ln;uy/=ln
                        cx=bx+rc*(math.cos(phi)*ux-math.sin(phi)*uy)
                        cy=by+rc*(math.cos(phi)*uy+math.sin(phi)*ux)
                    else:
                        cx=ax+rc*(math.cos(phi)*ux-math.sin(phi)*uy)
                        cy=ay+rc*(math.cos(phi)*uy+math.sin(phi)*ux)
                    if fam==6:
                        gx=L0;gy=0.;L5,L6,pu,pv=x[6:10];pos=False
                    else:
                        gx,gy,L5,L6,pu,pv=x[6:12];pos=False
                else:
                    L5,L6,pu,pv=x[4:8]
                    if fam==4:cx=ax;cy=ay;gx=L0;gy=0.;pos=False
                    else:cx=bx;cy=by;gx=0.;gy=0.;pos=True
                ok,ddx,ddy=intersection(cx,cy,L5,gx,gy,L6,pdx,pdy,i==0,pos)
                if not ok:return False,joints,0.
                ux=ddx-cx;uy=ddy-cy;ln=math.hypot(ux,uy);ux/=ln;uy/=ln
                xx=cx+pu*ux-pv*uy;yy=cy+pu*uy+pv*ux
                delta=min(delta,abs((ddx-cx)*(gy-ddy)-(ddy-cy)*(gx-ddx)))
        joints[i,0,0]=ax;joints[i,0,1]=ay
        joints[i,1,0]=bx;joints[i,1,1]=by
        joints[i,2,0]=cx;joints[i,2,1]=cy
        joints[i,3,0]=ddx;joints[i,3,1]=ddy
        joints[i,4,0]=xx;joints[i,4,1]=yy
        pbx=bx;pby=by;pdx=ddx;pdy=ddy;margin=min(margin,delta)
    return True,joints,margin

@njit(cache=True)
def resample(P,n):
    m=len(P)
    if math.hypot(P[0,0]-P[-1,0],P[0,1]-P[-1,1])<1e-9:m-=1
    s=np.zeros(m+1)
    for i in range(m):
        j=(i+1)%m;s[i+1]=s[i]+math.hypot(P[j,0]-P[i,0],P[j,1]-P[i,1])
    ans=np.zeros((n,2));j=0
    for i in range(n):
        t=i*s[m]/n
        while j<m-1 and s[j+1]<t:j+=1
        w=(t-s[j])/max(s[j+1]-s[j],1e-300);k=(j+1)%m
        ans[i]=P[j]*(1-w)+P[k]*w
    return ans

@njit(cache=True)
def aligned_metrics(P,Y,shift):
    n=len(P);pm=np.zeros(2);ym=np.zeros(2)
    for i in range(n):pm+=P[i];ym+=Y[i]
    pm/=n;ym/=n;a=0.;b=0.
    for i in range(n):
        k=(i-shift)%n
        p0=P[i,0]-pm[0];p1=P[i,1]-pm[1]
        y0=Y[k,0]-ym[0];y1=Y[k,1]-ym[1]
        a+=p0*y0+p1*y1;b+=p0*y1-p1*y0
    angle=math.atan2(b,a);c=math.cos(angle);s=math.sin(angle)
    aligned=np.empty_like(P);ss=0.;mx=0.
    for i in range(n):
        px=P[i,0]-pm[0];py=P[i,1]-pm[1]
        aligned[i,0]=c*px-s*py+ym[0];aligned[i,1]=s*px+c*py+ym[1]
        j=(i-shift)%n;e=(aligned[i,0]-Y[j,0])**2+(aligned[i,1]-Y[j,1])**2
        ss+=e;mx=max(mx,e)
    return math.sqrt(ss/n),math.sqrt(mx),aligned

@njit(cache=True)
def phase_alignment(P,Y,stride=4,allow=True):
    n=len(P);best=1e100;bs=0
    for s in range(0,n,stride if allow else n):
        rms,mx,ap=aligned_metrics(P,Y,s)
        if rms<best:best=rms;bs=s
    coarse=bs;best=1e100
    mask=np.zeros(n,np.bool_)
    if allow:
        for ds in range(-max(1,stride),max(1,stride)+1):mask[(coarse+ds)%n]=True
    else:mask[0]=True
    bm=0.;bp=P.copy()
    for s in range(n):
        if mask[s]:
            rms,mx,ap=aligned_metrics(P,Y,s)
            if rms<best:best=rms;bs=s;bm=mx;bp=ap
    return best,bm,bs,bp
