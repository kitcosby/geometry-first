from sage.all import *
# E77: Verify Kit's double-angle identities.
# ID1: sin(2x)/4 =?= 1/(urx + uxp)
# ID2: cos(2x) =?= (srx + crx - sxp - cxp)/(urx + uxp)

import math

def ops(deg):
    r=math.radians(deg); s=math.sin(r); c=math.cos(r)
    srx=abs(1/c)+abs(c/s) if False else (1+c)/s  # QI: csc+cot
    # Actually: srx = |csc x| + cot x. In QI, csc>0, so = csc+cot = (1+cos)/sin.
    srx = (1+c)/s
    sxp = 1/srx
    # cxp = |sec x| + tan x. In QI = sec+tan = (1+sin)/cos.
    cxp = (1+s)/c
    crx = 1/cxp
    urx = srx-crx; uxp = cxp-sxp
    return srx,sxp,cxp,crx,urx,uxp

print("Verifying ID1: sin(2x)/4 =?= 1/(urx+uxp)")
print("Verifying ID2: cos(2x) =?= (srx+crx-sxp-cxp)/(urx+uxp)")
print()
max_err1=0; max_err2=0
for deg in [5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85]:
    r=math.radians(deg)
    srx,sxp,cxp,crx,urx,uxp=ops(deg)
    # ID1
    lhs1=math.sin(2*r)/4
    rhs1=1/(urx+uxp)
    err1=abs(lhs1-rhs1)
    max_err1=max(max_err1,err1)
    # ID2
    lhs2=math.cos(2*r)
    rhs2=(srx+crx-sxp-cxp)/(urx+uxp)
    err2=abs(lhs2-rhs2)
    max_err2=max(max_err2,err2)
    if err1>1e-10 or err2>1e-10:
        print(f"  x={deg}: ID1 err={err1:.2e}, ID2 err={err2:.2e} <-- LARGE")

print(f"  Max ID1 error: {max_err1:.2e}")
print(f"  Max ID2 error: {max_err2:.2e}")
if max_err1<1e-10 and max_err2<1e-10:
    print("  BOTH IDENTITIES CHECKED (numeric, 17 points in QI).")
print()
# Also check: sin^2(2x) + cos^2(2x) = 1 via these forms.
print("Consistency: sin^2+cos^2 via operator forms:")
for deg in [15,30,45,60]:
    r=math.radians(deg)
    srx,sxp,cxp,crx,urx,uxp=ops(deg)
    s2x_form = 1/(urx+uxp)  # this is sin(2x)/4
    c2x_form = (srx+crx-sxp-cxp)/(urx+uxp)
    # sin(2x) = 4*s2x_form, cos(2x) = c2x_form
    check = (4*s2x_form)**2 + c2x_form**2
    print(f"  x={deg}: (4*sin_form)^2 + cos_form^2 = {check:.10f} (expect 1)")
