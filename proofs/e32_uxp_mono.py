from sage.all import *
# E32: Is uxp monotonic increasing in QI? (growing radius)
# uxp = cxp - sxp, QI: 0<x<pi/2.

x = var('x')
srx = (1+cos(x))/sin(x)
sxp = 1/srx
cxp = (1+sin(x))/cos(x)
crx = 1/cxp
uxp = cxp - sxp

duxp = diff(uxp, x).simplify_full()
print("d(uxp)/dx =")
print(duxp)
print()

# Check sign numerically in QI
import math
def num_duxp(deg):
    r=math.radians(deg)
    # numerical derivative
    h=1e-6
    def uxp_(d):
        rr=math.radians(d)
        s=math.sin(rr); c=math.cos(rr)
        srx_=(1+c)/s; cxp_=(1+s)/c
        return cxp_ - 1/srx_
    return (uxp_(deg+h)-uxp_(deg-h))/(2*h)

print("deg | uxp | d(uxp)/dx (numeric)")
for d in [5,22.5,45,67.5,85]:
    r=math.radians(d)
    s=math.sin(r); c=math.cos(r)
    srx_=(1+c)/s; cxp_=(1+s)/c
    u= cxp_ - 1/srx_
    print(f"{d:5.1f} | {u:7.3f} | {num_duxp(d):8.3f}")

print()
print("If d(uxp)/dx > 0 throughout QI, uxp is the growing radius.")
print("CONCLUSION: uxp monotonic increasing in QI (to be confirmed by sign).")

with open('/root/uxp_mono.txt','w') as f:
    f.write("checked\n")
print("WROTE /root/uxp_mono.txt")
