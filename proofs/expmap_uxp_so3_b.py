import numpy as np
from scipy.linalg import expm
rt3 = np.sqrt(3)
w = np.array([2-rt3, rt3, -2+2*rt3])  # (sxp, cxp, uxp) at 30deg
Y = np.array([[0, -w[2], w[1]],[w[2], 0, -w[0]],[-w[1], w[0], 0]])
th = np.linalg.norm(w)
R_rod = np.eye(3) + (np.sin(th)/th)*Y + ((1-np.cos(th))/th**2)*(Y@Y)
R_exp = expm(Y)
print('theta =', th, ' axis =', w/th)
print('Rodrigues == expm:', np.allclose(R_rod, R_exp))
print('R orthogonal:', np.allclose(R_rod@R_rod.T, np.eye(3)), '| det R =', round(np.linalg.det(R_rod), 12))
print('tr R = 1 + 2cos(theta):', np.allclose(np.trace(R_rod), 1 + 2*np.cos(th)))
print('Y^3 = -theta^2 Y (so(3) identity):', np.allclose(Y@Y@Y, -th**2*Y))
np.set_printoptions(precision=4, suppress=True)
print('exp(Y) at 30deg =\n', R_rod)
