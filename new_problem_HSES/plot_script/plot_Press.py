import numpy as np
import matplotlib.pyplot as plt

# dimensionless parameters: set G = M = a = 1
G = 1.0
M = 1.0
a = 1.0

def P_hernquist(r, G=1.0, M=1.0, a=1.0):
    ap = r + a
    term_log = np.log(ap / r)
    term1 = a / ap
    term2 = a**2 / (2.0 * ap**2)
    term3 = a**3 / (3.0 * ap**3)
    term4 = a**4 / (4.0 * ap**4)
    return (G * M**2) / (2.0 * np.pi * a**4) * (term_log - term1 - term2 - term3 - term4)

r = np.logspace(-2, 1, 400)  # 0.01 to 10 (create 400 points between 0.01 to 10)
P = P_hernquist(r, G, M, a)

plt.figure()
plt.loglog(r, P)
plt.xlabel(r'$r/a$')
plt.ylabel(r'$P(r)\;[\;G=M=a=1\;]$')
plt.title('Pressure profile for Hernquist sphere (dimensionless)')
plt.grid(True, which='major')
plt.show()
