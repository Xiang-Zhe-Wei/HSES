import numpy as np
import matplotlib.pyplot as plt

# dimensionless parameters: set M = a = 1
M = 1.0
a = 1.0

def rho_hernquist(r, M=1.0, a=1.0):
    return M * a / (2.0 * np.pi * r * (r + a)**3)


r = np.logspace(-2, 1, 400)  # 0.01 to 10 (create 400 points between 0.01 to 10)
rho = rho_hernquist(r, M, a)

plt.figure()
plt.loglog(r, rho)
plt.xlabel(r'$r/a$')
plt.ylabel(r'$\rho(r)\;[\;M=a=1\;]$')
plt.title('Density profile for Hernquist sphere (dimensionless)')
plt.grid(True, which='major')
plt.show()
