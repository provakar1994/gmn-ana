import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def fit_fn(x, a0, a1, a2):
    return a0 + a1*(x-a2)*(x-a2)

data = pd.read_csv("chi2_vs_R_sbs50_.82T_ap.csv")
R = data['R']
chi2 = data['chi2']

param, covar = curve_fit(fit_fn, R, chi2)
chi2m = param[0]
error = param[1]
Rm = param[2]
print("Rm = {}, chi2m = {}, Err = {}".format(Rm, chi2m, 1./np.sqrt(error))) 

fit_y = fit_fn(R, chi2m, error,  Rm)
plt.plot(R, chi2, 'o', label='data')
plt.plot(R, fit_y, '-', label='fit: {} + {} * (R-{})^2'.format(round(chi2m,3),round(error,3),round(Rm,3)))

# plt.scatter(R, chi2)
plt.xlabel('R')
plt.ylabel(r'$\chi^2$')
plt.legend()
plt.grid()
plt.show()
