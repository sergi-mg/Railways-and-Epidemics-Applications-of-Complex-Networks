# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 08:55:50 2026

@author: ASUS
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit




def read_plus_stats(lambda_lists,N_sim,N_list):
    
    data = []
    directory="../data_results/"

    z=1.96
    
    for N_i in range(len(N_list)):
        N=N_list[N_i]
        for lam in lambda_lists[N_i]:
            name="../data_results/N_"+str(N)+"_gamma_"+str(round(gamma,1))\
            +"_lambda_"+str(round(lam,4))+"_Nsim_"+str(2000)+".dat"
            filename = directory+name
                
            #read the file
            try:
                arr = np.loadtxt(filename)
            except OSError:
                print(f"Missing file:",round(lam,3),N,gamma)
                continue

            tau = arr[:, 0]
            end = arr[:, 1]

            # =========================
            # tau only when endemic=0
            # =========================
            tau0 = tau[end == 0]
            N0 = len(tau0)

            if N0 > 1:
                tau_mean = np.mean(tau0)
                tau_std = np.std(tau0, ddof=1)
                d_tau = z * tau_std / np.sqrt(N0)
                tau2_mean = np.mean(tau0**2)
                tau2_std = np.std(tau0**2, ddof=1)
                d_tau2 = z * tau2_std / np.sqrt(N0)
            else:
                tau_mean = np.nan
                d_tau = np.nan
                tau2_mean = np.nan
                d_tau2 = np.nan

            # =========================
            # P_end
            # =========================
            p_end = np.mean(end)
            p_end_std = np.std(end, ddof=1)
            Ntot = len(end)

            dP_end = z * p_end_std / np.sqrt(Ntot)

            data.append({
                "N": N,
                "lambda": lam,

                "tau_mean": tau_mean,
                "d_tau": d_tau,

                "P_end": p_end,
                "dP_end": dP_end,
                
                "tau2_mean": tau2_mean,
                "d_tau2": d_tau2,
            })
            
    df = pd.DataFrame(data)
    df = df.sort_values(["N", "lambda"])
    
    return df

def parabola(x, a, b, c):
    return a*x**2 + b*x + c

def quadratic_peak(x, y, window=5):
    """
    Fits a quadratic function around the maximum of the data and returns
    the vertex coordinates together with their 1σ uncertainties.

    Parameters
    ----------
    x, y : array_like
        Experimental data.
    window : int, optional
        Number of points taken on each side of the approximate maximum.

    Returns
    -------
    x_peak : float
        x-coordinate of the fitted maximum.
    y_peak : float
        y-coordinate of the fitted maximum.
    dx_peak : float
        1σ uncertainty of x_peak.
    dy_peak : float
        1σ uncertainty of y_peak.
    popt : ndarray
        Fitted parameters (a, b, c).
    pcov : ndarray
        Covariance matrix of the fitted parameters.
    """

    # Approximate maximum
    i = np.argmax(y)

    # Fitting window
    i0 = max(0, i - window)
    i1 = min(len(x), i + window + 1)

    x_fit = x[i0:i1]
    y_fit = y[i0:i1]

    # Quadratic fit
    popt, pcov = curve_fit(parabola, x_fit, y_fit)
    a, b, c = popt

    # Vertex coordinates
    x_peak = -b / (2 * a)
    y_peak = c - b**2 / (4 * a)

    # Error propagation
    grad_x = np.array([
        b / (2 * a**2),
        -1 / (2 * a),
        0
    ])

    grad_y = np.array([
        b**2 / (4 * a**2),
        -b / (2 * a),
        1
    ])

    dx_peak = np.sqrt(grad_x @ pcov @ grad_x)
    dy_peak = np.sqrt(grad_y @ pcov @ grad_y)

    return x_peak, y_peak, dx_peak, dy_peak

def calc_reg(x,y):
    import numpy as np
    #y=mx+b
    N=np.size(x)
    #calcul mitjanes
    xmed=np.sum(x)/N
    x2med=np.sum(x**2)/N
    ymed=np.sum(y)/N
    y2med=np.sum(y**2)/N
    xymed=np.sum(x*y)/N
    #calcul sigmes
    sigx2=x2med-xmed**2
    sigy2=y2med-ymed**2
    sigxy=xymed-xmed*ymed
    #calcul regressio
    m=sigxy/sigx2
    b=ymed-m*xmed
    r=sigxy/((sigx2*sigy2)**0.5)
    #calcul incerteses
    dyreg=((sigy2*(1-r**2)*N)/(N-2))**0.5
    dm=dyreg/(N*sigx2)**0.5
    db=dyreg*(x2med/(N*sigx2))**0.5
    #retornem els valors
    return m,dm,b,db,r
 

#%%
N_list = [10000, 30000, 50000,100000, 300000, 500000, 1000000]
#gamma=2.5
gamma=2.5

#list general
lambda_min=0.
lambda_max=0.1
d_lambda=0.01
lambda_list_g=np.arange(lambda_min,lambda_max+0.001,d_lambda)

#list 1
lambda_min_i=np.array([0.02,0.02,0.01,0.01,0.01,0.00,0.00])
lambda_max_i=lambda_min_i+0.01
lambda_max_i[0]=0.05
d_lambda_i=0.001

#list 2
lambda_min_2=np.array([0.029,0.018,0.016,0.012,0.010,0.008,0.007])
lambda_max_2=np.array([0.035,0.026,0.020,0.017,0.013,0.010,0.009])
d_lambda_2=0.0002

lambda_lists = []

for N_i in range(len(N_list)):
    lambda_list_i = np.arange(lambda_min_i[N_i],
                              lambda_max_i[N_i] + 0.0001,
                              d_lambda_i)

    lambda_list_2 = np.arange(lambda_min_2[N_i],
                              lambda_max_2[N_i] + 0.0001,
                              d_lambda_2)

    lambda_list = np.concatenate((lambda_list_g,
                                  lambda_list_i,
                                  lambda_list_2))
    
    lambda_list = np.unique(lambda_list)

    lambda_lists.append(lambda_list)

Nsim = 2000

#%%
data_df_1=read_plus_stats(lambda_lists, Nsim, N_list)



#%%
#gamma=3.5
gamma=3.5
lambda_min=0.05
lambda_max=0.17
d_lambda=0.01

lambda_list_g=np.arange(lambda_min,lambda_max+0.001,d_lambda)
lambda_min_i=np.array([0.13,0.12,0.12,0.11,0.11,0.10,0.10])
lambda_max_i=lambda_min_i+0.01
lambda_max_i[0]=0.15
lambda_min_2=np.array([0.138,0.127,0.120,0.117,0.108,0.106,0.101])
lambda_max_2=np.array([0.16,0.14,0.130,0.130,0.115,0.110,0.105])

lambda_lists = []

for N_i in range(len(N_list)):
    lambda_list_i = np.arange(lambda_min_i[N_i],
                              lambda_max_i[N_i] + 0.0001,
                              d_lambda_i)

    lambda_list_2 = np.arange(lambda_min_2[N_i],
                              lambda_max_2[N_i] + 0.0001,
                              d_lambda_2)

    lambda_list = np.concatenate((lambda_list_g,
                                  lambda_list_i,
                                  lambda_list_2))
    
    lambda_list = np.unique(lambda_list)

    lambda_lists.append(lambda_list)

Nsim = 2000

#%%
data_df_2=read_plus_stats(lambda_lists, Nsim, N_list)

#%%

fig, ax = plt.subplots(3, 2, figsize=(8,15))

plt.subplots_adjust(
    left=0.10,
    right=0.98,
    top=0.88,
    bottom=0.12,
    hspace=0.30,
    wspace=0.30
)

handles = []
labels = [r"$N=10^4$", r"$N=3\times10^4$", r"$N=5\times10^4$",
          r"$N=10^5$", r"$N=3\times10^5$", r"$N=5\times10^5$",
          r"$N=10^6$"]

cmap = plt.get_cmap("viridis")

# =====================================================
# γ = 2.5 (columna izquierda)
# =====================================================

count = 0
for N, dfN in data_df_1.groupby("N"):

    dfN = dfN.sort_values("lambda")

    line, = ax[0,0].plot(dfN["lambda"], dfN["tau_mean"],
                         "o--", markersize=2, color=cmap(count/7))
    ax[0,0].fill_between(
        dfN["lambda"],
        dfN["tau_mean"]-dfN["d_tau"],
        dfN["tau_mean"]+dfN["d_tau"],
        alpha=0.2,
        color=cmap(count/7)
    )

    handles.append(line)

    ax[1,0].plot(dfN["lambda"], dfN["tau2_mean"],
                 "o--", markersize=2, color=cmap(count/7))
    ax[1,0].fill_between(
        dfN["lambda"],
        dfN["tau2_mean"]-dfN["d_tau2"],
        dfN["tau2_mean"]+dfN["d_tau2"],
        alpha=0.2,
        color=cmap(count/7)
    )

    ax[2,0].plot(dfN["lambda"], dfN["P_end"],
                 "o--", markersize=2, color=cmap(count/7))
    ax[2,0].fill_between(
        dfN["lambda"],
        dfN["P_end"]-dfN["dP_end"],
        dfN["P_end"]+dfN["dP_end"],
        alpha=0.2,
        color=cmap(count/7)
    )

    count += 1


# =====================================================
# γ = 3.5 (columna derecha)
# =====================================================

count = 0
for N, dfN in data_df_2.groupby("N"):

    dfN = dfN.sort_values("lambda")

    ax[0,1].plot(dfN["lambda"], dfN["tau_mean"],
                 "o--", markersize=2, color=cmap(count/7))
    ax[0,1].fill_between(
        dfN["lambda"],
        dfN["tau_mean"]-dfN["d_tau"],
        dfN["tau_mean"]+dfN["d_tau"],
        alpha=0.2,
        color=cmap(count/7)
    )

    ax[1,1].plot(dfN["lambda"], dfN["tau2_mean"],
                 "o--", markersize=2, color=cmap(count/7))
    ax[1,1].fill_between(
        dfN["lambda"],
        dfN["tau2_mean"]-dfN["d_tau2"],
        dfN["tau2_mean"]+dfN["d_tau2"],
        alpha=0.2,
        color=cmap(count/7)
    )

    ax[2,1].plot(dfN["lambda"], dfN["P_end"],
                 "o--", markersize=2, color=cmap(count/7))
    ax[2,1].fill_between(
        dfN["lambda"],
        dfN["P_end"]-dfN["dP_end"],
        dfN["P_end"]+dfN["dP_end"],
        alpha=0.2,
        color=cmap(count/7)
    )

    count += 1


# =====================================================
# Límites γ=2.5
# =====================================================

ax[0,0].set_xlim(0,0.05)
ax[0,0].set_ylim(1,1.8)

ax[1,0].set_xlim(0,0.05)
ax[1,0].set_ylim(0,450)

ax[2,0].set_xlim(0,0.05)
ax[2,0].set_ylim(0,0.125)


# =====================================================
# Límites γ=3.5
# =====================================================

ax[0,1].set_xlim(0.09,0.17)
ax[0,1].set_ylim(1.5,3.5)

ax[1,1].set_xlim(0.09,0.17)
ax[1,1].set_ylim(0,550)

ax[2,1].set_xlim(0.09,0.17)
ax[2,1].set_ylim(0,0.12)


# =====================================================
# Etiquetas
# =====================================================

ylabel = [
    r"$\langle\tau\rangle$",
    r"$\langle\tau^2\rangle$",
    r"$P_{\rm end}$"
]

for i in range(3):

    ax[i,0].set_ylabel(ylabel[i], fontsize=16)

    ax[i,0].tick_params(axis="both", labelsize=14)
    ax[i,1].tick_params(axis="both", labelsize=14)

    ax[i,0].grid(True, linestyle="--", alpha=0.4)
    ax[i,1].grid(True, linestyle="--", alpha=0.4)

ax[2,0].set_xlabel(r"$\lambda$", fontsize=16)
ax[2,1].set_xlabel(r"$\lambda$", fontsize=16)


# =====================================================
# Leyenda
# =====================================================

fig.legend(
    handles,
    labels,
    loc="upper center",
    bbox_to_anchor=(0.5,1.0),
    ncol=4,
    fontsize=15,
    markerscale=3
)


# =====================================================
# Subtítulos
# =====================================================

ax[2,0].set_title(r"(a) $\gamma=2.5$", y=-0.55, fontsize=16)
ax[2,1].set_title(r"(b) $\gamma=3.5$", y=-0.55, fontsize=16)

plt.savefig("../plots/plot_1_2.pdf", bbox_inches="tight")
plt.show()

#%%
#gamma=3.5
print("------------------------------------------")
print("gamma=3.5")
print("------------------------------------------")
lambda_p_35=np.zeros((len(N_list)))

lp_35=np.zeros((7,2))
tau_lp_35=np.zeros((7,2))

lp_2_35=np.zeros((7,2))
tau2_lp_35=np.zeros((7,2))

count=0
for N, dfN in data_df_2.groupby("N"):
    dfN = dfN.sort_values("lambda")
    y=dfN["tau_mean"]
    y_2=dfN["tau2_mean"]
    x=dfN["lambda"]
    lp_35[count,0],tau_lp_35[count,0],lp_35[count,1],tau_lp_35[count,1]\
        =quadratic_peak(x, y, window=2)
    lp_2_35[count,0],tau2_lp_35[count,0],lp_2_35[count,1],tau2_lp_35[count,1]\
        =quadratic_peak(x, y_2, window=5)
    count+=1
    
#%% 
#find the position of lambda_p
  
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Finite-size scaling ansatz
def lambda_peak(N, lambda_c, A, nu):
    return lambda_c + A * N**(-1/nu)

# N, lambda_p y sus errores
lambda_p_35 = lp_35[:,0]
lambda_p_err_35 = lp_35[:,1]

lambda_p_2_35 = lp_2_35[:,0]
lambda_p_2_err_35 = lp_2_35[:,1]

N=np.array(N_list)

# Initial guess
p0 = [lambda_p_35[-1], 1.0, 3.0]


# Fit
popt, pcov = curve_fit(
    lambda_peak,
    N,
    lambda_p_35,
    sigma=lambda_p_err_35,
    absolute_sigma=True,
    p0=p0
)

lambda_c_35, A_35, nu_35 = popt
dlambda_c_35, dA_35, dnu_35 = np.sqrt(np.diag(pcov))


p0 = [lambda_p_2_35[-1], 1.0, 3.0]
# Fit
popt_2, pcov_2 = curve_fit(
    lambda_peak,
    N,
    lambda_p_2_35,
    sigma=lambda_p_2_err_35,
    absolute_sigma=True,
    p0=p0
)

lambda_c_2_35, A_2_35, nu_2_35 = popt_2
dlambda_c_2_35, dA_2_35, dnu_2_35 = np.sqrt(np.diag(pcov_2))



print("tau")
print("lambda_c:", lambda_c_35, "+-", dlambda_c_35)

print("tau^2")
print("lambda_c:", lambda_c_2_35, "+-", dlambda_c_2_35)

lambda_c_35_f = (lambda_c_2_35 + lambda_c_35) / 2
d_lambda_c_35_f = 0.5 * np.sqrt(dlambda_c_2_35**2 + dlambda_c_35**2)
print("lambda_c:", lambda_c_35_f, "+-", d_lambda_c_35_f)

print("tau")
print("nu:", nu_35, "+-", dnu_35)

print("tau^2")
print("nu:", nu_2_35, "+-", dnu_2_35)

nu_35_f = (nu_2_35 + nu_35) / 2
d_nu_35_f = 0.5 * np.sqrt(dnu_2_35**2 + dnu_35**2)
print("nu:", nu_35_f, "+-", d_nu_35_f)

#%%
#find the gamma_n
# Logarithms
logN = np.log(N)
log_tau_35 = np.log(tau_lp_35[:,0])
log_tau2_35 = np.log(tau2_lp_35[:,0])

# Linear regressions: y = m x + b
m1_35, dm1_35, b1_35, db1_35, r1_35 = calc_reg(logN, log_tau_35)
m2_35, dm2_35, b2_35, db2_35, r2_35 = calc_reg(logN, log_tau2_35)
g1_nu_35=m1_35 
dg1_nu_35=dm1_35
g2_nu_35=m2_35 
dg2_nu_35=dm2_35
print("gamma_1/nu:", g1_nu_35,"+-",dg1_nu_35)
print("gamma_2/nu:", g2_nu_35,"+-",dg2_nu_35)

# Smooth curves
Nfit = np.logspace(np.log10(N.min()), np.log10(N.max()), 500)

tau_fit_35 = np.exp(b1_35) * Nfit**m1_35
tau2_fit_35 = np.exp(b2_35) * Nfit**m2_35




#%%
#P_end
from scipy.interpolate import interp1d

P_end_lc_35=np.zeros((7))
P_end_lc_2_35=np.zeros((7))
count=0
for N_i,dfN in data_df_2.groupby("N"):
    dfN = dfN.sort_values("lambda")
    lambda_values=dfN["lambda"]
    P_end_values=dfN["P_end"]
    f = interp1d(
        lambda_values,
        P_end_values,
        kind="cubic",      # o "linear" si hay pocos puntos
        bounds_error=False,
        fill_value="extrapolate"
    )
    
    P_end_lc_35[count] = f(lambda_p_35[count])
    P_end_lc_2_35[count] = f(lambda_p_2_35[count])
    count+=1
    
# Logarithms
logN = np.log(N)
log_Pend_35 = np.log(P_end_lc_35)
log_Pend2_35 = np.log(P_end_lc_2_35)

# Linear regressions
m1_35, dm1_35, b1_35, db1_35, r1_35 = calc_reg(logN, log_Pend_35)
m2_35, dm2_35, b2_35, db2_35, r2_35 = calc_reg(logN, log_Pend2_35)
b1_nu_35=-m1_35
db1_nu_35=dm1_35
b2_nu_35=-m2_35 
db2_nu_35=dm2_35

print("tau")
print("beta/nu:",b1_nu_35,"+-",db1_nu_35)

print("tau^2")
print("beta/nu:",b2_nu_35,"+-",db2_nu_35)

b_nu_35=(b1_nu_35+b2_nu_35)/2
db_nu_35=0.5*np.sqrt((db1_nu_35)**2+(db2_nu_35)**2)
print("beta/nu:", b_nu_35,"+-",db_nu_35)

# Smooth curves
Nfit = np.logspace(np.log10(N.min()), np.log10(N.max()), 500)

Pend_fit_35 = np.exp(b1_35) * Nfit**m1_35
Pend2_fit_35 = np.exp(b2_35) * Nfit**m2_35


#%%
#gamma=2.5
print("------------------------------------------")
print("gamma=2.5")
print("------------------------------------------")
lambda_p_25=np.zeros((len(N_list)))

lp_25=np.zeros((7,2))
tau_lp_25=np.zeros((7,2))

lp_2_25=np.zeros((7,2))
tau2_lp_25=np.zeros((7,2))

count=0
for N, dfN in data_df_1.groupby("N"):
    dfN = dfN.sort_values("lambda")
    y=dfN["tau_mean"]
    y_2=dfN["tau2_mean"]
    x=dfN["lambda"]
    lp_25[count,0],tau_lp_25[count,0],lp_25[count,1],tau_lp_25[count,1]\
        =quadratic_peak(x, y, window=2)
    lp_2_25[count,0],tau2_lp_25[count,0],lp_2_25[count,1],tau2_lp_25[count,1]\
        =quadratic_peak(x, y_2, window=2)
    count+=1
    
#%% 
#find the position of lambda_p
  
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Finite-size scaling ansatz
def lambda_peak(N, lambda_c, A, nu):
    return lambda_c + A * N**(-1/nu)

# N, lambda_p y sus errores
lambda_p_25 = lp_25[:,0]
lambda_p_err_25 = lp_25[:,1]

lambda_p_2_25 = lp_2_25[:,0]
lambda_p_2_err_25 = lp_2_25[:,1]

N=np.array(N_list)

# Initial guess
p0 = [lambda_p_25[-1], 1.0, 3.0]


# Fit
popt, pcov = curve_fit(
    lambda_peak,
    N,
    lambda_p_25,
    sigma=lambda_p_err_25,
    absolute_sigma=True,
    p0=p0
)

lambda_c_25, A_25, nu_25 = popt
dlambda_c_25, dA_25, dnu_25 = np.sqrt(np.diag(pcov))


p0 = [lambda_p_2_25[-1], 1.0, 3.0]
# Fit
popt_2, pcov_2 = curve_fit(
    lambda_peak,
    N,
    lambda_p_2_25,
    sigma=lambda_p_2_err_25,
    absolute_sigma=True,
    p0=p0
)

lambda_c_2_25, A_2_25, nu_2_25 = popt_2
dlambda_c_2_25, dA_2_25, dnu_2_25 = np.sqrt(np.diag(pcov_2))

print("tau")
print("lambda_c:",lambda_c_25,"+-",dlambda_c_25)

print("tau^2")
print("lambda_c:",lambda_c_2_25,"+-",dlambda_c_2_25)

lambda_c_25_f=(lambda_c_2_25+lambda_c_25)/2
d_lambda_c_25_f=0.5*np.sqrt((dlambda_c_2_25)**2+(dlambda_c_25)**2)
print("lambda_c:",lambda_c_25_f,"+-",d_lambda_c_25_f)


print("tau")
print("nu:",nu_25,"+-",dnu_25)

print("tau^2")
print("nu:",nu_2_25,"+-",dnu_2_25)

nu_25_f=(nu_2_25+nu_25)/2
d_nu_25_f=0.5*np.sqrt((dnu_2_25)**2+(dnu_25)**2)
print("nu:",nu_25_f,"+-",d_nu_25_f)

#%%
#find the gamma_n
# Logarithms
logN = np.log(N)
log_tau_25 = np.log(tau_lp_25[:,0])
log_tau2_25 = np.log(tau2_lp_25[:,0])

# Linear regressions: y = m x + b
m1_25, dm1_25, b1_25, db1_25, r1_25 = calc_reg(logN, log_tau_25)
m2_25, dm2_25, b2_25, db2_25, r2_25 = calc_reg(logN, log_tau2_25)
g1_nu_25=m1_25 
dg1_nu_25=dm1_25
g2_nu_25=m2_25 
dg2_nu_25=dm2_25
print("gamma_1/nu:", g1_nu_25,"+-",dg1_nu_25)
print("gamma_2/nu:", g2_nu_25,"+-",dg2_nu_25)

# Smooth curves
Nfit = np.logspace(np.log10(N.min()), np.log10(N.max()), 500)

tau_fit_25 = np.exp(b1_25) * Nfit**m1_25
tau2_fit_25 = np.exp(b2_25) * Nfit**m2_25


#%%
#P_end
from scipy.interpolate import interp1d

P_end_lc_25=np.zeros((7))
P_end_lc_2_25=np.zeros((7))
count=0
for N_i,dfN in data_df_1.groupby("N"):
    dfN = dfN.sort_values("lambda")
    lambda_values=dfN["lambda"]
    P_end_values=dfN["P_end"]
    f = interp1d(
        lambda_values,
        P_end_values,
        kind="cubic",      # o "linear" si hay pocos puntos
        bounds_error=False,
        fill_value="extrapolate"
    )
    
    P_end_lc_25[count] = f(lambda_p_25[count])
    P_end_lc_2_25[count] = f(lambda_p_2_25[count])
    count+=1
    
# Logarithms
logN = np.log(N)
log_Pend_25 = np.log(P_end_lc_25)
log_Pend2_25 = np.log(P_end_lc_2_25)

# Linear regressions
m1_25, dm1_25, b1_25, db1_25, r1_25 = calc_reg(logN, log_Pend_25)
m2_25, dm2_25, b2_25, db2_25, r2_25 = calc_reg(logN, log_Pend2_25)
b1_nu_25=-m1_25
db1_nu_25=dm1_25
b2_nu_25=-m2_25
db2_nu_25=dm2_25

print("tau")
print("beta/nu:",b1_nu_25,"+-",db1_nu_25)

print("tau^2")
print("beta/nu:",b2_nu_25,"+-",db2_nu_25)

b_nu_25=(b1_nu_25+b2_nu_25)/2
db_nu_25=0.5*np.sqrt((db1_nu_25)**2+(db2_nu_25)**2)
print("beta/nu:", b_nu_25,"+-",db_nu_25)

# Smooth curves
Nfit = np.logspace(np.log10(N.min()), np.log10(N.max()), 500)

Pend_fit_25 = np.exp(b1_25) * Nfit**m1_25
Pend2_fit_25 = np.exp(b2_25) * Nfit**m2_25




#%%
cmap = plt.get_cmap("viridis")
titles = [r"(a) $\gamma=2.5$", r"(b) $\gamma=3.5$"]

fig, ax = plt.subplots(1, 2, figsize=(14, 7))
fig.subplots_adjust(wspace=0.4)

x = np.arange(0, 0.1001, 0.001)

# -------------------- (a) gamma=2.5 --------------------

ax[0].errorbar(N**(-1/nu_25), lambda_p_25, yerr=lambda_p_err_25, fmt="o",
               capsize=3,color=cmap(0), label=r"Data: $\langle\tau\rangle$")

ax[0].plot(x, lambda_c_25 + A_25*x, "--", color=cmap(0),
           label=r"Fit: $\langle\tau\rangle$")

ax[0].errorbar(N**(-1/nu_2_25), lambda_p_2_25, yerr=lambda_p_2_err_25, fmt="s",
               capsize=3,color=cmap(0.5), label=r"Data: $\langle\tau^2\rangle$")

ax[0].plot(x, lambda_c_2_25 + A_2_25*x, "--", color=cmap(0.5),
           label=r"Fit: $\langle\tau^2\rangle$")

ax[0].set_xlim([0, 0.04])
ax[0].set_ylim([0, 0.04])
ax[0].set_xticks(np.arange(0, 0.041, 0.01))
ax[0].set_yticks(np.arange(0, 0.041, 0.01))

# -------------------- (b) gamma=3.5 --------------------

ax[1].errorbar(N**(-1/nu_35), lambda_p_35, yerr=lambda_p_err_35, fmt="o", capsize=3,
               color=cmap(0), label=r"Data: $\langle\tau\rangle$")

ax[1].plot(x, lambda_c_35 + A_35*x, "--", color=cmap(0),
           label=r"Fit: $\langle\tau\rangle$")

ax[1].errorbar(N**(-1/nu_2_35), lambda_p_2_35, yerr=lambda_p_2_err_35, fmt="s",
               capsize=3,
               color=cmap(0.5), label=r"Data: $\langle\tau^2\rangle$")

ax[1].plot(x, lambda_c_2_35 + A_2_35*x, "--", color=cmap(0.5),
           label=r"Fit: $\langle\tau^2\rangle$")

ax[1].set_xlim([0, 0.10])

# -------------------- Estilo --------------------

for i in range(2):
    ax[i].set_xlabel(r"$N^{-1/\nu}$", fontsize=18)
    ax[i].set_ylabel(r"$\lambda_p$", fontsize=18)
    ax[i].tick_params(axis="both", labelsize=16)
    ax[i].grid(True, which="both", linestyle="--", alpha=0.4)
    ax[i].set_title(titles[i], y=-0.3, fontsize=20)

handles, labels = ax[0].get_legend_handles_labels()

fig.legend(handles, labels, fontsize=18, loc="lower center",
           bbox_to_anchor=(0.5, .9), ncol=4, markerscale=1)

plt.tight_layout(rect=[0, 0, 1, 0.90])
plt.show()

name="plot_3"
plt.savefig("../plots/"+name+".pdf",bbox_inches="tight")
plt.show()

#%%
cmap = plt.get_cmap("viridis")
titles = [r"(a) $\gamma=2.5$", r"(b) $\gamma=3.5$"]

fig, ax = plt.subplots(1, 2, figsize=(14, 7))
fig.subplots_adjust(wspace=0.4)

# -------------------- (a) gamma=2.5 --------------------

ax[0].errorbar(N, tau_lp_25[:,0], yerr=tau_lp_25[:,1], fmt="o", capsize=3,
               markersize=6, color=cmap(0), label=r"Data: $\langle\tau\rangle_{\max}$")

ax[0].errorbar(N, tau2_lp_25[:,0], yerr=tau2_lp_25[:,1], fmt="s", capsize=3,
               markersize=6, color=cmap(0.5), label=r"Data: $\langle\tau^2\rangle_{\max}$")

ax[0].plot(Nfit, tau_fit_25, "--", color=cmap(0),
           label=r"Fit: $\langle\tau\rangle$")

ax[0].plot(Nfit, tau2_fit_25, "--", color=cmap(0.5),
           label=r"Fit: $\langle\tau^2\rangle$")

ax[0].set_xscale("log")
ax[0].set_yscale("log")

# -------------------- (b) gamma=3.5 --------------------

ax[1].errorbar(N, tau_lp_35[:,0], yerr=tau_lp_35[:,1], fmt="o", capsize=3,
               markersize=6, color=cmap(0), label=r"Data: $\langle\tau\rangle_{\max}$")

ax[1].errorbar(N, tau2_lp_35[:,0], yerr=tau2_lp_35[:,1], fmt="s", capsize=3,
               markersize=6, color=cmap(0.5), label=r"Data: $\langle\tau^2\rangle_{\max}$")

ax[1].plot(Nfit, tau_fit_35, "--", color=cmap(0),
           label=r"Fit: $\langle\tau\rangle$")

ax[1].plot(Nfit, tau2_fit_35, "--", color=cmap(0.5),
           label=r"Fit: $\langle\tau^2\rangle$")

ax[1].set_xscale("log")
ax[1].set_yscale("log")

# -------------------- Estilo --------------------

for i in range(2):
    ax[i].set_xlabel(r"$N$", fontsize=18)
    ax[i].set_ylabel(r"$\langle\tau^n\rangle_{\max}$", fontsize=18)
    ax[i].tick_params(axis="both", labelsize=16)
    ax[i].grid(True, which="both", linestyle="--", alpha=0.4)
    ax[i].set_title(titles[i], y=-0.3, fontsize=20)

handles, labels = ax[0].get_legend_handles_labels()

fig.legend(handles, labels, fontsize=18, loc="lower center",
           bbox_to_anchor=(0.5, .9), ncol=4, markerscale=1)

plt.tight_layout(rect=[0, 0, 1, 0.90])
plt.show()

name="plot_4"
plt.savefig("../plots/"+name+".pdf",bbox_inches="tight")
plt.show()

#%%
cmap = plt.get_cmap("viridis")
titles = [r"(a) $\gamma=2.5$", r"(b) $\gamma=3.5$"]

fig, ax = plt.subplots(1, 2, figsize=(14, 7))
fig.subplots_adjust(wspace=0.4)

# -------------------- (a) gamma=2.5 --------------------

ax[0].plot(N, P_end_lc_25, "o", markersize=6, color=cmap(0),
           label=r"Data: $\langle\tau\rangle$")

ax[0].plot(N, P_end_lc_2_25, "s", markersize=6, color=cmap(0.5),
           label=r"Data: $\langle\tau^2\rangle$")

ax[0].plot(Nfit, Pend_fit_25, "--", color=cmap(0),
           label=r"Fit: $\langle\tau\rangle$")

ax[0].plot(Nfit, Pend2_fit_25, "--", color=cmap(0.5),
           label=r"Fit: $\langle\tau^2\rangle$")

ax[0].set_xscale("log")
ax[0].set_yscale("log")

# -------------------- (b) gamma=3.5 --------------------

ax[1].plot(N, P_end_lc_35, "o", markersize=6, color=cmap(0),
           label=r"Data: $\langle\tau\rangle$")

ax[1].plot(N, P_end_lc_2_35, "s", markersize=6, color=cmap(0.5),
           label=r"Data: $\langle\tau^2\rangle$")

ax[1].plot(Nfit, Pend_fit_35, "--", color=cmap(0),
           label=r"Fit: $\langle\tau\rangle$")

ax[1].plot(Nfit, Pend2_fit_35, "--", color=cmap(0.5),
           label=r"Fit: $\langle\tau^2\rangle$")

ax[1].set_xscale("log")
ax[1].set_yscale("log")

# -------------------- Estilo --------------------

for i in range(2):
    ax[i].set_xlabel(r"$N$", fontsize=18)
    ax[i].set_ylabel(r"$P_{\rm end}(\lambda_p)$", fontsize=18)
    ax[i].tick_params(axis="both", labelsize=16)
    ax[i].grid(True, which="both", linestyle="--", alpha=0.4)
    ax[i].set_title(titles[i], y=-0.3, fontsize=20)

handles, labels = ax[0].get_legend_handles_labels()

fig.legend(handles, labels, fontsize=18, loc="lower center",
           bbox_to_anchor=(0.5, .9), ncol=4, markerscale=1)

plt.tight_layout(rect=[0, 0, 1, 0.90])
plt.show()

name="plot_5"
plt.savefig("../plots/"+name+".pdf",bbox_inches="tight")
plt.show()


#%%
#plots collapsed


fig, ax = plt.subplots(3, 2, figsize=(8, 15), sharex=False)

plt.subplots_adjust(
    left=0.10,
    right=0.98,
    top=0.88,
    bottom=0.12,
    hspace=0.30,
    wspace=0.30
)

handles = []
labels = [r"$N=10^4$", r"$N=3\times10^4$", r"$N=5\times10^4$",
          r"$N=10^5$", r"$N=3\times10^5$", r"$N=5\times10^5$",
          r"$N=10^6$"]

cmap = plt.get_cmap("viridis")

# =====================================================
# Columna izquierda: gamma = 2.5
# =====================================================

count = 0
for N, dfN in data_df_1.groupby("N"):

    dfN = dfN.sort_values("lambda")
    x = (dfN["lambda"] - lambda_c_25_f) * N**(1 / nu_25)

    line, = ax[0,0].plot(
        x,
        N**(-g1_nu_25) * dfN["tau_mean"],
        "o--",
        markersize=2,
        color=cmap(count/7)
    )

    handles.append(line)

    ax[1,0].plot(
        x,
        N**(-g2_nu_25) * dfN["tau2_mean"],
        "o--",
        markersize=2,
        color=cmap(count/7)
    )

    ax[2,0].plot(
        x,
        N**(b_nu_25) * dfN["P_end"],
        "o--",
        markersize=2,
        color=cmap(count/7)
    )

    count += 1


# =====================================================
# Columna derecha: gamma = 3.5
# =====================================================

count = 0
for N, dfN in data_df_2.groupby("N"):

    dfN = dfN.sort_values("lambda")
    x = (dfN["lambda"] - lambda_c_35_f) * N**(1 / nu_35)

    ax[0,1].plot(
        x,
        N**(-g1_nu_35) * dfN["tau_mean"],
        "o--",
        markersize=2,
        color=cmap(count/7)
    )

    ax[1,1].plot(
        x,
        N**(-g2_nu_35) * dfN["tau2_mean"],
        "o--",
        markersize=2,
        color=cmap(count/7)
    )

    ax[2,1].plot(
        x,
        N**(b_nu_35) * dfN["P_end"],
        "o--",
        markersize=2,
        color=cmap(count/7)
    )

    count += 1


# =====================================================
# Límites gamma = 2.5
# =====================================================

ax[0,0].set_xlim(0, 1.5)
ax[0,0].set_ylim(1.75, 2.75)

ax[1,0].set_xlim(0, 1.5)
ax[1,0].set_ylim(0, 0.5)

ax[2,0].set_xlim(0, 1.5)
ax[2,0].set_ylim(0, 6)


# =====================================================
# Límites gamma = 3.5
# =====================================================

ax[0,1].set_xlim(0, 1)
ax[0,1].set_ylim(4, 7.5)

ax[1,1].set_xlim(0, 1)
ax[1,1].set_ylim(0, 35)

ax[2,1].set_xlim(0, 1)
ax[2,1].set_ylim(0, 1000)


# =====================================================
# Etiquetas
# =====================================================

ylabel = [
    r"$\langle\tau\rangle N^{\gamma_1/\nu}$",
    r"$\langle\tau^2\rangle N^{\gamma_2/\nu}$",
    r"$P_{\rm end}N^{\beta/\nu}$"
]

for i in range(3):

    ax[i,0].set_ylabel(ylabel[i], fontsize=16)

    ax[i,0].tick_params(axis="both", labelsize=14)
    ax[i,1].tick_params(axis="both", labelsize=14)

    ax[i,0].grid(True, linestyle="--", alpha=0.4)
    ax[i,1].grid(True, linestyle="--", alpha=0.4)

ax[2,0].set_xlabel(r"$(\lambda-\lambda_c)N^{1/\nu}$", fontsize=16)
ax[2,1].set_xlabel(r"$(\lambda-\lambda_c)N^{1/\nu}$", fontsize=16)


# =====================================================
# Leyenda única
# =====================================================

fig.legend(
    handles,
    labels,
    loc="upper center",
    bbox_to_anchor=(0.5, 1.0),
    ncol=4,
    fontsize=15,
    markerscale=3
)


# =====================================================
# Subtítulos
# =====================================================

ax[2,0].set_title(r"(a) $\gamma=2.5$", y=-0.55, fontsize=16)
ax[2,1].set_title(r"(b) $\gamma=3.5$", y=-0.55, fontsize=16)

plt.savefig("../plots/plot_collapsed.pdf", bbox_inches="tight")
plt.show()