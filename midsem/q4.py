from math import log, e, sqrt
import numpy as np

# --------------------* Linear Fit *--------------------
def linear_fit(X, Y, variance):
    n = len(X)  # number of datapoints

    s, sx, sy, sxx, sxy = 0, 0, 0, 0, 0
    for i in range(n):
        s += 1 / variance[i] ** 2
        sx += X[i] / variance[i] ** 2
        sy += Y[i] / variance[i] ** 2
        sxx += X[i] ** 2 / variance[i] ** 2
        sxy += X[i] * Y[i] / variance[i] ** 2

    delta = s * sxx - sx**2
    a = (sxx * sy - sx * sxy) / delta
    b = (s * sxy - sx * sy) / delta

    # calculate chi^2 / dof
    dof = n - 2
    chi2 = 0
    for i in range(n):
        chi2 += (Y[i] - a - b * X[i]) ** 2 / variance[i] ** 2

    delA2 = sxx / delta
    delB2 = s / delta
    cov = -sx / delta
    return a, b, delA2, delB2, cov, chi2


# --------------------* Reading the data  *--------------------
time = []
counts = []
uncertainty = []

with open("midsem/data/msfit.txt", "r") as file:
    for line in file:
        a, b, c = line.split()
        time.append(int(a))
        counts.append(int(b))
        uncertainty.append(int(c))


# --------------------* Average lifetime is given by A = A0 exp(-t/tau) *--------------------
X = time.copy()
Y = counts.copy()
sigma = uncertainty.copy()
variance = [0 for i in range(len(X))]

for i in range(len(X)):
    Y[i] = log(Y[i], e)
    variance[i] = log(sigma[i] ** 2, e)

a, b, delA2, delB2, cov, chi2 = linear_fit(X, Y, variance)

dof = len(X) - 2
lifetime = -1 / b
lifetime_error = lifetime**2 * delB2

t_crit = 1.860  # at 95% significance
N = len(X)
n = 5
ybar = 0
avg = 0
for i in range(n):
    ybar += counts[i] / n

for j in range(N):
    avg += counts[i] / N

sigma_full = sqrt(N)
t = (ybar - avg) / (sigma_full / sqrt(n))


print("lifetime:\t", lifetime)
print("lifetime_error:\t", lifetime_error)
print("t:\t\t", t)
print("t_crit:\t\t", t_crit)


# --------------------* OUTPUT *--------------------
# lifetime:        91.08264534723905
# lifetime_error:  6.993032401147983
# t:               8.909545442950504
# t_crit:          1.86
