from math import sqrt, pi
import numpy as np
from numpy import cos, sin

# --------------------* Fucntion to read from file *--------------------
def read_matrix(file):
    with open(file, "r") as f:
        a = [[int(num) for num in line.split()] for line in f]
    return a


# --------------------* Function to normaliz vector *--------------------
def normalize(vec):
    norm = np.linalg.norm(vec, 2)
    return [i / (sqrt(norm)) for i in vec]


# --------------------* Function for power method *--------------------
def power_method(A, x, tol, eignum=1) -> tuple:
    """r
    Function to evaluate the eigenvalues and corresponding eigenvectors for a
    given matrix `A`, a random vector `x` of same dimension, a tolerance `tol`
    and number of eigenvalues `eignum` (either 1 or 2).
    """
    n = len(A)
    x = x / np.linalg.norm(x)
    y = x.copy()
    if eignum == 1:
        diff = 1
        while diff > tol:
            xnew = A @ x
            eigval = np.dot(xnew, x) / np.dot(x, x)
            xnew = xnew / np.linalg.norm(xnew)
            diff = np.linalg.norm(xnew - x)
            x = xnew.copy()

        vec = xnew

        return eigval, vec

    elif eignum == 2:
        diff = 1
        while diff > tol:
            xnew = A @ x
            eigval1 = np.dot(xnew, x) / np.dot(x, x)
            xnew = xnew / np.linalg.norm(xnew)
            diff = np.linalg.norm(xnew - x)
            x = xnew.copy()

        vec1 = xnew

        A = A - eigval1 * np.outer(vec1, vec1.T)
        diff = 1
        while diff > tol:
            ynew = A @ y
            eigval2 = np.dot(ynew, y) / np.dot(y, y)
            ynew = ynew / np.linalg.norm(ynew)
            diff = np.linalg.norm(ynew - y)
            y = ynew.copy()

        vec2 = ynew

        return eigval1, eigval2, vec1, vec2


# --------------------* Reading the file *--------------------
A = read_matrix("midsem/data/mstrimat.txt")

eigval1, eigval2, eigvec1, eigvec2 = power_method(A, np.random.rand(len(A)), 1e-4, 2)

# Given values of eigenvalues and eigenvectors
b = 2
a, c = -1, -1
n = 5

k = 1
given_eigval1 = b + 2 * sqrt(a * c) * cos(k * pi / (n + 1))
given_eigvec1 = []
for i in range(5):
    given_eigvec1.append(2 * (sqrt(c / a)) ** k * sin(k * pi * i / (n + 1)))

k = 2
given_eigval2 = b + 2 * sqrt(a * c) * cos(k * pi / (n + 1))
given_eigvec2 = []
for i in range(5):
    given_eigvec2.append(2 * (sqrt(c / a)) ** k * sin(k * pi * 1j / (n + 1)))

print("----------* First Eigenvalue*--------------------")
print("Obtained:\t", eigval1)
print("Given:   \t", given_eigval1)
print("----------* First EigneVector*--------------------")
print("Obtained:\t", normalize(eigvec1))
print("Given:   \t", normalize(given_eigvec1))
print("----------* Second Eigenvalue*--------------------")
print("Obtained:\t", eigval2)
print("Given:   \t", given_eigval2)
print("----------* Second EigneVector*--------------------")
print("Obtained:\t", normalize(eigvec2))
print("Given:   \t", normalize(given_eigvec2))


# --------------------* OUTPUT *--------------------

# ----------* First Eigenvalue*--------------------
# Obtained:        3.732050643185673
# Given:           3.7320508075688776
# ----------* First EigneVector*--------------------
# Obtained:        [-0.28886556887642506, 0.5001904228733041, -0.5773502312264027, 0.49980950457731926, -0.2884846505804403]
# Given:           [0.0, 0.5491004867761125, 0.9510699415570292, 1.098200973552225, 0.9510699415570294]
# ----------* Second Eigenvalue*--------------------
# Obtained:        3.000000090774257
# Given:           3.0
# ----------* Second EigneVector*--------------------
# Obtained:        [-0.49978486284529383, 0.499763004870703, 0.00019530599045070705, -0.5002368736616228, 0.5002150156863868]
# Given:           [1.0571035245083658j, 1.0571035245083658j, 1.0571035245083658j, 1.0571035245083658j, 1.0571035245083658j]
