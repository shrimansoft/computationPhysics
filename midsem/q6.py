import numpy as np

# Read matrix from a file given as a string (space separated file)


def read_matrix(file):
    with open(file, "r") as f:
        a = [[float(num) for num in line.split()] for line in f]
    return a


# ------------------------* Gauss-Seidel *-------------------------


def gauss_seidel(A, b, tol):
    n = len(A)
    x = np.zeros(n)
    x0 = np.ones(n)
    iterations = []
    residue = []
    count = 0  # counts the number of iterations
    while np.linalg.norm(x - x0) > tol:
        iterations.append(count)
        count += 1
        for i in range(n):
            s1, s2 = 0, 0
            for j in range(i):
                s1 += A[i][j] * x[j]
            for j in range(i + 1, n):
                s2 += A[i][j] * x0[j]
            x[i] = 1 / A[i][i] * (b[i] - s1 - s2)
        residue.append(np.linalg.norm(x - x0))
        x0 = x.copy()
    return x, iterations, residue


# --------------------* Jacobi Method *--------------------


def jacobi(A: np.ndarray, b: np.ndarray, tol: float) -> np.ndarray:
    n = len(A)
    x = np.ones(n)  # define a dummy vector for storing solution vector
    xold = np.zeros(n)
    iterations = []
    residue = []
    count = 0
    while np.linalg.norm(xold - x) > tol:
        iterations.append(count)
        count += 1
        residue.append(np.linalg.norm(xold - x))
        xold = x.copy()
        for i in range(n):
            total = 0
            for j in range(n):
                if i != j:
                    total += A[i][j] * x[j]

            x[i] = 1 / A[i][i] * (b[i] - total)

    return x, iterations, residue


# --------------------* Reading the date and seting the inital values *--------------------
A = read_matrix("midsem/data/A.txt")
b = read_matrix("midsem/data/B.txt")
b = b[0]
tol = 1e-5

# --------------------* calcultaion by calling the funcitons *--------------------
x_jacobi = jacobi(A, b, tol)[0]
x_gs = gauss_seidel(A, b, tol)[0]

print("--------------------* Solutions *---------------------")
print("Jacobi: \t", x_jacobi)
print("Gauss-Seidel: \t", x_gs)


# --------------------* OUTPUT *-----------------------------

# --------------------* Solutions *---------------------
# Jacobi:          [ 1.49999882 -0.5         2.         -2.49999941  1.         -1.        ]
# Gauss-Seidel:    [ 0.25       -0.375       1.95833333 -0.875       1.05       -0.73333333]
