# MA5060: NUMERICAL ANALYSIS
# TANMAY GOYAL
# AI20BTECH11021

# NOTE: 1. The Lower Triangular matrix is assumed to have 1 on the primary diagonal.
# NOTE: 2. We calculate the LU decomposition only for square matrices.

# Some simple calculations show:
# l_{i1} = a_{i1} / a_{11}
# l_{ij} = [a_{i_j} - sum from k=1 to j-1 l_{ik}u_{kj}] / u_{jj}
# u_{1i} = a_{1i}
# u_{ij} = a_{ij} -  sum from k=1 to i-1 l_{ik}u_{kj}

import numpy as np

n = -1 
while(n<=0):
    n = int(input("Enter the dimension of the square matrix: "))
    if(n <= 0):
        print("Incorrect Dimensions. Please try again.")

A = []

for i in range(n*n):
    A.append(int(input("Enter a_{}{}: ".format(i//n + 1 , i%n + 1))))

A = np.array(A).reshape((n , n))

L = np.zeros_like(A)
U = np.zeros_like(A)

# initialization
for i in range(n):
    for j in range(n):
        if i == j:
            L[i][j] = 1
        elif j == 0:
            L[i][j] = A[i][j] / A[0][0]

U[0][:] = A[0][:]

for i in range(1 , n):

    # update ith row in U
    for j in range(i , n):
        U[i][j] = A[i][j]
        for k in range(i):
            U[i][j] -= L[i][k] * U[k][j]
    
    # update ith column in L
    for j in range(i+1,n):
        L[j][i] = A[j][i]
        for k in range(i):
            L[j][i] -= L[j][k] * U[k][i]
        L[j][i] /= U[i][i]

print("A = \n{}".format(A))
print("Lower Triangular Matrix = \n{}".format(L))
print("Upper Triangular Matrix = \n{}".format(U))