# MA5060: NUMERICAL ANALYSIS
# TANMAY GOYAL
# AI20BTECH11021
# NOTE: 1. We assume A is a square matrix 

# example matrix can be [1,2,4,2,1,2,4,2,1]

import numpy as np

n = -1 
while(n<=0):
    n = int(input("Enter the dimension of the square matrix: "))
    if(n <= 0):
        print("Incorrect Dimensions. Please try again.")

A = []

for i in range(n*n):
    A.append(int(input("Enter a_{}{}: ".format(i//n + 1 , i%n + 1))))

A = np.array(A , dtype = np.float64).reshape((n,n))

print("A = \n{}".format(A))

x_prev = np.array([1 for i in range(n)] , dtype = np.float64)

tolerance = 0.0001
result = -1

for k in range(1000):

    y = A @ x_prev
    eig = y.max()
    x = y / eig

    if np.sqrt(np.sum((x - x_prev)**2)) <= tolerance:
        print("Largest eigenvalue in modulus is: {}".format(eig))
        print("By comparision, the largest eigenvalue of A in modulus is {}".format(np.max(np.abs(np.linalg.eigvals(A)))))
        print("Dominant eigenvector is: {}".format(x))
        result = 1
        break

    else:
        x_prev = x

if result < 0:
    print("Convergence doesnt happen.")        
