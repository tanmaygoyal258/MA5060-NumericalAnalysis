# MA5060: NUMERICAL ANALYSIS
# TANMAY GOYAL
# AI20BTECH11021

# NOTE: 1. We take the input as Ax = b, and output x
# NOTE: 2. We assume A is a square matrix 
# NOTE : 3. A is diagonally dominant ->  can be rearranged by user

import numpy as np
import warnings

warnings.filterwarnings('ignore')

def gauss_jacobi(A , b):
    x_prev = np.array([0 for i in range(n)] , np.float)
    x = np.array([0 for i in range(n)] , np.float)

    tolerance = 0.0001
    iters = 0
    result = -1
    for k in range(1000):

        for i in range(n):
            
            if A[i][i] == 0:
                flag = -1
                for k in range(n):
                    if A[k][i] != 0:
                        flag = k
                        break
                temp = A[flag][:]
                A[flag][:] = A[i][:]
                A[i][:] = temp
                temp = b[flag]
                b[flag] = b[i]
                b[i] = temp
                
            x[i] = -1 * b[i]    

            for j in range(n):
                if j == i:
                    continue
                else:
                    x[i] += A[i][j] * x_prev[j]
            x[i] /= -A[i][i]
            
            if np.sqrt(np.sum((x_prev - x) ** 2)) < tolerance:
                result = 1
                iters = k
                break
            else:
                x_prev = x.copy()

    if result < 0:
        print("No convergence occurs.")
    else:
        print("x = {}. It was found in {} iterations by Gauss- Jacobi method. ".format(x , iters))

def gauss_seidel(A,b):
    x_prev = np.array([0 for i in range(n)] , np.float)
    x = np.array([0 for i in range(n)] , np.float)

    tolerance = 0.0001
    iters = 0
    result = -1
    for k in range(1000):

        for i in range(n):
            
            if A[i][i] == 0:
                flag = -1
                for k in range(n):
                    if A[k][i] != 0:
                        flag = k
                        break
                temp = A[flag][:]
                A[flag][:] = A[i][:]
                A[i][:] = temp
                temp = b[flag]
                b[flag] = b[i]
                b[i] = temp
                
            x[i] = -1 * b[i]    

            for j in range(n):
                if j == i:
                    continue
                else:
                    if j >= i:
                        x[i] += A[i][j] * x_prev[j]
                    else:
                        x[i] += A[i][j] * x[j]
            x[i] /= -A[i][i]
            
            if np.sqrt(np.sum((x_prev - x) ** 2)) < tolerance:
                result = 1
                iters = k
                break
            else:
                x_prev = x.copy()

    if result < 0:
        print("No convergence occurs.")
    else:
        print("x = {}. It was found in {} iterations by Gauss-Seidel method. ".format(x , iters))



if __name__ == '__main__':

    n = -1 
    while(n<=0):
        n = int(input("Enter the dimension of the square matrix: "))
        if(n <= 0):
            print("Incorrect Dimensions. Please try again.")

    A = []
    b = []

    for i in range(n*n):
        A.append(int(input("Enter a_{}{}: ".format(i//n + 1 , i%n + 1))))
    for i in range(n):
        b.append(int(input("Enter b_{}: ".format(i))))


    A = np.array(A , dtype = np.float).reshape((n,n))
    b = np.array(b , dtype = np.float)

    print("A = \n{}".format(A))
    print("b = {}".format(b))
    gauss_jacobi(A,b)
    gauss_seidel(A,b)