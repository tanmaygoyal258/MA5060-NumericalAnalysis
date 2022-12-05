# MA5060: NUMERICAL ANALYSIS
# TANMAY GOYAL
# AI20BTECH11021

# NOTE: 1. We take the input as Ax = b, and output x
# NOTE: 2. We assume A is a square matrix

import numpy as np

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

A = np.array(A , dtype = np.float64).reshape((n,n))
b = np.array(b,  dtype = np.float64)

print("A = \n{}".format(A))
print("b = {}".format(b))

result = 1

# we will do row operations n times
for t in range(n):
    for i in range(t+1 , n):
        
        if A[t][t] == 0:
            # we switch with a row which does not have 0 in diagonal element
            flag = -1
            for k in range(t+1 , n):
                if A[k][t] != 0:
                    flag = k
                    break
            temp = A[flag][:]
            A[flag][:] = A[t][:]
            A[t][:] = temp
            temp = b[flag]
            b[flag] = b[t]
            b[t] = temp
            
        update = A[i][t] / A[t][t]
        A[i][:] -= update * A[t][:]
        b[i] -= update * b[t]

        if np.linalg.det(A) == 0:
            print("No solution exists")
            result = -1
            break
    
    if result < 0:
        break

if result > 0:

    # we obtain upper triangular matrix, we now do back substitution:
    x = np.zeros_like(b)
    for i in range(n-1 , -1 , -1):
        x[i] = b[i]
        for j in range(i+1 , n):
            x[i] -= A[i][j] * x[j]
        x[i] /= A[i][i]
        
    print("x = {}".format(x))


