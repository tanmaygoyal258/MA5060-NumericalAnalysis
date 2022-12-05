# MA5060: NUMERICAL ANALYSIS
# TANMAY GOYAL
# AI20BTECH11021

# NOTE: 1. We assume A is a square matrix

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
result = 1
if np.linalg.det(A) == 0:
    print("No inverse possible")
    result = -1


print("A = \n{}".format(A))

I = np.identity(n)

# we will do row operations n times
for t in range(n):
    
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
        temp = I[flag][:]
        I[flag][:] = I[t][:]
        I[t][:] = temp
        
    divide = A[t][t]
    A[t][:] /= divide
    I[t][:] /= divide
    for i in range(n):
        if i == t:
            continue
        else:
            update = A[i][t] / A[t][t]
            A[i][:] -= update * A[t][:]
            I[i][:] -= update * I[t][:]


print("Inverse of A is \n{}".format(I))


        