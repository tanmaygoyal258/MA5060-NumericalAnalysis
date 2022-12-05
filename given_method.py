# MA5060: NUMERICAL ANALYSIS
# TANMAY GOYAL
# AI20BTECH11021
# NOTE: 1. We assume A is a square matrix 
# NOTE: 2. We assume eigenvalues lie between -10 and 10 

# example matrix can be [1,2,4,2,1,2,4,2,1]

import numpy as np

# finding value of polynommial
def val(coeff , a):
    sum = 0
    for i in range(len(coeff)):
        sum += coeff[i] * pow(a , i)

    return sum

# secant method to find roots
def secant(coeff , x_k_ , x_k , err_tol = 0.00001):
    
    n_iters = 0
    flag = False
    
    while not flag and  abs(x_k_ - x_k) >= err_tol and n_iters <= 10000:

        n_iters += 1
        x_m = x_k - ((x_k - x_k_) * val(coeff , x_k) / (val(coeff , x_k) - val(coeff , x_k_)))
        
        if val(coeff , x_m) == 0:
            print("The EXACT eigenvalue is {} ".format(x_k ))
            flag = True

        x_k_ = x_k
        x_k = x_m

    if not flag:
        print("The APPROXIMATE eigenvalue is : {}".format(x_k))

# taking inputs
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


# producing tridiagonal matrix
for i in range(1 , n):
    for j in range(i+1 , n):
        theta = np.arctan(A[i-1][j] / A[i-1][i])
        S = np.identity(n)
        S[i][i] = np.cos(theta)
        S[i][j] = -np.sin(theta)
        S[j][i] = np.sin(theta)
        S[j][j] = np.cos(theta)
        A = S.T @ A @ S

# getting sturm sequence
sturm = []
f0 = [1.0] 
f1 = [1.0 , -A[0][0]]
sturm.append(f0)
sturm.append(f1)

for i in range(1 , n):
    
    # (lamda-b_r)*f_{r-1}
    factor = np.array([1 , -A[i][i]], dtype = np.float64)
    temp = [s * factor[0] for s in f1]
    temp.append(0)
    term1 = np.array(temp, dtype = np.float64)
    term2 = [0]
    for s in f1:
        term2.append(s * factor[1])
    
    term2 = np.array(term2, dtype = np.float64)
    f = term1 + term2
    c = A[i-1][i]
    
    term3 = [0,0]
    for s in f0:
        term3.append(c**2 * s)
    term3 = np.array(term3, dtype = np.float64)
    
    # -c_{r-1}^2 f_{r-2}
    f = f - term3

    # updates
    f0 = f1
    f1 = f.tolist()
    sturm.append(f1)

characteristic_polynomial  = sturm[-1].copy()

# sign changes in Sturm Sequence
sturm_sign_table = []

for x in range(-10 , 11):
    signs = []
    for k in sturm:
        sum = 0 
        for j in range(len(k)):
            sum += k[j] * np.power(x , len(k) - j - 1)

        if sum < 0:
            signs.append('-')
        elif sum > 0:
            signs.append('+')
        else:
            signs.append('0')
            
    sturm_sign_table.append(signs)

# number of sign changes
sturm_sign_changes = []
for i in sturm_sign_table:
    sign_changes = 0
    
    for j in range(1 , len(i)):
        if i[j] != i[j-1] and i[j]!= '0':
            sign_changes += 1

    sturm_sign_changes.append(sign_changes)

# location of roots
location = []
for i in range(-9,11):
    if sturm_sign_changes[i+10] != sturm_sign_changes[i+9]:
            location.append([i-1 , i])

# finding the roots
for pair in location:
    exact = False
    for t in range(1):
        if val(characteristic_polynomial[::-1] , pair[t]) == 0:
            print("Exact eigenvalue is {}".format(pair[t]))
            exact = True

    if not exact:
        secant(characteristic_polynomial[::-1] , pair[0] , pair[1])

print("By comparision, the eigenvalues of A are: {}".format(np.linalg.eigvals(A)))