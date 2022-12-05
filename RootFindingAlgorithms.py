# MA5060: NUMERICAL ANALYSIS
# TANMAY GOYAL
# AI20BTECH11021


# NOTE: THIS CODE IS ONY INTENDED TO WORK FOR POLYNOMIAL FUNCTIONS

import numpy as np

MAX_ITERS = 100000

def val(coeff , a):
    sum = 0
    for i in range(len(coeff)):
        sum += coeff[i] * pow(a , i)

    return sum


def derv_coeff(coeff):
    derv = []

    for i in range(len(coeff)-1 , 1 , -1):
        derv.append(coeff[i] * i)
    
    derv = derv[::-1]

    return derv


def bisection(coeff , a , b , err_tol = 0.00001):
    
    n_iters = 0
    flag = False
    
    while not flag and  b - a >= err_tol and n_iters <= MAX_ITERS:

        n_iters += 1
        mid = a + (b-a)/2
        
        if val(coeff,mid) == 0:
            print("The EXACT root is : {} and was found in {} iterations".format(mid , n_iters))
            flag = True
        
        elif val(coeff,mid) * val(coeff,a) < 0:
            b = mid
        
        else:
             a = mid

    if not flag:
        print("The APPROXIMATE root is : {} and was found in {} iterations".format(a+(b-a)/2 , n_iters))


def secant(coeff , x_k_ , x_k , err_tol = 0.00001):
    
    n_iters = 0
    flag = False
    
    while not flag and  abs(x_k_ - x_k) >= err_tol and n_iters <= MAX_ITERS:

        n_iters += 1
        x_m = x_k - ((x_k - x_k_) * val(coeff , x_k) / (val(coeff , x_k) - val(coeff , x_k_)))
        
        if val(coeff , x_m) == 0:
            print("The EXACT root is {} and was found in {} iterations".format(x_k , n_iters))
            flag = True

        x_k_ = x_k
        x_k = x_m

    if not flag:
        print("The APPROXIMATE root is : {} and was found in {} iterations".format(x_k , n_iters))


def regula_false(coeff , x_k_ , x_k , err_tol = 0.00001):

    n_iters = 0
    flag = False
    
    while not flag and  abs(val(coeff , x_k)) >= err_tol and n_iters <= MAX_ITERS:

        n_iters += 1
        x_m = x_k - ((x_k - x_k_) * val(coeff , x_k) / (val(coeff , x_k) - val(coeff , x_k_)))

        if val(coeff , x_m) == 0:
            print("The EXACT root is {} and was found in {} iterations".format(x_k , n_iters))
            flag = True

        elif val(coeff , x_k_) * val(coeff , x_m):
            x_k = x_m
        
        elif val(coeff , x_k) * val(coeff , x_m):
            x_k_ = x_m

        else:
            print("Root not found")
            flag = True
        

    if not flag:
        print("The APPROXIMATE root is : {} and was found in {} iterations".format(x_k , n_iters))


def newton_raphson(coeff , x_k , err_tol = 0.00001):
    
    n_iters = 0
    flag = False
    derv = derv_coeff(coeff)

    while not flag and abs(val(coeff , x_k)) >= err_tol and n_iters <= MAX_ITERS:
        
        n_iters += 1
        x_m = x_k - (val(coeff , x_k) / val(derv , x_k))

        if val(coeff , x_m) == 0:
            print("The EXACT root is {} and was found in {} iterations".format(x_k , n_iters))
            flag = True

        x_k = x_m

    if not flag:
        print("The APPROXIMATE root is : {} and was found in {} iterations".format(x_k , n_iters))


def phi(coeff , a):
    degree = len(coeff) - 1
    value = -val(coeff[:-1] , a)/coeff[-1]
    if degree % 2 == 0:
        if value < 0:
            return (np.copysign(np.abs(value) ** (1. / degree), value))
        else:
            return (np.copysign(np.abs(value) ** (1. / degree), value))
    else:
        return (np.copysign(np.abs(value) ** (1. / degree), value))


def iterative(coeff , x_k , err_tol = 0.00001):
    n_iters = 0
    flag = False

    while not flag and n_iters <= MAX_ITERS:
        
        n_iters += 1
        
        v = phi(coeff , x_k)
            
        if v == x_k:
            print("The EXACT root is {} and was found in {} iterations".format(x_k , n_iters))
            flag = True

        elif abs(v - x_k) < err_tol:
            print("The APPROXIMATE root is : {} and was found in {} iterations".format(x_k , n_iters))
            flag = True

        x_k = v

    if not flag:
        print("The APPROXIMATE root is : {} and was found in {} iterations".format(x_k , n_iters))
        



if __name__=='__main__':

    print("Choose the method you wish to apply: ")
    print("1. Bisection Method")
    print("2. Regula- False Method")
    print("3. Secant Method")
    print("4. Newton-Raphson Method")
    print("5. General Iteration Scheme")
    choice = int(input("Enter choice: "))

    print("Enter the degree of the function you wish to approximate: ")
    n = int(input("Enter degree: "))
    assert(n >= 0)

    coeff = []
    for i in range(n+1):
        coeff.append(float(input("Enter coefficient of x^{}: ".format(i))))

    if choice==1:
        a = float(input("Enter lower bound for root: "))
        b = float(input("Enter upper bound for root: "))
        if val(coeff , a) == 0:
            print("The EXACT root for the function is {}".format(a))
        elif val(coeff , b) == 0:
            print("The EXACT root for the function is {}".format(b))
        elif val(coeff , a) * val(coeff , b) > 0:
            print("Incorrect bounds given.")
        else:
            bisection(coeff , a , b)

    elif choice == 2:
        x_k_ = float(input("Enter lower bound for root: "))
        x_k = float(input("Enter upper bound for root: "))
        if val(coeff , x_k_) == 0:
            print("The EXACT root for the function is {}".format(x_k_))
        elif val(coeff , x_k) == 0:
            print("The EXACT root for the function is {}".format(x_k))
        elif val(coeff , x_k_) * val(coeff , x_k) > 0:
            print("Incorrect bounds given.")
        else:
            regula_false(coeff , x_k_ , x_k)

    elif choice == 3:
        x_k_ = float(input("Enter lower bound for root: "))
        x_k = float(input("Enter upper bound for root: "))
        if val(coeff , x_k_) == 0:
            print("The EXACT root for the function is {}".format(x_k_))
        elif val(coeff , x_k) == 0:
            print("The EXACT root for the function is {}".format(x_k))
        elif val(coeff , x_k_) * val(coeff , x_k) > 0:
            print("Incorrect bounds given.")
        else:
            secant(coeff , x_k_ , x_k)

    elif choice == 4:
        x_0 = float(input("Enter approximation for root: "))
        if val(coeff , x_0) == 0:
            print("The EXACT root for the function is {}".format(x_0))
        else:
            newton_raphson(coeff , x_0)

    else:
        x_0 = float(input("Enter approximation for root: "))
        if abs(val(derv_coeff(coeff) , x_0)) >= 1:
            print("The iteration scheme will not converge")
        elif val(coeff , x_0) == 0:
            print("The EXACT root for the function is {}".format(x_0))
        else:
            iterative(coeff , x_0)
