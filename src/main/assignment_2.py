import numpy as np

#Function that creates Neville's method
def NevilleDomo(x,val) :

    #W is the desired f(w), or the position you want to approximate
    w = 3.7
    #Create the Array the size of your values
    length = len(x)
    neville = np.zeros((length,length))

    #Make first slot of each row the f(x) value
    for i in range(length) :
        neville[i][0] = val[i]

    #Caculation for each approximation
    for i in range(length):
        for j in range(1,i+1) :
             termOne = (w-x[i-j]) * neville[i][j-1]
             termTwo = (w-x[i]) * neville[i-1][j-1]

             neville[i][j] = (termOne - termTwo) / (x[i] - x[i-j])

    #Return the value at the 2nd degree caluclated by the function
    return neville[2][length-1]

#Function that prints the Derived Data
def print_table(dd_table):
    #This formatting puts it in scientific notation, and matches the spacing of expected output
    for row in dd_table:
        formatted_row = ["{:.8e}".format(val) for val in row]
        print(f"[ {' '.join(formatted_row)}]")

#Function that prints the cubline matrix as desired
def print_cublic(matrix):
    for row in matrix:
        formatted_row = [("{:.0f}"+".").format(val) for val in row]
        print(f"[ {' '.join(formatted_row)}]")
        
#Function that solves the NewtonForward Approximation given x and f(x)
def NewtonForward() :
    xi = [7.2, 7.4, 7.5,7.6] 
    fxi = [23.5492,25.3913,26.8224,27.4589]
 
    lim = len(xi)
    diffs = np.zeros((4,4))

    for i in range(lim) :
        diffs[i][0] = fxi[i]

    for i in range(1, lim):
        for j in range(1, i + 1):
            diffs[i][j] = (diffs[i][j-1] - diffs[i-1][j-1]) / (xi[i] - xi[0])

    #Print the matrix created
    print(diffs[1][1])
    print(diffs[2][2])
    print(diffs[3][3])
    print()
  
    #Solves the approx with the given f(x)
    #Solves it to the 3rd degree
    p1 = fxi[0] +(7.3-xi[0])*diffs[1][1]
    p2 = p1 - diffs[2][2]*(7.3-xi[0])*(7.3-xi[1])
    p3 = p2 + diffs[3][3]*(7.3-xi[0])*(7.3-xi[1])*(7.3-xi[2])
    print(p3,"\n")
                    
    return p3

#Functio that solves the Hermite Polynomial approximation
def hermite_Approx():
    
    #Given x, f(x) and fprime, note that these could be asked through a user input
    x = [3.6,3.6,3.8,3.8,3.9,3.9] #To make coding it easier just duplicate all x and fx values in the array
    fx = [1.675,1.675,1.436,1.436,1.318,1.318]
    fprime = [-1.195,-1.188,-1.182]
    
    # Create the divided difference table, note it is x rows with x-1 columns each
    n = len(x)
    table = np.zeros((6,5))

    #Fill first column with x values
    for i in range(6):
        table[i][0] = x[i]

    # Fill second column with f(x) values
    for i in range(6):
        table[i][1] = fx[i]
    #Fill third column with either f'(x), or the calulated value
    for i in range(6):
        if (i==0) :
            table[i][2] = 0
        elif (i % 2 == 1) :
            table[i][2] = fprime[i // 2]
        else :
            table[i][2] = (table[i][1] - table[i-1][1]) / (table[i][0] -table[i-1][0])
        
    #Fill forth and fifth column with the calulated difference, note a nested loop could do this for however many columns
    for i in range(2,6):
           table[i][3] = (table[i][3-1] - table[i-1][3-1]) /(table[i][0] -table[i-2][0])
    for i in range(3,6):
            table[i][4] = (table[i][4-1] - table[i-1][4-1]) /(table[i][0] -table[i-2][0])

    print_table(table)
    print("")
    return table

#Function to solve cublic spline interpolation
def cublicSpline():
    x = [2,5,8,10]
    fx = [3,5,7,9]

    #Solve for h values for matrix A, along with create the matrix
    hZero = x[1] -x[0]
    hOne = x[2] -x[1]
    hTwo = x[3] -x[2]
    matrixA = np.zeros((4,4))

    #Solve matrix A
    matrixA[0, 0] = 1
    for i in range (2,4) :
        matrixA[1, i] = 0
    matrixA[1,0] = hZero
    matrixA[1,1] = 2*(hZero+hOne)
    matrixA[1,2] = hOne
    matrixA[2,1] = hOne
    matrixA[2,2] = 2*(hOne+hTwo)
    matrixA[2,3] = hTwo
    matrixA[3,3] = 1
    #Print the matrix
    print_cublic(matrixA)

    #Vector B time, create it and find a values
    vectorB = np.zeros((4))
    aZero = fx[0]
    aOne = fx[1]
    aTwo = fx[2]
    aThree = fx[3]

    #Solve vector B, note the first and last will always be 0
    vectorB[0] = 0. 
    vectorB[1] = (3/hOne) * (aTwo - aOne) - (3/hZero) *(aOne - aZero)
    vectorB[2] = (3/hTwo) *(aThree - aTwo) - (3/hOne) * (aTwo- aOne)
    vectorB[3] = 0.
    print(vectorB) #Print the vector

    #Solve vector x using numpy(bless it exist)
    vectorX = np.linalg.solve(matrixA,vectorB)

    #Print VectorX without having huge decimal gaps
    print("["+f"{vectorX[0]:.0f}"+".",f"{vectorX[1]:.8f}",f"{vectorX[2]:.8f}",f"{vectorX[3]:.0f}"+"."+"]") #This was overly difficult
    

#Main function that calls all calculations to be run
def main() :

    #Run Question 1, the Neville Method with hardcoded values
    x_points = [3.6, 3.8, 3.9] 
    val = [1.675, 1.436, 1.318]
    Neville_Solved = NevilleDomo(x_points, val)
    print(Neville_Solved,"\n")
    

    #Run question 2 and 3 using the NewtonForward method
    NewtonForward_solved = NewtonForward()

    #Run question 4, the Hermite Polynomial Approximation
    hermite_table = hermite_Approx()

    #Run question 5, cublic spline interpolation
    cublicSpline()
    
#Make main run 
if __name__=="__main__":
    main()
