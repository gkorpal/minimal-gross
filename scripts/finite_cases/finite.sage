##########################################################################
# Get the current working directory                                      #
# Get the absolute path of the parent directory of the current directory #
# Add the parent directory to sys.path                                   #
##########################################################################
import sys
import os
current_dir = os.getcwd()
parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir, os.pardir))
sys.path.insert(0, parent_dir)

from functions.fp2.eichler.maximal import get_maximal_order
from functions.fp2.eichler.bounds import d3_bound, minimal_basis_matrix


# Check if the correct number of arguments is provided
if len(sys.argv) != 3:
    print("Usage: sage finite.sage <M> <N>")
    sys.exit(1)

M = abs(Integer(sys.argv[1])) # p >= M
N = abs(Integer(sys.argv[2])) # p <= N

with open(os.path.join(current_dir, f"cases_{M}_{N}.txt"), "w") as file:
    file.write(f"The Eisenstein-reduced Gram matrices of all maximal orders for p lying between {M} and {N}.\n\n")
    # loop over the primes between [M,N]
    counter = 0
    for p in prime_range(M, N+1):
        file.write(f"{p = }.\n")
        counter+=1
        print(f"working on prime number {counter}...")
        B.<i,j,k> = QuaternionAlgebra(p)
        a = i^2
        b = j^2
        file.write(f"Quaternion Algebra: i^2 = {a} and j^2 = {b}.\n")
        
        print(f"finding the maximal order using Magma...")
        # magma output [(a,b), base0, (base, matrix), ...]
        #upper = ceil(2*p^(2/3)) #Elkies
        #for D1 in range(1, upper+1):

        m = get_maximal_order(a, b, 0)           

        n = len(m)
        if n <= 2:
            print(f"No maximal order found.")
            file.write(f"No maximal order found.\n\n")
        else: 
            for idx in range(2,n):
                # need to find the corresponding j-invariant
                base = sage_eval(str(m[idx][0]), locals={'i': i, 'j': j, 'k': k})
                O = B.quaternion_order(basis=base)
                
                Gred = Matrix(ZZ, [[Integer(x) for x in row.strip('[]').split()] for row in str(m[idx][1]).split('\n')])

                Qform = TernaryQF([Gred[0,0], Gred[1,1], Gred[2,2], 2*Gred[1,2], 2*Gred[0,2], 2*Gred[0,1]])
                Ered = (ZZ(1) / ZZ(2)) * Qform.reduced_form_eisenstein(matrix=False).matrix()

                Mred = minimal_basis_matrix(Ered)
                # test if in Fp
                bound = d3_bound(Mred,p)
                if Mred[2][2] < p:
                    file.write(f"Maximal order: {base}.\n")
                    #file.write(f"LLL-reduced Gram matrix:\n{Gred}\n")
                    file.write(f"Reduced Gram matrix:\n{Ered}\n")
                    file.write(f"Rearranged Gram matrix with D1 <= D2 <= D3 and angles flipped if needed:\n{Mred}\n")
                    file.write(f"The j-invariant does not belong to Fp since D3<{p}. Moreover, D3 is less than {bound[1]}.\n\n")
                else:
                    file.write(f"Maximal order found, but over Fp.\n\n")

    
    if counter == 0:
        print(f"No prime number found.")
        file.write(f"No prime number found.\n")
