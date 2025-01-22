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
from functions.fp2.deuring.correspondence import constructive_deuring
from functions.fp.gross.type import get_matrix_type


# Check if the correct number of arguments is provided
if len(sys.argv) != 4:
    print("Usage: sage fam20.sage <c> <M> <N>")
    sys.exit(1)

c = abs(Integer(sys.argv[1])) # p = c mod 20
M = abs(Integer(sys.argv[2])) # p >= M
N = abs(Integer(sys.argv[3])) # p <= N

with open(os.path.join(current_dir, f"fam20_{c}_{M}_{N}.txt"), "w") as file:
    file.write(f"The LLL-reduced Gram matrices of all maximal orders with the first successive minima equal to 20 and lying in Bp for p = {c} (mod 20) lying between {M} and {N}.\n\n")
    # loop over the primes c mod 20 and between [M,N]
    counter = 0
    for p in prime_range(M, N+1):
        if p % 20 == c:
            file.write(f"{p = }.\n")
            counter+=1
            print(f"working on prime number {counter}...")
            B.<i,j,k> = QuaternionAlgebra(p)
            a = i^2
            b = j^2
            file.write(f"Quaternion Algebra: i^2 = {a} and j^2 = {b}.\n")
            
            print(f"finding the maximal order using Magma...")
            # magma output [(a,b), base0, (base, matrix), ...]
            m = get_maximal_order(a, b, 20)              

            n = len(m)
            if n <= 2:
                print(f"No maximal order with D1=20 found.")
                file.write(f"No maximal order with D1=20 found.\n\n")
            else: 
                for idx in range(2,n):
                    # need to find the corresponding j-invariant
                    base = sage_eval(str(m[idx][0]), locals={'i': i, 'j': j, 'k': k})
                    O = B.quaternion_order(basis=base)
                    file.write(f"Maximal order with D1=20: {base}.\n")
                    
                    Gred = Matrix(ZZ, [[Integer(x) for x in row.strip('[]').split()] for row in str(m[idx][1]).split('\n')])
                    file.write(f"LLL-reduced Gram matrix:\n{Gred}\n")

                    Mred = minimal_basis_matrix(Gred)
                    file.write(f"Rearranged Gram matrix with D1 <= D2 <= D3 and angles flipped if needed:\n{Mred}\n")
                    # test if in Fp
                    bound = d3_bound(Mred,p)
                    if bound[0] == True:
                        file.write(f"The j-invariant belongs to Fp, since D3 is between {p} and {bound[1]}.\n")
                        matrix_type = get_matrix_type(Mred, p)
                        file.write(f"Type {matrix_type[0]} matrix with parameters {matrix_type[1]}.\n\n")
                    else:
                        file.write(f"The j-invariant does not belong to Fp, since D3 is less than {bound[1]}.\n\n")

    
    if counter == 0:
        print(f"No prime number found.")
        file.write(f"No prime number found.\n")
