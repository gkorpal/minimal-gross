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
from functions.fp2.eichler.bounds import minimal_basis_matrix
from functions.fp2.deuring.correspondence import constructive_deuring


# Check if the correct number of arguments is provided
if len(sys.argv) != 4:
    print("Usage: sage fam2mod3.sage <c> <M> <N>")
    sys.exit(1)

c = abs(Integer(sys.argv[1])) # p = c mod 20 and 2 mod 3
M = abs(Integer(sys.argv[2])) # p >= M
N = abs(Integer(sys.argv[3])) # p <= N

with open(os.path.join(current_dir, f"fam2mod3_{c}_{M}_{N}.txt"), "w") as file:
    file.write(f"The j-invariants of all maximal orders with the first successive minima equal to 20 and lying in Bp for p = {c} (mod 20), 2 (mod 3) lying between {M} and {N}.\n\n")
    # loop over the primes c mod 20, 2 mod 3 and between [M,N]
    counter = 0
    for p in prime_range(M, N+1):
        if p % 20 == c and p % 3 == 2:
            file.write(f"{p = }.\n")
            counter+=1
            print(f"working on prime number {counter}...")
            # magma output [(a,b), base0, (base, matrix), ...]
            print(f"finding the maximal order using Magma...")
            m = get_maximal_order(-3, -p, 20)  
            a = m[0][0]
            b = m[0][1]
            B.<i,j,k> = QuaternionAlgebra(a, b)
            file.write(f"Quaternion Algebra: i^2 = {i^2} and j^2 = {j^2}.\n")
            assert p == B.discriminant(), "j(E)=0 does not belong to Bp"
            # the maximal order with j-invariant 0
            base0 = sage_eval(str(m[1]), locals={'i': i, 'j': j, 'k': k})
            O0 = B.quaternion_order(basis=base0)
            file.write(f"Maximal order with j(E)=0: {base0}.\n")

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

                    file.write(f"Rearranged Gram matrix with D1 <= D2 <= D3 and angles flipped if needed:\n{minimal_basis_matrix(Gred)}\n")
                    

                    # maximal orders are Eichler orders
                    O0O = O0.intersection(O)
                    N = O0O.discriminant()/B.discriminant()
                    I = N*O0*O 

                    ## Endomorphism ring ##
                    F.<alpha> = GF(p^2, 'alpha', modulus=x^2+3)
                    E0 = EllipticCurve(F,[0,1])
                    # https://math.stackexchange.com/questions/4370619/
                    # endomorphism iota corresponding to i
                    ker = E0(0,1)
                    phi = E0.isogeny(ker)
                    iso = phi.codomain().isomorphism_to(E0)
                    iota = iso * phi 

                    ## Deuring correspondence ##
                    while True:
                        try:
                            print("finding the corresponding j-invariant using Deuring for the People...")
                            E1, phi, _ = constructive_deuring(I, E0, iota)
                            newj = E1.j_invariant()
                            file.write(f"Corresponding j-invariants (alpha^2=-3): {newj}, {newj.conjugate()}.\nNorm and trace: {newj.norm()} and {newj.trace()}.\n\n")
                            print("")
                            break
                        except Exception as e:
                            #print(f"Error: {e}")  # Print out the error
                            print("will retry Deuring...")

    
    if counter == 0:
        print(f"No prime number found.")
        file.write(f"No prime number found.\n")
