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

from functions.fp.ibukiyama.maximal import IbukiyamaMaximalOrder
from functions.fp.gross.type import get_matrix_type 


# Check if the correct number of arguments is provided
if len(sys.argv) != 1:
    print("Usage: sage csidh419.sage")
    sys.exit(1)

p = 419
F = GF(p)
# disc : [F(jj), "j = jj", N] 
cm_data = {
    3: [F(0), "j = 0", 5],
    4: [F(12^3), "j = 12^3", 7],
    7: [F(-15^3), "j = -15^3", 13],
    8: [F(20^3), "j = 20^3", 23],
    11: [F(-32^3), "j = -32^3", 29],
    12: [F(2*30^3), "j = 2*30^3", 41],
    16: [F(66^3), "j = 66^3", 67],
    19: [F(-96^3), "j = -96^3", 79],
    27: [F(-3*160^3), "j = -3*160^3", 167],
    28: [F(255^3), "j = 255^3", 181],
    43: [F(-960^3), "j = -960^3", 433],
    67: [F(-5280^3), "j = -5280^3", 1103],
    163: [F(-640320^3), "j = -640320^3", 6691]
}

SSjp2 = SupersingularModule(p).supersingular_points()[0]
SSj = Set(SSjp2).intersection(Set(F)).symmetric_difference(Set([F(1728)]))
ssj = list(SSj)
D1_values = []

with open(os.path.join(current_dir, f"csidh419.txt"), "w") as file:
    file.write(f"Toy CSIDH prime p = 419 = 4*3*5*7-1\n\n")
    jj = F(1728) # special case; p = 3 (mod 8)
    # E : y^2 = x^3 + 0*x + x
    file.write(f"j-invariant {jj} leads to Montgomery curve with label 0.\n")
    
    file.write(f"**This is one of the 13 CM curves over Q. {cm_data[4][1]} with D1 = 4.**\n")

    B.<i,j,k> = QuaternionAlgebra(-1, -p)
    order_base = [1, i, (1+k)/2, (i-j)/2] # minimal basis; see the paper.
    file.write(f"Maximal order: O({1}) = {order_base}\n")
    
    OT_base = [2*aa - aa.reduced_trace() for aa in order_base[1:]]
    file.write(f"Gross lattice: {OT_base}\n")
    
    G = IbukiyamaMaximalOrder.gram_matrix(OT_base)
    U = G.LLL_gram()
    Gred = U.transpose() * G * U
    file.write(f"LLL-reduced Gram matrix:\n{Gred}\n")
    
    Mred = IbukiyamaMaximalOrder.minimal_basis_matrix(Gred)
    file.write(f"Rearranged Gram matrix with D1 <= D2 <= D3 and angles flipped if needed:\n{Mred}\n")
    D1_values.append(Mred[0, 0])
    
    matrix_type, parameters = get_matrix_type(Mred, p)
    file.write(f"Type {matrix_type} matrix with parameters {parameters}.\n\n")

    for jj in ssj:
        E = EllipticCurve_from_j(F(jj))
        try:
            E = E.montgomery_model()
            A = E.a2()
            file.write(f"j-invariant {jj} leads to Montgomery curves with labels {A} and {p-A}.\n")

            d_values = [key for key in cm_data if cm_data[key][0] == jj]
            d = min(d_values) if d_values else None
            if d and p >= cm_data[d][2]:
                file.write(f"**This is one of the 13 CM curves over Q. {cm_data[d][1]} with D1 = {d} since p >= {cm_data[d][2]}.**\n")
            elif d and p < cm_data[d][2]:
                file.write(f"**This is one of the 13 CM curves over Q. {cm_data[d][1]} with D1 <= {d} since p < {cm_data[d][2]}.**\n")

            order_type, (q, r), order_base, OT_base, Gred, Mred = IbukiyamaMaximalOrder.get_matrix(jj, p)
            
            B.<i,j,k> = QuaternionAlgebra(-q,-p)
            if order_type == True:
                file.write(f"Maximal order: O({q},{r}) = {order_base}\n")
            else:
                file.write(f"Maximal order is O'({q},{r}) = {order_base}\n")
            
            file.write(f"Gross lattice: {OT_base}\n")
            
            file.write(f"LLL-reduced Gram matrix:\n{Gred}\n")
            
            file.write(f"Rearranged Gram matrix with D1 <= D2 <= D3 and angles flipped if needed:\n{Mred}\n")
            D1_values.append(Mred[0, 0])

            matrix_type, parameters = get_matrix_type(Mred, p)
            file.write(f"Type {matrix_type} matrix with parameters {parameters}.\n\n")
        except ValueError:
            pass
    file.write(f"D1 values: {sorted(D1_values)}.")
