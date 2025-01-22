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

# CM data
# |disc| : [[a1,a2,...], N, jj, "j(E) = "]
# here p = ai mod |disc| is a supersingular prime family.
# N <= p < 37
cm_data = {
    3: [[2], 5, 0, "j(E) = 0"],
    4: [[3], 7, 12**3, "j(E) = 1728 = 12^3"],
    7: [[3, 5, 6], 13, -15**3, "j(E) = -15^3"],
    8: [[5, 7], 23, 20**3, "j(E) = 20^3"],
    11: [[2, 6, 7, 8, 10], 29, -32**3, "j(E) = -32^3"],
}

# Check if the correct number of arguments is provided
if len(sys.argv) != 2:
    print("Usage: sage cm_finite.sage <d>")
    sys.exit(1)

# Get the discriminant from the command line arguments
d = abs(Integer(sys.argv[1]))

# Check if the discriminant exists in the cm_data dictionary
if d not in cm_data:
    print(f"Error: enter discriminant of one of the d = 3,4,7,8,11.")
    sys.exit(1)

param = cm_data[d]
a_values = param[0]
N = param[1]
jj = param[2]


with open(os.path.join(current_dir, f"cases_{d}.txt"), "w") as file:
    file.write(f"{param[3]}\n\n")
    file.write(f"Obtaining Gram matrix of the minimal Gross lattice using LLL-reduction.\n\n")
    for p in primes(N,37):
        n = p%d
        if n in a_values:
            file.write(f"p = {p} is {n} mod {d}.\n")
            order_type, (q, r), order_base, OT_base, Gred, Mred = IbukiyamaMaximalOrder.get_matrix(jj, p)
            B.<i,j,k> = QuaternionAlgebra(-q,-p)
            if order_type == True:
                file.write(f"Maximal order: O({q},{r}) = {order_base}\n")
            else:
                file.write(f"Maximal order is O'({q},{r}) = {order_base}\n")
            
            file.write(f"Gross lattice: {OT_base}\n")
            file.write(f"LLL-reduced Gram matrix:\n{Gred}\n")
            file.write(f"Rearranged Gram matrix with D1 <= D2 <= D3 and angles flipped if needed:\n{Mred}\n")
            matrix_type, parameters = get_matrix_type(Mred, p)
            file.write(f"Type {matrix_type} matrix with parameters {parameters}.\n\n")
