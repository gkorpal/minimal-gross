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

from functions.fp.gross.numeric import NumericGross
from functions.fp.gross.type import get_matrix_type

p = 419 # >= 37

# computed using Kaneko bound
D1_bound = ceil(4*sqrt(p/3))

# Check if the correct number of arguments is provided
if len(sys.argv) != 1:
    print("Usage: sage csidh419gross.sage")
    sys.exit(1)

with open(os.path.join(current_dir, f"csidh419kaneko.txt"), "w") as file:
    file.write(f"Toy CSIDH prime p = 419 = 4*3*5*7-1\n\n")
    for D1 in range(1, D1_bound):
        file.write(f"**********\n")
        file.write(f"D1={D1}\n")
        D1 = Integer(D1)
        if  p%3 == 2 and D1 == 3:
            x = 1
            y = 1
            z = -(2*p-1)/3
            D2 = (4*p+1)/3
            D3 = (4*p+1)/3
            result = Matrix(ZZ, [[D1, x, y], [x, D2, z], [y, z, D3]])
            matrix_type = get_matrix_type(result, p)
            file.write(f"{result}\n")
            file.write(f"Type {matrix_type[0]} matrix with parameters {matrix_type[1]}.\n")
            file.write(f"Belongs to CSIDH graph.\n\n")
        else:
            ng = NumericGross(p)
            try:
                results = ng.gross(D1)
                if len(results) == 0:
                    file.write(f"No Gross lattice found.\n\n")
                else:
                    for result in results:
                        matrix_type = get_matrix_type(result, p)
                        file.write(f"{result}\n")
                        file.write(f"Type {matrix_type[0]} matrix with parameters {matrix_type[1]}.\n")
                        # we can directly tell from Gross lattice if the maximal order
                        # is of type O(q,r) or O'(q,r).
                        # only O(q,r) belong to CSIDH graph.
                        if result[0][0] != 4 and result[2][2] == p: 
                            file.write(f"Does not belong to CSIDH graph.\n\n")
                        else:
                            file.write(f"Belongs to CSIDH graph.\n\n")
            except Exception as e:
                # Code that runs if an exception occurs
                file.write(f"An error occurred: {e}\n\n")

        
