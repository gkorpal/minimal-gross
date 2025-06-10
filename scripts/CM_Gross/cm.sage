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

from functions.fp.gross.symbolic import SymbolicGross
from functions.fp.gross.type import get_matrix_type

# CM data
# disc : [[a1,a2,...], b, D1, N, "j(E) = "]
# here p = ai mod b is a supersingular prime family.
cm_data = {
    3: [[2], 3, 3, 5, "j(E) = 0"],
    4: [[3], 4, 4, 7, "j(E) = 1728 = 12^3"],
    7: [[3, 5, 6], 7, 7, 13, "j(E) = -15^3"],
    8: [[5, 7], 8, 8, 23, "j(E) = 20^3"],
    11: [[2, 6, 7, 8, 10], 11, 11, 29, "j(E) = -32^3"],
    12: [[5, 11], 12, 12, 41, "j(E) = 2*30^3"],
    16: [[3, 7, 11, 15], 16, 16, 67, "j(E) = 66^3"],
    19: [[2, 3, 8, 10, 12, 13, 14, 15, 18], 19, 19, 79, "j(E) = -96^3"],
    27: [[2, 5, 8, 11, 14, 17, 20, 23, 26], 27, 27, 167, "j(E) = -3*160^3"],
    28: [[3, 5, 13, 17, 19, 27], 28, 28, 181, "j(E) = 255^3"],
    43: [[2, 3, 5, 7, 8, 12, 18, 19, 20, 22, 26, 27, 28, 29, 30, 32, 33, 34, 37, 39, 42], 43, 43, 433, "j(E) = -960^3"],
    67: [[2, 3, 5, 7, 8, 11, 12, 13, 18, 20, 27, 28, 30, 31, 32, 34, 38, 41, 42, 43, 44, 45, 46, 48, 50, 51, 52, 53, 57, 58, 61, 63, 66], 67, 67, 1103, "j(E) = -5280^3"],
    163: [[2, 3, 5, 7, 8, 11, 12, 13, 17, 18, 19, 20, 23, 27, 28, 29, 30, 31, 32, 37, 42, 44, 45, 48, 50, 52, 59, 63, 66, 67, 68, 70, 72, 73, 75, 76, 78, 79, 80, 82, 86, 89, 92, 94, 98, 99, 101, 102, 103, 105, 106, 107, 108, 109, 110, 112, 114, 116, 117, 120, 122, 123, 124, 125, 127, 128, 129, 130, 137, 138, 139, 141, 142, 147, 148, 149, 153, 154, 157, 159, 162], 163, 163, 6691, "j(E) = -640320^3"]
}

# Check if the correct number of arguments is provided
if len(sys.argv) != 2:
    print("Usage: sage cm.sage <d>")
    sys.exit(1)

# Get the discriminant from the command line arguments
d = abs(Integer(sys.argv[1]))

# Check if the discriminant exists in the cm_data dictionary
if d not in cm_data:
    print(f"Error: enter discriminant of one of the 13 CM curves over Q.")
    sys.exit(1)

param = cm_data[d]
a_values = param[0]
b = param[1]
D1 = param[2]
N = param[3]

# open a file called "cm_d.txt" in current_dir; overwrite exiting file with this name
with open(os.path.join(current_dir, f"cm_{d}.txt"), "w") as file:
    file.write(f"{param[4]}\n")
    file.write(f"N_E = {N}\n\n")
    for a in a_values:
        file.write(f"p = {a} mod {b}\n")
        if a == 2 and b == 3 and D1 == 3:
            p = var('p')
            x = 1
            y = 1
            z = -(2*p-1)/3
            D2 = (4*p+1)/3
            D3 = (4*p+1)/3
            result = matrix([[D1, x, y], [x, D2, z], [y, z, D3]])
            matrix_type = get_matrix_type(result)
        else:
            sg = SymbolicGross(a, b)
            result = sg.gross(D1, N)
            matrix_type = get_matrix_type(result)
        file.write(f"{result}\n")
        file.write(f"Type {matrix_type[0]} matrix with parameters {matrix_type[1]}.\n\n")
