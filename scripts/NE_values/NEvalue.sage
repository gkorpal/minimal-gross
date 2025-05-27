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

# CM data
# |disc| : [[a1,a2,...], jj, "j(E) = "]
# here p = ai mod |disc| is a supersingular prime family.
cm_data = {
    3: [[2], 0, "j(E) = 0"],
    4: [[3], 12**3, "j(E) = 1728 = 12^3"],
    7: [[3, 5, 6], -15**3, "j(E) = -15^3"],
    8: [[5, 7], 20**3, "j(E) = 20^3"],
    11: [[2, 6, 7, 8, 10], -32**3, "j(E) = -32^3"],
    12: [[5, 11], 2*30**3, "j(E) = 2*30^3"],
    16: [[3, 7, 11, 15], 66**3, "j(E) = 66^3"],
    19: [[2, 3, 8, 10, 12, 13, 14, 15, 18], -96**3, "j(E) = -96^3"],
    27: [[2, 5, 8, 11, 14, 17, 20, 23, 26], -3*160**3, "j(E) = -3*160^3"],
    28: [[3, 5, 13, 17, 19, 27], 255**3, "j(E) = 255^3"],
    43: [[2, 3, 5, 7, 8, 12, 18, 19, 20, 22, 26, 27, 28, 29, 30, 32, 33, 34, 37, 39, 42], -960**3, "j(E) = -960^3"],
    67: [[2, 3, 5, 7, 8, 11, 12, 13, 18, 20, 27, 28, 30, 31, 32, 34, 38, 41, 42, 43, 44, 45, 46, 48, 50, 51, 52, 53, 57, 58, 61, 63, 66], -5280**3, "j(E) = -5280^3"],
    163: [[2, 3, 5, 7, 8, 11, 12, 13, 17, 18, 19, 20, 23, 27, 28, 29, 30, 31, 32, 37, 42, 44, 45, 48, 50, 52, 59, 63, 66, 67, 68, 70, 72, 73, 75, 76, 78, 79, 80, 82, 86, 89, 92, 94, 98, 99, 101, 102, 103, 105, 106, 107, 108, 109, 110, 112, 114, 116, 117, 120, 122, 123, 124, 125, 127, 128, 129, 130, 137, 138, 139, 141, 142, 147, 148, 149, 153, 154, 157, 159, 162], -640320**3, "j(E) = -640320^3"]
}


# Check if the correct number of arguments is provided
if len(sys.argv) != 2:
    print("Usage: sage NEvalue.sage <d>")
    sys.exit(1)

# Get the discriminant from the command line arguments
d = abs(Integer(sys.argv[1]))

# Check if the discriminant exists in the cm_data dictionary
if d not in cm_data:
    print(f"Error: enter discriminant of one of the 13 CM curves over Q.")
    sys.exit(1)

param = cm_data[d]
a_values = param[0]
jj = param[1]

with open(os.path.join(current_dir, f"NE_{d}.txt"), "w") as file:
    file.write(f"{param[2]}\n\n")
    file.write(f"N_E is the smallest supersingular prime (greater than 3) such that for p >= N_E we get D1 = {d}.\n\n")

    p = 3
    counter = 0
    while True:
        p = next_probable_prime(p+1) # p+1 is always even; avoiding Pari related indexing issue.
        n = p%d
        if n in a_values:
            Mred = IbukiyamaMaximalOrder.get_matrix(jj, p)[-1]
            if Mred[0, 0] == d:
                counter += 1
                if counter == 1:
                    N = p
                elif counter == 10:
                    break
            else:
                counter = 0
    file.write(f"We get N_E = {N} by going through the list of supersingular primes and looking for a continuous sequence of 10 Eisenstein-reduced Gram matrices with D1 = {d}.\n")


