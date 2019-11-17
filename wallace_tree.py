import argparse
import math
from numpy import binary_repr
import random

from verilog_code import verilog_code

from wallace_add import wallace_add, wallace_simple_add


def calculate_max_length(M, N):     # page 5
    bit_number_floor = math.floor(math.log(M, 2))
    max_range = 2 ** bit_number_floor
    bit_number_floor_temp = bit_number_floor

    while bit_number_floor_temp - N >= 0:
        max_range += 2 ** (bit_number_floor_temp - N)
        bit_number_floor_temp -= N

    if M <= max_range:
        return N + bit_number_floor
    else:
        return N + bit_number_floor + 1


def wallace_adder(M, N):
    """This is the main function of the program, and it's responsible for generation of Verilog code."""
    max_length = calculate_max_length(M, N)  # number of bits result will have
    operands_len = M * N
    operands_decimal = []

    print("\nMax length is: ", max_length)

    verilog_code.append(
    f"""\
module wallace_tree(
    output reg [{max_length - 1}:0] result, 
    input [{M * N - 1}:0] initial_operands 
    );
    
    reg [{M * N - 1}:0] operands;       // Better name would be "intermediate_operands", but for simplicity's sake, it's "operands"
    reg [{M * N - 1}:0] operands_old;   // ...And this one has the values of the previous (the one before) iteration
    reg carry = 1'b0;
    
    initial begin
        operands_old = initial_operands;\
    """)
    # setting up test fixtures
    for m in range(M):
        operands_decimal.append(random.randrange(-2**(N-1), 2**(N-1)-1))

    # operands_decimal = [12, -4, 13, 7, -13, 9, 7, -13, -3, 9, -12, 9]
    # operands_decimal = [-12, 10, -13, 5]
    # operands_decimal = [-109, 105, 22, -66, 103, -40, -110, 8, 22, 117, -87, -125, 106, 22, 52, 21, 37, -85, 122, 3, 72, 52, 56, 25, 125, -106, 111, -105, 110, 108, -75, 0, 86, 120, -24, -62, 42, -122, 109, -79, 25, -85, 62, 116, -70, 48, 1, -49, -39, 82, -8, -53, 115, 114, 49, -101, -84, 51, 85, 104, -8, 112, 59, -3]
    # operands_decimal = [-34, -90, 23, 39, 113, -65, 107, -103, -123, -70, -90, -119, -63, 57, 41, 127, -44, 45, 82, -42, 16, 86, -2, 94, -52, 48, 16, -88, -67, -128, -93, 26, 6, -66, -21, 16, 78, -73, -90, 93, 64, -93, -78, -74, 37, -13, 110, -54, 66, 94, 42, 76, 56, 9, 37, 73, 11, 87, 24, 7, 88, -94, 30, -123]



    print(operands_decimal)

    operands_verilog = "".join([binary_repr(operand_decimal, N) for operand_decimal in operands_decimal]) #[::-1]]
    operands = [binary_repr(operand_decimal, N)[::-1] for operand_decimal in operands_decimal]
    operands_flat = list(map(int, list(''.join(operands))))

    print(operands_verilog)

    operands = operands_flat

    iter = 1
    (operands, n, operands_len, mut_one_len_new, mut_one_calib_new, mut_two_len_new, mut_two_calib_new) = wallace_add(operands, M, N, iter, max_length)
    del operands[operands_len:]

    verilog_code.append(
        f"""
        operands_old = operands;\
        """)

    while n >= 3:
        iter += 1
        (operands, n, operands_len, mut_one_len_new, mut_one_calib_new, mut_two_len_new, mut_two_calib_new) = wallace_add(operands, n, min(max_length, N - 1 + iter), iter, max_length, mut_one_len_new, mut_one_calib_new, mut_two_len_new, mut_two_calib_new)
        del operands[operands_len:]
        verilog_code.append(
            f"""
        operands_old = operands;\
            """)


    if (mut_two_len_new):
        result = wallace_simple_add(operands, max_length, mut_one_len_new, mut_one_calib_new, mut_two_len_new, mut_two_calib_new)
    else:
        result = wallace_simple_add(operands, max_length, min(max_length, N - 1 + iter + 1), 0, mut_one_len_new, mut_one_calib_new)

    verilog_code.append(
    f"""
    end
endmodule
    """)
    #print(operands_decimal)
    real_result = sum(operands_decimal)
    #print("\nMax length is: ", max_length)
    if result != real_result:
        raise Exception("It doesn't work.")
    print("\n\nCorrect result: ", real_result, "\nCalculated one: ", result)



    with open("wallace_tree.txt", "w") as outfile:
        outfile.write("\n".join(verilog_code))

    # print(operands_decimal)
    # print(operands)
    # print(operands_flat)



def build_arg_parser():
    parser = argparse.ArgumentParser(description='Generates Verilog code for Wallace tree of \
                3-2 carry save adders using specified number of operands (M) and number of bits (N).')
    parser.add_argument("M", help="Specify the number of operands")
    parser.add_argument("N", help="What's the bit-length of these \
                                                           operands")
    return parser


if __name__ == '__main__':
    args = build_arg_parser().parse_args()
    M = int(args.M)
    N = int(args.N)
    # for M in range(5, 200):
    #     for N in range(2, 30):
    wallace_adder(M, N)

