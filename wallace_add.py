from numpy import binary_repr
import numpy as np
from pandas import *
from matrix_sum import matrix_sum
from verilog_code import verilog_code

split_display_ind: bool


def bin_format(result):
    return list(map(int, list(binary_repr(result, 2))))


def print_state(operands, n, length, mut_count, mut_one_calib, mut_two_calib, iter):
    split_array = get_split(iter, n, length, mut_count, mut_one_calib, mut_two_calib)

    global split_display_ind
    if split_display_ind:
        print("Split ARRAY", iter, ": ", split_array)
        split_display_ind = False

    oper1 = np.split(np.array(operands), split_array)
    oper1 = [oper[::-1] for oper in oper1]
    oper1 = oper1[:-1]
    a = DataFrame(oper1, columns=(length + 1 - saturation) * [''], index=n * [''])
    a = a.fillna('')
    print(a)
    print("The sum is: ", matrix_sum(oper1, n, length + 1 - saturation))



def get_split(iter, n_new, length, mut_count_new, mut_one_calib_new, mut_two_calib_new):
    if iter == 1:
        split_array = [length + i for i in [1, 0] * (int)((n_new - mut_count_new) / 2)]
    else:
        split_array_one = [length + 1 - saturation for _ in range(0, n_new - mut_count_new, 2)]
        split_array_two = [length - 1 - saturation + (i % 2) for _, i in zip(range(1, n_new - mut_count_new, 2), range(1, n_new - mut_count_new))]
        split_array = [item for pair in zip(split_array_one, split_array_two) for item in pair]

    if (n_new - mut_count_new) % 2 == 1:
        split_array.append(length + 1 - saturation)
    if mut_one_calib_new is not None:
        split_array += [length + 1 - saturation - mut_one_calib_new]
    if mut_two_calib_new is not None:
        split_array += [length + 1 - saturation -  mut_two_calib_new]
    return [sum(split_array[0:i]) for i in range(1, len(split_array) + 1)]


def wallace_add(operands_old, n, length, iter, max_length, mut_one_len=0, mut_one_calib=0, mut_two_len=0, mut_two_calib=0):
    """This function is responsible for one iteration of the wallace addition. Generates the logic,
    and also does the simulation

    Args:
        operands_old (list)(int): concatenated binary representations of our operands from the previous iteration
        n (number): number of distinct operands in this iteration
        length (number): the biggest length of an operand in this iteration
        iter (number): iteration number which is used for calibration

    Returns:
        operands_len (number): number of operands left
    """
    operands = [0] * len(operands_old)  # result of one iteration

    remainder = n % 3  # remainder is used to determine mutation process
    calibrator = 0  # calibrators are those little offsets in the operands
    calibrator_two = 0
    calibrator_sum = 0  # calibrator sums are used for calculating our position inside operands array(string)
    calibrator_sum_old = 0
    g = 0
    element_two_index = 0
    index_modulus = 0

    # mutation variables track mutated operands
    if mut_one_len and mut_two_len:
        mut_count = 2
    elif mut_one_len or mut_two_len:
        mut_count = 1
    else:
        mut_count = 0

    # these two variables are used for finding the calibrator value
    if (iter == 1):
        modulus = 1
        offset_base = 0
    elif (iter == 2):
        modulus = 1
        offset_base = 1
    else:
        modulus = 2
        offset_base = 1

    healthy_n = (n - remainder) if (remainder >= mut_count) else (n - remainder - 3)
    health_plus_one = 0         # do we have the healthy operand (as a result) inside the mutation part 2-1 or 1-0 combo
    global saturation
    saturation = 1 if length == max_length else 0

    #TODO this print is fine
    #print("Iteration ", iter)

    """
    HEALTHY PART
    """
    global split_display_ind
    split_display_ind = True

    for i in range(length):
        ind = False  # type of the group
        g = -1  # counts the group
        index_modulus = -1  # generates correct calibrator

        calibrator_sum = 0
        calibrator_sum_old = 0
        for k in range(0, healthy_n, 3):
            ind = not ind
            g += 1

            if ind:
                index_modulus += 1
                calibrator = offset_base + (index_modulus % modulus)
                calibrator_sum_old = calibrator_sum
                calibrator_sum += calibrator

                if k >= 6 and iter == 1:
                    element_one_index = int(g * 2 * (length + 1 - saturation) - g + i)
                else:
                    element_one_index = int(g * 2 * (length + 1 - saturation) - 1.5 * g + i)

                element_two_index = element_one_index + (length + 1 - saturation)

                if i < calibrator:
                    operands[element_two_index], operands[element_one_index] = bin_format(operands_old[k * length - calibrator_sum_old + i] + \
                                                                                          operands_old[(k + 2) * length - calibrator_sum + i])
                    verilog_code.append(
                    f"""\
        {{operands[{element_two_index}], operands[{element_one_index}]}} = operands_old[{k * length - calibrator_sum_old + i}] + operands_old[{(k + 2) * length - calibrator_sum + i}];\
                    """)
                elif i >= calibrator:
                    if (i != length - 1) or (not saturation):
                        operands[element_two_index], operands[element_one_index] = bin_format(operands_old[k * length - calibrator_sum_old + i] + \
                                                                                              operands_old[(k + 1) * length - calibrator_sum_old + (i - calibrator)] + \
                                                                                              operands_old[(k + 2) * length - calibrator_sum + i])
                        verilog_code.append(
                        f"""\
        {{operands[{element_two_index}], operands[{element_one_index}]}} = operands_old[{k * length - calibrator_sum_old + i}] + operands_old[{(k + 1) * length - calibrator_sum_old + (i - calibrator)}] + operands_old[{(k + 2) * length - calibrator_sum + i}];\
                        """)
                    else:
                        operands[element_one_index] = operands_old[k * length - calibrator_sum_old + i] ^ \
                                                      operands_old[(k + 1) * length - calibrator_sum_old + (i - calibrator)] ^ \
                                                      operands_old[(k + 2) * length - calibrator_sum + i]
                        verilog_code.append(
                        f"""\
        operands[{element_one_index}] = operands_old[{k * length - calibrator_sum_old + i}] ^ operands_old[{(k + 1) * length - calibrator_sum_old + (i - calibrator)}] ^ operands_old[{(k + 2) * length - calibrator_sum + i}];\
                        """)
                        element_two_index -= 1

            else:
                index_modulus += 1
                calibrator = offset_base + (index_modulus % modulus)
                index_modulus += 1
                calibrator_two = offset_base + (index_modulus % modulus)

                calibrator_sum_old = calibrator_sum
                calibrator_sum += calibrator

                if g == 1:
                    # 2 * (length + 1) -1 = 2 * length + 1
                    # 2 * (length + 1 - saturation) - 1
                    element_one_index = 2 * (length + 1 - saturation) - 1 + i
                else:
                    if k >= 6 and iter == 1:
                        element_one_index = int(g * 2 * (length + 1 - saturation) - g + i)
                    else:
                        element_one_index = int(g * 2 * (length + 1 - saturation) - 1.5 * (g - 1) - 1 + i)
                element_two_index = element_one_index + (length + 1 - saturation) - min(calibrator_two, calibrator)
                # Explanation for this above... "min(calibrator_two, calibrator)" is used to deduct from "i" so it would start at "element_one_index + (length + 1 - saturation) + 0"
                # Another system is being used for the mutated region

                if (i == 0) and (iter == 1):
                    operands[element_two_index], operands[element_one_index] = bin_format(operands_old[k * length - calibrator_sum_old + (i - calibrator)] + \
                                                                                          operands_old[(k + 1) * length - calibrator_sum + i] + \
                                                                                          operands_old[(k + 2) * length - calibrator_sum + (i - calibrator_two)])
                    verilog_code.append(
                    f"""\
        {{operands[{element_two_index}], operands[{element_one_index}]}} = operands_old[{k * length - calibrator_sum_old + (i - calibrator)}] + operands_old[{(k + 1) * length - calibrator_sum + i}] + operands_old[{(k + 2) * length - calibrator_sum + (i - calibrator_two)}];\
                    """)
                elif (i == 0):
                    operands[element_one_index] = operands_old[(k + 1) * length - calibrator_sum]
                    verilog_code.append(
                    f"""\
        operands[{element_one_index}] = operands_old[{(k + 1) * length - calibrator_sum}];\
                    """)
                elif (i == length - 1) and saturation:
                    operands[element_one_index] = operands_old[k * length - calibrator_sum_old + (i - calibrator)] ^ \
                                                             operands_old[(k + 1) * length - calibrator_sum + i] ^ \
                                                             operands_old[(k + 2) * length - calibrator_sum + (i - calibrator_two)]
                    verilog_code.append(
                    f"""\
        operands[{element_one_index}] = operands_old[{k * length - calibrator_sum_old + (i - calibrator)}] ^ operands_old[{(k + 1) * length - calibrator_sum + i}] ^ operands_old[{(k + 2) * length - calibrator_sum + (i - calibrator_two)}];\
                    """)
                    element_two_index -= 1  ## These are neccesarry if there is no interpolation to the left for the 'upper' operand
                elif i >= calibrator and i >= calibrator_two:
                    operands[element_two_index], operands[element_one_index] = bin_format(operands_old[k * length - calibrator_sum_old + (i - calibrator)] + \
                                                                                          operands_old[(k + 1) * length - calibrator_sum + i] + \
                                                                                          operands_old[(k + 2) * length - calibrator_sum + (i - calibrator_two)])
                    verilog_code.append(
                    f"""\
        {{operands[{element_two_index}], operands[{element_one_index}]}} = operands_old[{k * length - calibrator_sum_old + (i - calibrator)}] + operands_old[{(k + 1) * length - calibrator_sum + i}] + operands_old[{(k + 2) * length - calibrator_sum + (i - calibrator_two)}];\
                    """)
                elif calibrator < calibrator_two:
                    operands[element_two_index], operands[element_one_index] = bin_format(operands_old[k * length - calibrator_sum_old + (i - calibrator)] + \
                                                                                          operands_old[(k + 1) * length - calibrator_sum + i])
                    verilog_code.append(
                    f"""\
        {{operands[{element_two_index}], operands[{element_one_index}]}} = operands_old[{k * length - calibrator_sum_old + (i - calibrator)}] + operands_old[{(k + 1) * length - calibrator_sum + i}];\
                    """)
                elif calibrator_two < calibrator:
                    operands[element_two_index], operands[element_one_index] = bin_format(operands_old[(k + 1) * length - calibrator_sum + i] + \
                                                                                          operands_old[(k + 2) * length - calibrator_sum + (i - calibrator_two)])
                    verilog_code.append(
                    f"""\
        {{operands[{element_two_index}], operands[{element_one_index}]}} = operands_old[{(k + 1) * length - calibrator_sum + i}] + operands_old[{(k + 2) * length - calibrator_sum + (i - calibrator_two)}];\
                    """)
                calibrator_sum += calibrator_two

            if (i == length - 1) and (not saturation):
                operands[element_one_index + 1] = operands[element_one_index]
                verilog_code.append(
                f"""\
        operands[{element_one_index + 1}] = operands[{element_one_index}];\
                """)

            #TODO this print is fine
            #print_state(operands, (int)(healthy_n * 2 / 3), length, 0, None, None, iter)


    """
    MUTATED PART
    """
    normal_length = healthy_n * length - calibrator_sum
    split_display_ind = True

    if remainder < mut_count:
        """
        Complicated mutation
        """

        ind = not ind
        index_modulus += 1
        calibrator = offset_base + (index_modulus % modulus)

        if healthy_n != 0:
            element_one_index = element_two_index + 1
        else:
            element_one_index = 0

        if mut_count == 2 and remainder == 0:  # ind = 0
            if ind:  # this one is handling the ind = 1 (which can happen when you have had 2-1 previously)
                calibrator = 0
            mut_one_len_new = length - min(calibrator, mut_one_calib, mut_two_calib)
            element_two_index = element_one_index + mut_one_len_new

            if (length - (mut_one_calib + mut_one_len)):
                verilog_code.append(
                f"""\
        operands_old = {{{{operands_old[{len(operands_old) - 1}:{len(operands_old) - mut_two_len}]}}, {{{length - (mut_one_calib + mut_one_len)}{{operands_old[{len(operands_old) - mut_two_len - 1}]}}}}, {{operands_old[{len(operands_old) - mut_two_len - 1}:0]}}}};\
                """)

            operands_old[-mut_two_len:-mut_two_len] = [operands_old[-mut_two_len - 1]] * (length - (mut_one_calib + mut_one_len))

            if (length - (mut_two_calib + mut_two_len)):
                verilog_code.append(
                f"""\
            operands_old = {{{{ {length - (mut_two_calib + mut_two_len)}{{operands_old[{len(operands_old) - 1}]}} }}, operands_old[{len(operands_old)-1}:0]}};\
                """)

            operands_old += [operands_old[-1]] * (length - (mut_two_calib + mut_two_len))



            #print("Extended operands_old")
            # split_array = get_split(iter, n, length - 1, mut_count, mut_one_calib, mut_two_calib)
            #             # print("Extended array ", split_array)
            #             # oper1 = np.split(np.array(operands_old), split_array)
            #             # oper1 = [oper[::-1] for oper in oper1]
            #             # oper1 = oper1[:-1]
            #             # a = DataFrame(oper1, columns=length * [''], index=n * [''])
            #             # a = a.fillna('')
            #             # print(a)

            ### #print_state(operands_old, n, length - 1, mut_count, mut_one_calib, mut_two_calib, iter)

            for i in range(length):

                if i >= calibrator and i >= mut_one_calib and i >= mut_two_calib:
                    try:
                        if (i != length - 1) or (not saturation):
                            operands[element_two_index], operands[element_one_index] = bin_format(operands_old[normal_length + i - calibrator] + \
                                                                                                  operands_old[normal_length + length + i - calibrator - mut_one_calib] + \
                                                                                                  operands_old[normal_length + length + length + i - calibrator - mut_one_calib - mut_two_calib])
                            verilog_code.append(
                                f"""\
        {{operands[{element_two_index}], operands[{element_one_index}]}} = operands_old[{normal_length + i - calibrator}] + operands_old[{normal_length + length + i - calibrator - mut_one_calib}] + operands_old[{normal_length + length + length + i - calibrator - mut_one_calib - mut_two_calib}];\
                                """)
                        else:
                            operands[element_one_index] = operands_old[normal_length + i - calibrator] ^ \
                                                          operands_old[normal_length + length + i - calibrator - mut_one_calib] ^ \
                                                          operands_old[normal_length + length + length + i - calibrator - mut_one_calib - mut_two_calib]
                            verilog_code.append(
                            f"""\
        operands[{element_one_index}] = operands_old[{normal_length + i - calibrator}] ^ operands_old[{normal_length + length + i - calibrator - mut_one_calib}] ^ operands_old[{normal_length + length + length + i - calibrator - mut_one_calib - mut_two_calib}];\
                            """)
                            element_two_index -= 1
                        # operands_old[normal_length + length + mut_one_len + i - calibrator - mut_one_calib - mut_two_calib]
                    except:
                        print("error")
                    element_one_index += 1
                    element_two_index += 1
                elif i >= calibrator and i >= mut_one_calib:
                    operands[element_two_index], operands[element_one_index] = bin_format(operands_old[normal_length + i - calibrator] + \
                                                                                          operands_old[normal_length + length + i - calibrator - mut_one_calib])
                    verilog_code.append(
                        f"""\
        {{operands[{element_two_index}], operands[{element_one_index}]}} = operands_old[{normal_length + i - calibrator}] + operands_old[{normal_length + length + i - calibrator - mut_one_calib}];\
                        """)
                    element_one_index += 1
                    element_two_index += 1
                elif i >= calibrator and i >= mut_two_calib:
                    operands[element_two_index], operands[element_one_index] = bin_format(operands_old[normal_length + i - calibrator] + \
                                                                                          operands_old[normal_length + length + length + i - calibrator - mut_one_calib - mut_two_calib])
                    verilog_code.append(
                    f"""\
        {{operands[{element_two_index}], operands[{element_one_index}]}} = operands_old[{normal_length + i - calibrator}] + operands_old[{normal_length + length + length + i - calibrator - mut_one_calib - mut_two_calib}];\
                    """)
                    element_one_index += 1
                    element_two_index += 1
                elif i >= mut_one_calib and i >= mut_two_calib:
                    operands[element_two_index], operands[element_one_index] = bin_format(operands_old[normal_length + length + i - calibrator - mut_one_calib] + \
                                                                                          operands_old[normal_length + length + length + i - calibrator - mut_one_calib - mut_two_calib])
                    verilog_code.append(
                        f"""\
        {{operands[{element_two_index}], operands[{element_one_index}]}} = operands_old[{normal_length + length + i - calibrator - mut_one_calib}] + operands_old[{normal_length + length + length + i - calibrator - mut_one_calib - mut_two_calib}];\
                        """)
                    element_one_index += 1
                    element_two_index += 1
                elif i >= calibrator:
                    operands[element_one_index] = operands_old[normal_length + i - calibrator]
                    verilog_code.append(
                        f"""\
        operands[{element_one_index}] = operands_old[{normal_length + i - calibrator}];\
                        """)
                    element_one_index += 1
                elif i >= mut_one_calib:
                    operands[element_one_index] = operands_old[normal_length + length + i - calibrator - mut_one_calib]
                    verilog_code.append(
                        f"""\
        operands[{element_one_index}] = operands_old[{normal_length + length + i - calibrator - mut_one_calib}];\
                        """)
                    element_one_index += 1
                elif i >= mut_two_calib:
                    operands[element_one_index] = operands_old[normal_length + length + length + i - calibrator - mut_one_calib - mut_two_calib]
                    verilog_code.append(
                        f"""\
        operands[{element_one_index}] = operands_old[{normal_length + length + length + i - calibrator - mut_one_calib - mut_two_calib}];\
                        """)
                    element_one_index += 1

                # oper1 = np.split(np.array(operands), [8, 15, 23, 29])
                # oper1 = [oper[::-1] for oper in oper1]
                # oper1 = oper1[:-1]
                # a = DataFrame(oper1, columns=8 * [''], index=4 * [''])
                # a = a.fillna('')
                # print(a)

            mut_one_calib_new = min(calibrator, mut_one_calib, mut_two_calib)
            temp = [calibrator, mut_one_calib, mut_two_calib]
            temp.sort()
            mut_two_calib_new = temp[1] + 1
            mut_two_len_new = length + 1 - saturation - mut_two_calib_new

            # if length < max_length:
            #     operands[element_one_index] =  operands[element_one_index - 1]
            # oper1 = np.split(np.array(operands), [8, 15, 23, 29])
            # oper1 = [oper[::-1] for oper in oper1]
            # oper1 = oper1[:-1]
            # a = DataFrame(oper1, columns=8 * [''], index=4 * [''])
            # a = a.fillna('')
            # print(a)
        elif (mut_count == 1 and remainder == 0) or (mut_count == 2 and remainder == 1):
            mut_one_len_new = length - min(calibrator, mut_one_calib) - saturation     # only one mutagen in *here*
                                                                                       # and saturation is here for cutoff
            element_two_index = element_one_index + length + 1 - saturation            # interpolate to the left because of this length + 1

            health_plus_one = 1
            # Extending...
            if mut_count == 2:
                if length - (mut_one_calib + mut_one_len):
                    verilog_code.append(
                    f"""\
            operands_old = {{{{operands_old[{len(operands_old) - 1}:{len(operands_old) - mut_two_len}]}}, {{{length - (mut_one_calib + mut_one_len)}{{operands_old[{len(operands_old) -mut_two_len - 1}]}}}}, {{operands_old[{len(operands_old) - mut_two_len - 1}:0]}}}};\
                    """)
                operands_old[-mut_two_len:-mut_two_len] = [operands_old[-mut_two_len - 1]] * (length - (mut_one_calib + mut_one_len))
            else:
                if length - (mut_one_calib + mut_one_len):
                    verilog_code.append(
                        f"""\
            operands_old = {{{{ {length - (mut_one_calib + mut_one_len)}{{operands_old[{len(operands_old) - 1}]}} }}, operands_old[{len(operands_old) - 1}: 0]}};\
                        """)
                operands_old += [operands_old[-1]] * (length - (mut_one_calib + mut_one_len))

            #print("Extended operands_old")
            ### #print_state(operands_old + (length - mut_two_len - mut_two_calib) * [None] , n, length - 1, mut_count, mut_one_calib, mut_two_calib, iter)

            if ind:
                for i in range(length):
                    if i >= calibrator and i >= mut_one_calib:
                        if (i != length - 1) or (not saturation):
                            operands[element_two_index], operands[element_one_index] = bin_format(operands_old[normal_length + i] + \
                                                                                                  operands_old[normal_length + length + i - calibrator] + \
                                                                                                  operands_old[normal_length + 2 * length + i - calibrator - mut_one_calib])
                            verilog_code.append(
                                f"""\
        {{operands[{element_two_index}], operands[{element_one_index}]}} = operands_old[{normal_length + i}] + operands_old[{normal_length + length + i - calibrator}] + operands_old[{normal_length + 2 * length + i - calibrator - mut_one_calib}];\
                                """)
                        else:
                            operands[element_one_index] = operands_old[normal_length + i] ^ \
                                                          operands_old[normal_length + length + i - calibrator] ^ \
                                                          operands_old[normal_length + 2 * length + i - calibrator - mut_one_calib]
                            verilog_code.append(
                                f"""\
        operands[{element_one_index}] = operands_old[{normal_length + i}] ^ operands_old[{normal_length + length + i - calibrator}] ^ operands_old[{normal_length + 2 * length + i - calibrator - mut_one_calib}];\
                                """)
                            element_two_index -= 1  # so it would negate the following addition...
                    elif i >= calibrator:
                        operands[element_two_index], operands[element_one_index] = bin_format(operands_old[normal_length + i] + \
                                                                                              operands_old[normal_length + length + i - calibrator])
                        verilog_code.append(
                            f"""\
        {{operands[{element_two_index}], operands[{element_one_index}]}} = operands_old[{normal_length + i}] + operands_old[{normal_length + length + i - calibrator}];\
                            """)
                    elif i >= mut_one_calib:
                        operands[element_two_index], operands[element_one_index] = bin_format(operands_old[normal_length + i] + \
                                                                                              operands_old[normal_length + 2 * length + i - calibrator - mut_one_calib])
                        verilog_code.append(
                            f"""\
        {{operands[{element_two_index}], operands[{element_one_index}]}} = operands_old[{normal_length + i}] + operands_old[{normal_length + 2 * length + i - calibrator - mut_one_calib}];\
                            """)
                    else:
                        operands[element_one_index] = operands_old[normal_length + i]
                        verilog_code.append(
                            f"""\
        operands[{element_one_index}] = operands_old[{normal_length + i}];\
                            """)
                        element_two_index -= 1  # so it would negate the following addition...
                    element_one_index += 1
                    element_two_index += 1
            else:  # ind = 0
                for i in range(length):
                    if i >= calibrator and i >= mut_one_calib:
                        if (i != length - 1) or (not saturation):

                            operands[element_two_index], operands[element_one_index] = bin_format(operands_old[normal_length + i - calibrator] + \
                                                                                                  operands_old[normal_length + length + i - calibrator] + \
                                                                                                  operands_old[normal_length + 2 * length + i - calibrator - mut_one_calib])
                            verilog_code.append(
                                f"""\
        {{operands[{element_two_index}], operands[{element_one_index}]}} = operands_old[{normal_length + i - calibrator}] + operands_old[{normal_length + length + i - calibrator}] + operands_old[{normal_length + 2 * length + i - calibrator - mut_one_calib}];\
                                """)
                        else:
                            operands[element_one_index] = operands_old[normal_length + i - calibrator] ^ \
                                                          operands_old[normal_length + length + i - calibrator] ^ \
                                                          operands_old[normal_length + 2 * length + i - calibrator - mut_one_calib]
                            verilog_code.append(
                                f"""\
        operands[{element_one_index}] = operands_old[{normal_length + i - calibrator}] ^ operands_old[{normal_length + length + i - calibrator}] ^ operands_old[{normal_length + 2 * length + i - calibrator - mut_one_calib}];\
                                """)
                            element_two_index -= 1

                    elif i >= calibrator:
                        operands[element_two_index], operands[element_one_index] = bin_format(operands_old[normal_length + i - calibrator] + \
                                                                                              operands_old[normal_length + length + i - calibrator])
                        verilog_code.append(
                            f"""\
        {{operands[{element_two_index}], operands[{element_one_index}]}} = operands_old[{normal_length + i - calibrator}] + operands_old[{normal_length + length + i - calibrator}];\
                            """)
                    elif i >= mut_one_calib:
                        operands[element_two_index], operands[element_one_index] = bin_format(operands_old[normal_length + length + i - calibrator] + \
                                                                                              operands_old[normal_length + 2 * length + i - calibrator - mut_one_calib])
                        verilog_code.append(
                            f"""\
        {{operands[{element_two_index}], operands[{element_one_index}]}} = operands_old[{normal_length + length + i - calibrator}] + operands_old[{normal_length + 2 * length + i - calibrator - mut_one_calib}];\
                            """)
                    else:
                        operands[element_one_index] = operands_old[normal_length + length + i - calibrator]
                        verilog_code.append(
                            f"""\
        operands[{element_one_index}] = operands_old[{normal_length + length + i - calibrator}];\
                            """)
                        element_two_index -= 1  # so it would negate the following addition...
                    element_one_index += 1
                    element_two_index += 1

            if (not saturation):
                operands[element_one_index] = operands[element_one_index - 1]       # interpolation finished
                verilog_code.append(
                    f"""\
        operands[{element_one_index}] = operands[{element_one_index - 1}];\
                    """)
            mut_one_calib_new = min(calibrator, mut_one_calib) + 1

            if remainder == 1:
                for i in range(mut_two_len):
                    operands[element_two_index + i] = operands_old[-mut_two_len + i]

                verilog_code.append(
                    f"""\
        operands[{element_two_index + mut_two_len - 1}:{element_two_index}] = operands_old[{len(operands_old) - 1}:{len(operands_old) - mut_two_len}];\
                    """)
                # This one handles the 4th iteration (2 healthy + 2 mut)
                mut_two_len_new = mut_two_len
                mut_two_calib_new = mut_two_calib
                element_two_index += mut_two_len
                # if n == 4:
            else:
                mut_two_len_new = 0
                mut_two_calib_new = None

        #print("Length is: ", element_two_index)

        operands_len = element_two_index
    elif remainder >= mut_count:
        """
        Simple mutation
        """

        element_two_index += 1  # used for copy-pasting values from operands_old to operands
        index_modulus += 1
        calibrator = offset_base + (index_modulus % modulus)
        ind = not ind

        if mut_count == 0:
            for i in range(len(operands_old) - normal_length):
                operands[element_two_index + i] = operands_old[normal_length + i]

            if (len(operands_old) > normal_length):
                verilog_code.append(
                    f"""\
        operands[{len(operands_old) - normal_length + element_two_index - 1}:{element_two_index}] = operands_old[{len(operands_old) - 1}:{normal_length}];\
                    """)
            if remainder == 1:
                if ind:
                    if iter > 1:
                        # probably never happens
                        mut_one_len_new = length - calibrator
                        mut_one_calib_new = calibrator
                    else:
                        mut_one_len_new = length
                        mut_one_calib_new = 0

                else:
                    # if iter > 2:
                    #     mut_one_len_new = length - 2
                    #     mut_one_calib_new = 2
                    # elif iter == 2:
                    #     mut_one_len_new = length - 1
                    #     mut_one_calib_new = 1
                    # else:
                    #     mut_one_len_new = length
                    #     mut_one_calib_new = 0
                    if iter == 1:
                        mut_one_len_new = length
                        mut_one_calib_new = 0
                    else:
                        mut_one_len_new = length - calibrator
                        mut_one_calib_new = calibrator

                mut_two_len_new = 0
                mut_two_calib_new = None
            elif remainder == 2:
                if ind:
                    mut_one_len_new = length
                    mut_one_calib_new = 0

                    if iter > 1:
                        mut_two_len_new = length - calibrator
                        mut_two_calib_new = calibrator
                    else:
                        mut_two_len_new = length
                        mut_two_calib_new = 0
                else:
                    if iter == 1:
                        mut_one_len_new = length
                        mut_one_calib_new = 0
                        mut_two_len_new = length
                        mut_two_calib_new = 0
                    # mut_one_len_new = length - calibrator
                    # mut_one_calib_new = 0
                    #
                    # if iter > 2:
                    #     mut_two_len_new = length - 2
                    #     mut_two_calib_new = 2
                    # elif iter == 2:
                    #     mut_two_len_new = length - 1
                    #     mut_two_calib_new = 1
                    # else:
                    #     mut_two_len_new = length
                    #     mut_two_calib_new = 0
                    else:
                        raise Exception("this should not happen")

            else:
                mut_one_len_new = 0
                mut_one_calib_new = None
                mut_two_len_new = 0
                mut_two_calib_new = None
            # TODO check it

        elif mut_count == 1:
            if remainder == 2:
                for i in range(len(operands_old) - normal_length - mut_one_len):
                    operands[element_two_index + i] = operands_old[normal_length + i]

                verilog_code.append(
                    f"""\
        operands[{len(operands_old) - normal_length - mut_one_len + element_two_index - 1}:{element_two_index}] = operands_old[{len(operands_old) - mut_one_len - 1}:{normal_length}];\
                    """)

                mut_one_len_new = len(operands_old) - normal_length - mut_one_len
                mut_one_calib_new = calibrator if not ind else 0
                for i in range(mut_one_len):
                    operands[element_two_index + mut_one_len_new + i] = operands_old[-mut_one_len + i]

                    verilog_code.append(
                        f"""\
        operands[{element_two_index + mut_one_len_new + mut_one_len - 1}:{element_two_index + mut_one_len_new}] = operands_old[{len(operands_old) - 1}:{len(operands_old) - mut_one_len}];\
                        """)

                mut_two_len_new = mut_one_len
                mut_two_calib_new = mut_one_calib
            elif remainder == 1:
                for i in range(mut_one_len):
                    operands[element_two_index + i] = operands_old[normal_length + i]

                verilog_code.append(
                    f"""\
        operands[{element_two_index + mut_one_len - 1}:{element_two_index}] = operands_old[{normal_length + mut_one_len - 1}:{normal_length}];\
                    """)

                mut_one_len_new = mut_one_len
                mut_one_calib_new = mut_one_calib
                mut_two_len_new = 0
                mut_two_calib_new = None
        elif mut_count == 2:
            for i in range(mut_one_len):
                operands[element_two_index + i] = operands_old[normal_length + i]
            for i in range(mut_two_len):
                operands[element_two_index + mut_one_len + i] = operands_old[normal_length + mut_one_len + i]

            verilog_code.append(
                f"""\
        operands[{element_two_index + mut_one_len + mut_two_len - 1}:{element_two_index}] = operands_old[{normal_length + mut_one_len + mut_two_len - 1}:{normal_length}];\
                """)

            mut_one_len_new = mut_one_len
            mut_one_calib_new = mut_one_calib
            mut_two_len_new = mut_two_len
            mut_two_calib_new = mut_two_calib

        operands_len = element_two_index + mut_one_len_new + mut_two_len_new
        #TODO this print is fine
        #print("Length is: ", operands_len)

    mutated_operands = [] * operands_len
    mutated_operands[:] = operands[:operands_len]
    mut_count_new = 2 if (mut_one_len_new and mut_two_len_new) \
        else 1 if (mut_one_len_new or mut_two_len_new) \
        else 0
    n_new = (g + 1) * 2 + mut_count_new
    # if n == 4 and mut_count == 2:  # fix the second to last iteration
    #     n_new = 3
    n_new += health_plus_one
    if mut_count_new == 2:
        for _ in range(length + 1 - saturation - mut_one_len_new - mut_one_calib_new):
            mutated_operands.insert(-mut_two_len_new, None)
        for _ in range(length + 1 -saturation - mut_two_len_new - mut_two_calib_new):
            mutated_operands.append(None)
    elif mut_count_new == 1:
        for _ in range(length + 1 -saturation - mut_one_len_new):
            mutated_operands.append(None)

    # split_array = get_split(iter, n_new, length, mut_count_new, mut_one_calib_new, mut_two_calib_new)
    # print("Split ARRAY", iter, ": ", split_array)
    # oper1 = np.split(np.array(mutated_operands), split_array)
    # oper1 = [oper[::-1] for oper in oper1]
    # oper1 = oper1[:-1]
    # a = DataFrame(oper1, columns=(length + 1) * [''], index=n_new * [''])
    #
    # a = a.fillna('')
    # print(a)
    split_display_ind = True
    # TODO this print is fine
    #print_state(mutated_operands, n_new, length, mut_count_new, mut_one_calib_new, mut_two_calib_new, iter)

    return (operands, n_new, operands_len, mut_one_len_new, mut_one_calib_new, mut_two_len_new, mut_two_calib_new)


def wallace_simple_add(operands, max_length, mut_one_len, mut_one_calib, mut_two_len, mut_two_calib):
    verilog_code.append("""\t\tcarry = 1'b0;""");

    if (mut_one_calib + mut_one_len) >= (mut_two_calib + mut_two_len):
        length = mut_one_calib + mut_one_len
        verilog_code.append(
            f"""
        operands = {{{{ {length - mut_two_calib + mut_two_len}{{operands[{len(operands) - 1}]}} }}, {{operands[{len(operands)-1}:0]}}}};\
            """)
        for _ in range(length - mut_two_calib + mut_two_len):
            operands.append(operands[-1])
    else:
        length = mut_two_calib + mut_two_len
        verilog_code.append(
            f"""
        operands = {{{{operands[{len(operands)-1}:{mut_one_len}]}}, {{{length - mut_one_calib - mut_one_len}{{operands[{mut_one_len - 1}]}}}}, {{operands[{mut_one_len - 1}:0]}}}};\
            """)
        for _ in range(length - mut_one_calib - mut_one_len):
            operands.insert(mut_one_len, operands[mut_one_len - 1])

    result = [0] * max_length
    carry = 0
    for i in range(max_length):
        if i >= mut_one_calib and i >= mut_two_calib:
            carry, result[i] = bin_format(operands[i - mut_one_calib] + operands[length + i - mut_two_calib] + carry)
            verilog_code.append(
                f"""\
        {{carry, result[{i}]}} = operands[{i - mut_one_calib}] + operands[{length + i - mut_two_calib}] + carry;\
                """)
        elif i >= mut_one_calib:
            verilog_code.append(
                f"""\
        {{carry, result[{i}]}} = operands[{i - mut_one_calib}] + carry;\
                """)
            carry, result[i] = bin_format(operands[i - mut_one_calib] + carry)
        elif i >= mut_two_calib:
            verilog_code.append(
                f"""\
        {{carry, result[{i}]}} = operands[{mut_one_len + i - mut_two_calib}] + carry;\
                """)
            carry, result[i] = bin_format(operands[mut_one_len + i - mut_two_calib] + carry)
        else:
            verilog_code.append(
                f"""\
        result[{i}] = carry;\
                """)
            result[i] = carry
    result_dec = 0

    print(result)
    if result[max_length - 1] == 0:
        for i in range(max_length):
            result_dec += result[i] * 2 ** i
    else:
        for i in range(max_length - 1):
            result_dec += result[i] * 2 ** i
        result_dec = result_dec - 2 ** (max_length - 1)

    return result_dec
