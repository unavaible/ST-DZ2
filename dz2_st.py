import numpy as np
from sympy.utilities.iterables import multiset_permutations
from prettytable import PrettyTable

n = 7
k = 4
gen_polynom = 0b1011

information = 0b0001


# util to transform binary number to array of powers of 2
# Утилита для преобразования двоичного числа в массив степеней двойки
def transform_to_arr(bin_num):
    i = 0
    arr = np.delete(np.array([0]), 0)
    while (bin_num > 0):
        if bin_num % 2:
            arr = np.append(arr, i)
        i += 1
        bin_num //= 2
    return arr


# Утилита для преобразования массива нулей и единиц в число
def arr_to_num(arr):
    num = 0
    for i in range(len(arr)):
        num += arr[i] << i
    return num


# util to get syndrome
# Утилита для получения синдрома
def get_syndrome(shifted_information):
    if shifted_information == 0:
        return 0
    info_arr = transform_to_arr(shifted_information)
    polynom_arr = transform_to_arr(gen_polynom)
    polynom_arr_max = max(polynom_arr)
    while polynom_arr_max <= max(info_arr):
        quotient = max(info_arr) - max(polynom_arr)
        current = polynom_arr + quotient
        i = 0

        # xor realization
        # Реализация сложения по модулю 2
        while i < len(info_arr):
            element = info_arr[i]
            if (element in current):
                current = np.delete(current, np.where(current == element))
                info_arr = np.delete(info_arr, i)
                i -= 1
            i += 1
        i = 0
        for i in range(len(current)):
            if (current[i] not in info_arr):
                info_arr = np.append(info_arr, current[i])

        # if get_syndrome passed without remainder
        if (len(info_arr) == 0 and len(current) == 0):
            break
    if (len(info_arr) == 0):
        return 0
    syndrome = 0
    for i in info_arr:
        syndrome += 1 << i
    return syndrome


# Основная функция кодирования сигнала
def coding():
    shifted_information = information << n - k
    syndrome = get_syndrome(shifted_information)
    shifted_information += syndrome
    return shifted_information


# Утилита для получения всех ошибок в классе
def get_all_errors(mult):
    temp = [0] * 7
    for i in range(mult):
        temp[i] = 1
    all_errors = list(multiset_permutations(temp))
    return all_errors


# Утилита для получения вектора ошибки по синдрому
def get_error_vector(syndrome):
    error_vector = 0
    all_errors = get_all_errors(1)
    for vector in all_errors:
        vector_num = arr_to_num(vector)
        if (get_syndrome(vector_num) == syndrome):
            error_vector = vector_num
    return error_vector

# Определяем исправлена и/или обнаружена ошибка
def is_error_corrected(error_vector, all_errors):
    encoded_info = coding()
    error_info = encoded_info ^ error_vector
    syndrome = get_syndrome(error_info)
    supposed_error_vector = get_error_vector(syndrome)
    corrected_info = error_info ^ supposed_error_vector
    initial_info = corrected_info >> n - k
    if (initial_info == information):
        return 2
    elif syndrome != 0:
        return 1
    else:
        return 0

# Вывод результатов вычисления корректирующей способности в виде таблицы
def Ck_calc():
    table = PrettyTable()
    Nk_vector = [0] * 7
    No_vector = [0] * 7
    Co_vector = [0] * 7
    Ck_vector = [0] * 7
    i_vector = [0] * 7
    Cn_vector = [0] * 7
    for i in range(1, 8):
        i_vector[i - 1] = i
        Nk = 0
        No = 0
        all_errors = get_all_errors(i)
        Cn_vector[i - 1] = len(all_errors)
        for error_vector in all_errors:
            flag = is_error_corrected(arr_to_num(error_vector), all_errors)
            if flag == 2:
                Nk += 1
                No += 1
            elif flag == 1:
                No += 1
            else:
                print("error not found")
        Nk_vector[i - 1] = Nk
        No_vector[i - 1] = No
        Ck_vector[i - 1] = Nk / Cn_vector[i - 1] * 100
        Co_vector[i - 1] = No / Cn_vector[i - 1] * 100
    table.add_column("i", i_vector)
    table.add_column("Nk", Nk_vector)
    table.add_column("No", No_vector)
    table.add_column("Cni", Cn_vector)
    table.add_column("Ck", Ck_vector)
    table.add_column("Co", Co_vector)
    return table


def main():
    encoded_info = coding()
    print("------------------------------TEST------------------------------")
    print(f"encoded information: {bin(encoded_info)}")
    e_test = 0b100
    info_error = encoded_info ^ e_test
    print(f"encoded information with test error vector (0100): {bin(info_error)}")
    syndrome_test = get_syndrome(info_error)
    print(f"syndrome: {bin(syndrome_test)}")
    error_vector = get_error_vector(syndrome_test)
    print(f"computed error vector: {bin(error_vector)}")
    info_corrected = info_error ^ error_vector
    print(f"corrected information: {bin(info_corrected)}")
    print(f"initial information: {bin(info_corrected >> n - k)}")
    print("------------------------------TEST------------------------------\n\n")
    print("------------------------------TABLE---------------------------")
    print(Ck_calc())


if __name__ == "__main__":
    main()
