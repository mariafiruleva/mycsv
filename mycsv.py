import os
import re


def read_csv(path, sep=","):
    """
    This function takes two arguments:
    1st is the path to the file,
    2nd is the delimiter character.
    :param path: path to input file
    :param sep: separator which will be used for string separation (default: comma)
    :return: list of separated strings
    """
    result = list()
    if not os.path.isfile(path):
        print("Error, such file doesn't exist")
        return result
    with open(path, 'r') as in_file:
        for line in in_file:
            result.append(re.findall(fr'".+?"|[^"{sep} ]+', line.strip()))
    return result


def write_csv(path_to_csv_file, data, sep=","):
    """
    This function should save data from data variable to the output file.
    :param path_to_csv_file: path to output file
    :param data: list with input data
    :param sep: separator which will be used for string separation in output file (default: comma)
    :return: None
    """
    with open(path_to_csv_file, 'w') as out_file:
        for line in data:
            string = f"{sep}".join(map(lambda x: f'"{x}"' if ',' in x else x, line))
            out_file.write(f'{string}\n')
