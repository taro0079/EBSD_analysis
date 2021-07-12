from os import path, supports_follow_symlinks
import pandas as pd
import numpy as np


class Txt_data:
    def __init__(self, path) -> None:
        self.path: str = path

    def _all_lines_read(self):
        f = open(self.path)
        return f.readlines()

    def _delete_header(self, lines):
        return [i for i in lines if i[0] != "#"]

    def _get_header_section(self, lines):
        return [line for line in lines if line[0] == "#"][2::]

    def _get_header(self, lines):
        return [line[12:-1] for line in lines]

    def _delete_blank(self, lines):
        return [line.replace(" ", "") for line in lines]

    def _line_split(self, lines):
        return [i.split() for i in lines]

    def _convert_to_float(self, lines):
        return [map(float, i) for i in lines]

    def read_data(self):
        lines = self._all_lines_read()
        header_deleted_lines = self._delete_header(lines)
        header_range = self._get_header_section(lines)
        header = self._get_header(header_range)
        header_blank_deleted = self._delete_blank(header)
        splited_lines = self._line_split(header_deleted_lines)
        float_list = self._convert_to_float(splited_lines)
        df = pd.DataFrame(float_list)
        df.columns = header_blank_deleted
        return Data(df)


class Data:
    def __init__(self, data) -> None:
        self.data = data


def main():
    textdata = Txt_data(path="./Nb2Ti1Hf-Cu14Sn_grain.txt")
    data = textdata.read_data()
    print(data.data)


if __name__ == "__main__":
    main()
