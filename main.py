from os import path, supports_follow_symlinks
from typing import List
import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.stats.stats import sem


class Txt_data:
    def __init__(self, path) -> None:
        self.path: str = path

    # text fileの全ての行を読み込む
    def _all_lines_read(self) -> List[str]:
        f = open(self.path)
        return f.readlines()

    # "#"部分をheaderとしてその部分を削除したリスト配列を出力
    def _delete_header(self, lines) -> List[str]:
        return [i for i in lines if i[0] != "#"]

    # header領域を抽出
    def _get_header_section(self, lines, row_start) -> List[str]:
        return [line for line in lines if line[0] == "#"][row_start::]

    # header領域からheader部分のみを抽出
    def _get_header(self, lines) -> List[str]:
        return [line[12:-1] for line in lines]

    # 文字列中の空白を削除
    def _delete_blank(self, lines) -> List[str]:
        return [line.replace(" ", "") for line in lines]

    # 文字列を分割
    def _line_split(self, lines) -> List[str]:
        return [i.split() for i in lines]

    # floatにマップする
    def _convert_to_float(self, lines) -> List[float]:
        return [map(float, i) for i in lines]

    # データ読み込みの一連の処理 結晶粒径などの微細組織ファイル用
    def read_morphological_data(self):
        lines = self._all_lines_read()
        header_deleted_lines = self._delete_header(lines)
        header_range = self._get_header_section(lines, row_start=2)
        header = self._get_header(header_range)
        header_blank_deleted = self._delete_blank(header)
        splited_lines = self._line_split(header_deleted_lines)
        float_list = self._convert_to_float(splited_lines)
        df = pd.DataFrame(float_list)
        df.columns = header_blank_deleted
        return morphological_data(df)


class morphological_data:
    def __init__(self, data) -> None:
        self.data = data
        self.grainarea = data.iloc[:, 1]
        # grain areaからそれぞれのgrain sizeを算出
        self.grainsize = 2 * np.sqrt(data.iloc[:, 1] / np.pi)
        self.aspectratio = data.iloc[:, 3]
        self.theshold = 0.5  # theshold for columar of equiaxed grains

    # Number fractionの算出関数（単なる平均）
    def _calculation_of_number_fraction(self, inputdata):
        return np.average(inputdata)

    # Area fractionの算出関数（grain areaで重み付けした加重平均）
    def _calclation_of_area_fraction(self, inputdata):
        return np.average(inputdata, weights=self.grainarea)

    # 結晶粒径(um)の平均値を出し1000倍してnmに変換 (in number fraction)
    def calc_grainsize_numberfraction(self):
        return self._calculation_of_number_fraction(self.grainsize) * 1000

    # 結晶粒径(um)の平均値を出し1000倍してnmに変換 (in area fraction)
    def calc_grainsize_areafraction(self):
        return self._calclation_of_area_fraction(self.grainsize) * 1000

    # 平均アスペクト比を算出し100をかけて%に変換 (in number fraction)
    def calc_aspectratio_numberfraction(self):
        return self._calculation_of_number_fraction(self.aspectratio) * 100

    # 平均アスペクト比を算出し100をかけて%に変換 (in area fraction)
    def calc_aspectratio_areafraction(self):
        return self._calclation_of_area_fraction(self.aspectratio) * 100

    #  結晶面積の合計（= 分析面積）
    def _get_grainarea_sum(self):
        return np.sum(self.grainarea)

    # 等軸結晶の基準を設定
    def _get_equiaxed_condition(self) -> List[bool]:
        return (self.aspectratio >= self.theshold).values

    # 等軸結晶領域の面積合計
    def _sum_of_equiaxed(self) -> float:
        return np.sum(self.grainarea[self._get_equiaxed_condition()])

    # 等軸結晶領域比（%）の算出関数
    def get_equiaxedgrain_areafraction(self) -> float:
        return self._sum_of_equiaxed() / self._get_grainarea_sum() * 100

    # 柱状結晶の基準を設定
    def _get_columnar_condition(self) -> List[bool]:
        return (self.aspectratio < self.theshold).values

    # 柱状結晶領域の面積合計を算出
    def _sum_of_columnar(self) -> float:
        return np.sum(self.grainarea[self._get_columnar_condition()])

    # 等軸結晶領域比（%）を算出
    def _get_columnargrain_areafraction(self) -> float:
        return self._sum_of_columnar() / self._get_grainarea_sum() * 100

    def create_data(self):
        df = pd.DataFrame({
            'avegsnum': [self.calc_grainsize_numberfraction()],
            'avegsarea': [self.calc_grainsize_areafraction()],
            'gsstd': [np.std(self.grainsize)*1000],
            'gssem': [stats.sem(self.grainsize)*1000],
            'garnum': [self.calc_aspectratio_numberfraction()],
            'gararea': [self.calc_aspectratio_areafraction()],
            'garstd': [np.std(self.aspectratio)*100],
            'garsem': [stats.sem(self.aspectratio)*100],
            'equiaxedArea': [self.get_equiaxedgrain_areafraction()],
            'columnarArea': [self._get_columnargrain_areafraction()]
        })
        return df


def main():
    textdata = Txt_data(path="./Nb2Ti1Hf-Cu14Sn_grain.txt")
    data = textdata.read_morphological_data()
    gs = data.calc_grainsize_areafraction()
    gar = data.calc_aspectratio_areafraction()
    eq = data.get_equiaxedgrain_areafraction()
    print(data.create_data())


if __name__ == "__main__":
    main()
