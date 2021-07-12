import os
from typing import List
import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.stats.stats import sem


class Txt_data:
    def __init__(self, path) -> None:
        self.path: str = path

    def get_basename(self):
        basename = os.path.basename(self.path)
        return basename.split(sep=".")[0]

    def _all_lines_read(self):
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

    # misorientation用読み込み関数
    def read_misorientation_data(self):
        lines = self._all_lines_read()
        header_deleted_lines = self._delete_header(lines)
        splited_lines = self._line_split(header_deleted_lines)
        float_list = self._convert_to_float(splited_lines)
        df = pd.DataFrame(float_list)[:, :5]
        return misorientation(df)


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
        return csvRepo(df)


class misorientation:
    def __init__(self, data) -> None:
        self._phi1 = data.iloc[:, 0]
        self._PHI = data.iloc[:, 1]
        self._phi2 = data.iloc[:, 2]
        self._xdata = data.iloc[:, 3]
        self._ydata = data.iloc[:, 4]
        self._vc = np.array([0,0,1])

    # オイラー角 -> デカルト座標系変換用行列1
    def _make_g1(self):
        return np.array([[np.cos(self._phi1), np.sin(self._phi1), 0], [-np.sin(self._phi1), np.cos(self._phi1), 0], [0, 0, 1]])

    # オイラー角 -> デカルト座標系変換用行列2
    def _make_g2(self):
        return np.array([[1, 0, 0], [0, np.cos(self._PHI), np.sin(self._PHI)], [0, -np.sin(self._PHI), np.cos(self._PHI)]])

    # オイラー角 -> デカルト座標系変換用行列3
    def _make_g3(self):
        return np.array([np.cos(self._phi2), np.sin(self._phi2), 0], [-np.sin(self._phi2), np.cos(self._phi2), 0], [0, 0, 1])
    
    # オイラー角 -> デカルト座標系変換用行列
    def make_gmat(self):
        return np.dot(np.dot(self._make_g1(), self._make_g2()), self._make_g3())


class csvRepo:
    def __init__(self, data) -> None:
        self.data = data

    def _make_path(self, writepath):
        return writepath + '.csv'

    def write(self, writepath):
        path = self._make_path(writepath)
        self.data.to_csv(path)


def main():
    textdata = Txt_data(path="./Nb4Ta1Hf-Cu14Sn_grain.txt")
    data = textdata.read_morphological_data()
    writedata = data.create_data()
    writedata.write(textdata.get_basename())


if __name__ == "__main__":
    main()
