import numpy as np
from scipy.interpolate import interp1d


class SpectrumData:
    def __init__(self, file_name=None):
        self.data: np.ndarray = None
        self.range = None
        self.step = None
        self.num_points = None
        self.normalised = False
        self.load_from_file(file_name)

    def load_from_file(self, file_name):
        self.data = np.loadtxt(file_name, delimiter=",")
        self.range = (self.data[0, 0], self.data[-1, 0])
        self.step = self.data[1, 0] - self.data[0, 0]
        self.num_points = self.data.shape[0]
        self.data = self.data[:, 1:]
        
    def normalise(self):
        self.data = self.data / np.max(self.data)
        self.normalised = True
        return self

    def reshape(self, target_range):
        original_range = np.linspace(self.range[0], self.range[1], self.num_points)
        new_range = np.linspace(target_range[0], target_range[1], target_range[2])
        expanded_data = np.zeros((target_range[2], self.data.shape[1]))
        for i in range(self.data.shape[1]):
            f = interp1d(original_range, self.data[:, i], kind="linear", fill_value='extrapolate')
            expanded_data[:, i] = f(new_range)
        self.data = expanded_data
        self.range = target_range[:2]
        self.num_points = target_range[2]
        self.step = (self.range[1] - self.range[0]) / (self.num_points - 1)
        return self


def main():
    file_name = (
        "/Users/jackchou/Desktop/Code/colour-notebook/data/CIE_std_illum_D65.csv"
    )
    std_illumi_D65 = SpectrumData(file_name).reshape((380, 780, 401)).normalise()
    print(std_illumi_D65.range)
    print(std_illumi_D65.step)
    print(std_illumi_D65.data.shape)


if __name__ == "__main__":
    main()
