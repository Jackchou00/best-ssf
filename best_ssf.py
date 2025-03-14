import numpy as np
from data_init import SpectrumData
from matplotlib import pyplot as plt

CIEXYZ2015_10 = (
    SpectrumData("data/CIE_xyz_1931_2deg.csv")
    .reshape((400, 700, 301))
    .normalise()
    .data
)
PMCC_Reflectance = SpectrumData("data/PMC_R.csv").reshape((400, 700, 301)).normalise().data
# Reflectance = Reflectance[:, 19]

def virtual_spectra(wavelengths, center=500, width=10):
    return np.exp(-0.5 * ((wavelengths - center) / width) ** 2)


wavelengths = np.linspace(400, 700, 301)
Reflectance = np.zeros((301, 31))
Reflectance[:, :30] = PMCC_Reflectance
Reflectance[:, 30] = virtual_spectra(wavelengths, 500, 10)

plt.figure()
plt.plot(wavelengths, Reflectance)
plt.xlabel("wavelength (nm)")
plt.ylabel("reflectance")
plt.show()

max_response = wavelengths / 700

lms_samples = Reflectance.T @ CIEXYZ2015_10
# lms_samples: (30, 3)


def generate_CCM(LMS, lms_samples, max_response, max_iter=10000):
    best_CCM = None
    best_rgb_sum_ratio = 0
    for i in range(max_iter):
        # random initialisation
        current_ccm = np.random.rand(3, 3)
        # judge if the matrix is invertible
        if np.linalg.matrix_rank(current_ccm) != 3:
            continue
        # evaluate ccm
        rgb_samples = lms_samples @ current_ccm
        rgb_sum = rgb_samples.sum()
        ssf = LMS @ current_ccm
        ratio = np.max(ssf, axis=1) / max_response
        max_ratio = np.max(ratio)
        rgb_sum_ratio = rgb_sum / max_ratio
        if rgb_sum_ratio > best_rgb_sum_ratio:
            best_rgb_sum_ratio = rgb_sum_ratio
            best_CCM = current_ccm / max_ratio
    return best_CCM


# 生成CCM矩阵
CCM = generate_CCM(CIEXYZ2015_10, lms_samples, max_response)

# 验证结果
if CCM is not None:  # 如果成功生成CCM矩阵
    print("找到可行的CCM矩阵：")
    print(CCM)
    print("行列式：", np.linalg.det(CCM))
    # 检查SSF约束
    SSF = CIEXYZ2015_10 @ CCM
    constraint_met = True
    for i in range(len(max_response)):
        if any(SSF[i] > max_response[i]):
            constraint_met = False
            break
    print("约束条件满足：", constraint_met)
    # 检查RGB相等
    V = CIEXYZ2015_10.sum(axis=0)
    RGB = V @ CCM
    print("RGB值：", RGB)
    print("RGB是否相等：", np.allclose(RGB[0], RGB[1]) and np.allclose(RGB[1], RGB[2]))
else:
    print("未能生成有效的CCM矩阵")


SSF = CIEXYZ2015_10 @ CCM

# CCM 是否可逆
print("CCM是否可逆：", np.linalg.matrix_rank(CCM) == 3)

plt.plot(wavelengths, SSF)
plt.plot(wavelengths, max_response, "--")
plt.xlabel("wavelegnth (nm)")
plt.ylabel("relative response")
plt.title("Spectral Sensitivity Function")
plt.show()
