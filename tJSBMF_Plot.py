import numpy as np
import matplotlib.pyplot as plt
import json

file_path = "ResultData/s_wave.json"
with open(file_path, 'r') as fi:
    data = json.load(fi)

N = int(data["N"])
T = float(data["T"])
wave = data["wave"]
atol = float(data["atol"])
rtol = float(data["rtol"])
xList = np.array(data["xList"])
JList = np.array(data["JList"])
DeltaList = np.array(data["DeltaList"], dtype=float)
# print(xList[40], JList[45], DeltaList[40,45])

#plot
plt.imshow(DeltaList.transpose(),extent=[xList.min(), xList.max(), JList.min(), JList.max()],
        cmap="hot", origin="lower")
plt.xlabel("x(dopping)")
plt.ylabel("J/t")
plt.colorbar(label=r"$\Delta$")
plt.title("d-wave")
plt.show()