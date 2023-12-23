import numpy as np
import matplotlib.pyplot as plt
import json

file_path = "ResultData/xJBd_wave_N200.json"
with open(file_path, 'r') as fi:
    data = json.load(fi)

def xJ_plot(data):
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

def xT_plot(data):
    N = int(data["N"])
    J = float(data["J"])
    wave = data["wave"]
    atol = float(data["atol"])
    rtol = float(data["rtol"])
    xList = np.array(data["xList"])
    TList = np.array(data["TList"])
    DeltaList = np.array(data["DeltaList"])
    # print(xList[40], JList[45], DeltaList[40,45])

    #plot
    plt.imshow(DeltaList.transpose(),extent=[xList.min(), xList.max(), TList.min(), TList.max()],
            cmap="hot", origin="lower", aspect='auto')
    plt.xlabel("x(dopping)")
    plt.ylabel("T/t")
    plt.colorbar(label=r"$\Delta$")
    plt.title("d-wave")
    plt.show()

def JT_plot(data):
    N = int(data["N"])
    x = float(data["x"])
    wave = data["wave"]
    atol = float(data["atol"])
    rtol = float(data["rtol"])
    JList = np.array(data["JList"])
    TList = np.array(data["TList"])
    DeltaList = np.array(data["DeltaList"])
    # print(xList[40], JList[45], DeltaList[40,45])

    #plot
    plt.imshow(DeltaList.transpose(),extent=[JList.min(), JList.max(), TList.min(), TList.max()],
            cmap="hot", origin="lower", aspect='auto')
    plt.xlabel("J/t")
    plt.ylabel("T/t")
    plt.colorbar(label=r"$\Delta$")
    plt.title("d-wave")
    plt.show()

# JT_plot(data)
xList = np.array(data["xList"])
JList = np.array(data["JList"])
DeltaList = np.array(data["DeltaList"], dtype=float)
DeltaSCList = np.array(data["DeltaSCList"])
BList = np.array(data["BList"], dtype=float)
fig = plt.figure()
ax = fig.gca()
ax.plot(xList, DeltaList[:,40])
ax.plot(xList, BList[:,40]*0.40*2)
ax.plot(xList, DeltaSCList[:,40])
# ax.set_xbound(0, 0.4)
print(JList[40])
plt.show()