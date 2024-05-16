import numpy as np

def getFixedPoints(lam):
    if lam**2 < 6:
        print("B:\tx = " + str(x := lam/np.sqrt(6)) +   "\ty = " + str(y := np.sqrt(1 - (lam**2)/6)) + "\tz = " + str(0) + "\tgamma = " + str((2*x**2)/(x**2 + y**2)))

    if lam**2 > 3:
        print("C:\tx = " + str(x := np.sqrt(3/2)/lam) + "\ty = " + str(y := np.sqrt(3/2)/lam)        + "\tz = " + str(0) + "\tgamma = " + str((2*x**2)/(x**2 + y**2)))

    if lam**2 > 4:
        print("E:\tx = " + str(x := np.sqrt(8/3)/lam) + "\ty = " + str(y := 2/(lam * np.sqrt(3)))   + "\tz = " + str(np.sqrt(1 - 4/(lam**2))) + "\tgamma = " + str((2*x**2)/(x**2 + y**2)))
