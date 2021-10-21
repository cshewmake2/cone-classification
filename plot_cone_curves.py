import numpy as np
import matplotlib.pyplot as plt
import cone

lmbda = np.arange(400,710,10)
cL = list(map(cone.cone_L,lmbda))
cM = list(map(cone.cone_M,lmbda))
cS = list(map(cone.cone_S,lmbda))

plt.figure()
plt.plot(lmbda, cL)
plt.plot(lmbda, cM)
plt.plot(lmbda, cS)
plt.show()
