import numpy as np
import feasst as fst
import pyfeasst

table = fst.Table3D().deserialize(pyfeasst.read_checkpoint("table"))
data = pyfeasst.vector3d_to_list(table.data())
zslice = list()
for bn in range(table.num2()):
    zslice.append(table.bin_to_value(2, bn)*10)

import matplotlib.pyplot as plt
plt.scatter(zslice, data[20][0])
plt.show()
