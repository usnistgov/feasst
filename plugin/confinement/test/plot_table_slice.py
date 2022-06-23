import numpy as np
import feasst as fst

def read_checkpoint(filename):
    """Return contents of checkpoint file as a string"""
    with open (filename, "r") as myfile:
        checkpoint=myfile.readlines()
    assert(len(checkpoint) == 1)  # checkpoint files should have only one line
    return checkpoint[0]

def vector3d_to_list(vec):
    """ converts a swig stl vector to python list """
    lst = list()
    for _, vec1 in enumerate(vec):
        lst2 = list()
        for _, vec2 in enumerate(vec1):
            lst3 = list()
            for _, vec3 in enumerate(vec2):
                lst3.append(vec3)
            lst2.append(lst3)
        lst.append(lst2)
    return lst

table = fst.Table3D().deserialize(read_checkpoint("table"))
data = vector3d_to_list(table.data())
zslice = list()
for bn in range(table.num2()):
    zslice.append(table.bin_to_value(2, bn)*10)

import matplotlib.pyplot as plt
plt.scatter(zslice, data[20][0])
plt.show()
