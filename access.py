import numpy as np

data = []

for l in open("dccm-p.dat"):
    row = [x for x in l.split()]
    data.append(row)
        
cross_correlations = np.zeros((71,71))
cross_correlations[1:,1:] = data

numpy_list = np.asarray(cross_correlations)
print numpy_list.shape
np.savetxt("./dccm-p1.dat", numpy_list)
