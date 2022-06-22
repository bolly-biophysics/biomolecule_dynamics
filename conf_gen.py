from prody import *
from matplotlib.pylab import *

def gammaDistanceDependent(dist2, *args):
    if dist <= 16:
        return 10
    elif dist2 <= 100:
        return 2
    elif dist2 <= 225:
        return 1
    else:
        return 0

mol=parsePDB('1y26.pdb')
anm=ANM('1y26_nma')
anm.buildHessian(mol, cutoff=15, gamma=gammaDistanceDependent)
anm.calcModes()

writeNMD('1y26_anm.nmd', anm, mol)

ensemble=sampleModes(anm[:3], mol, n_confs=200, rmsd=2.0)
conf_gen=mol.copy()
conf_gen.addCoordset(ensemble)
writePDB('1y26_ensemble.pdb', conf_gen)
