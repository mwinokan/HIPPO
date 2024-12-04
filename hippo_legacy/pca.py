# https://github.com/rdkit/rdkit-tutorials/blob/master/notebooks/005_Chemical_space_analysis_and_visualization.ipynb

import numpy as np
from rdkit import DataStructs
from rdkit.Chem import AllChem, rdFingerprintGenerator


class FP:
    """
    Molecular fingerprint class, useful to pack features in pandas df

    Parameters
    ----------
    fp : np.array
        Features stored in numpy array
    names : list, np.array
        Names of the features
    """

    def __init__(self, fp, names):
        self.fp = fp
        self.names = names

    def __str__(self):
        return "%d bit FP" % len(self.fp)

    def __len__(self):
        return len(self.fp)


def get_cfps(
    mol, radius=2, nBits=1024, useFeatures=False, counts=False, dtype=np.float32
):
    """Calculates circular (Morgan) fingerprint.
    http://rdkit.org/docs/GettingStartedInPython.html#morgan-fingerprints-circular-fingerprints

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
    radius : float
        Fingerprint radius, default 2
    nBits : int
        Length of hashed fingerprint (without descriptors), default 1024
    useFeatures : bool
        To get feature fingerprints (FCFP) instead of normal ones (ECFP), defaults to False
    counts : bool
        If set to true it returns for each bit number of appearances of each substructure (counts). Defaults to false (fingerprint is binary)
    dtype : np.dtype
        Numpy data type for the array. Defaults to np.float32 because it is the default dtype for scikit-learn

    Returns
    -------
    ML.FP
        Fingerprint (feature) object
    """
    arr = np.zeros((1,), dtype)

    if counts is True:
        info = {}
        fp = AllChem.GetHashedMorganFingerprint(
            mol, radius, nBits, useFeatures=useFeatures
        )
        DataStructs.ConvertToNumpyArray(fp, arr)
    else:

        # https://greglandrum.github.io/rdkit-blog/posts/2023-01-18-fingerprint-generator-tutorial.html#additional-information-explaining-bits
        fmgen = rdFingerprintGenerator.GetMorganGenerator(
            radius=radius,
            fpSize=nBits,
            atomInvariantsGenerator=rdFingerprintGenerator.GetMorganFeatureAtomInvGen(),
        )

        assert not useFeatures

        DataStructs.ConvertToNumpyArray(
            fmgen.GetFingerprint(mol),
            # AllChem.GetMorganFingerprintAsBitVect(
            #     mol, radius, nBits=nBits, useFeatures=useFeatures
            # ),
            arr,
        )
    return FP(arr, range(nBits))
