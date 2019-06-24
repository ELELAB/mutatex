from Bio import PDB
from six import iteritems

ptm_residues = {"y": "PTR",
                "p": "TPO",
                "s": "SEP",
                "h": "HYP",
                "z": "TYS",
                "k": "MLZ",
                "m": "MLY",
                "l": "M3L",
                "o": "H1S",
                "e": "H2S",
                "f": "H3S"}
len_three2one = len(PDB.Polypeptide.d1_to_index)
idx=len_three2one
for k,v in iteritems(ptm_residues):
    PDB.Polypeptide.d1_to_index[k] = idx
    PDB.Polypeptide.dindex_to_1[idx] = k

    PDB.Polypeptide.d3_to_index[v] = idx
    PDB.Polypeptide.dindex_to_3[idx] = v
    idx += 1

from . import utils
