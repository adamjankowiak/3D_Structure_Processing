import Bio
from Bio import PDB
from Bio.PDB.PDBList import PDBList

pdbl = PDBList()
# fetch_pdb = pdbl.retrieve_pdb_file('4ywo', file_format='pdb')
pdb_parser = PDB.PDBParser()
structure = pdb_parser.get_structure("4ywo", "yw/pdb4ywo.ent")
def calculate_phi_psi(structure):
    phi_psi_angles = []
    for model in structure:
        for chain in model:
            polypeptides = PDB.PPBuilder().build_peptides(chain)
            phi_psi_angles.extend(polypeptides[0].get_phi_psi_list())
    return phi_psi_angles

angles = calculate_phi_psi(structure)
print(angles)
print("hi")