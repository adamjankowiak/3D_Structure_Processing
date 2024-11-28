import numpy as np
import math
import Bio
from Bio import PDB
from Bio.PDB.PDBList import PDBList

def degrees(rad_angle) :
    if rad_angle is None :
        return None
    angle = rad_angle * 180 / math.pi
    while angle > 180 :
        angle = angle - 360
    while angle < -180 :
        angle = angle + 360
    return angle

def calculate_phi_psi(structure):
    phi_psi_list = np.zeros([2])
    for model in structure:
        for chain in model:
            poly = Bio.PDB.Polypeptide.Polypeptide(chain)
            phi_psi = (degrees(poly.get_phi_psi_list()[1][0]),degrees(poly.get_phi_psi_list()[1][1]))
            v = np.array(phi_psi)
            phi_psi_list = np.vstack((phi_psi_list,v))
    return phi_psi_list[1:,:]

pdbl = PDBList()
# fetch_pdb = pdbl.retrieve_pdb_file('4ywo', file_format='pdb')
pdb_parser = PDB.PDBParser()
structure = pdb_parser.get_structure("4ywo", "yw/pdb4ywo.ent")
angles = calculate_phi_psi(structure)
print(angles)
print("hi")