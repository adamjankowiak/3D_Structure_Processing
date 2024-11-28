from Bio.PDB.PDBList import PDBList
from Bio import PDB
import numpy as np
import matplotlib.pyplot as plt

pdbl = PDBList()
fetch_pdb = pdbl.retrieve_pdb_file('4ywo', file_format='pdb')

mmcif_parser = PDB.MMCIFParser()
pdb_parser = PDB.PDBParser()

structure = pdb_parser.get_structure("4ywo", "yw/pdb4ywo.ent")

list_atoms = []
for model in structure:
  for chain in model:
      for res in chain:
          for atom in res.get_atoms():
              if atom.get_name() == 'CA':
                  list_atoms.append(atom.get_coord())
                  print(f"Nazwa reszty: {res.get_resname()} | Nazwa atomu: {atom.get_name()} | Koordynaty atomu: {atom.get_coord()}")

num_atoms = len(list_atoms)
matrix = np.zeros((num_atoms, num_atoms),dtype=int)

max_distance = 8.0

for i in range(num_atoms):
    for j in range(i+1, num_atoms):
        distance = np.linalg.norm(list_atoms[i] - list_atoms[j])
        # print(distance)
        if distance < max_distance:
            matrix[i][j] = 1
            matrix[j][i] = 1

x,y=np.where(matrix==1)

print(matrix)
plt.scatter(x,y)
plt.scatter(x,y)
plt.axline((0,0),slope=1)
plt.grid(True)
plt.show()