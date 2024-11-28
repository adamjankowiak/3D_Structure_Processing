from Bio.PDB.PDBList import PDBList
from Bio import PDB
from Bio.PDB import PDBIO,Select
import argparse

class CoarseGrainSelect(Select):
    def accept_atom(self, atom):
        residue = atom.get_parent()
        atom_name = atom.name.strip()
        residue_name = residue.get_resname().strip()

        # For purines (A, G): N9, C2, C6
        if residue_name in ["A", "G"]:
            if atom_name in ["N9", "C2", "C6"]:
                return True

        # For pyrimidines (C, U): N1, C2, C4
        if residue_name in ["C", "U"]:
            if atom_name in ["N1", "C2", "C4"]:
                return True

        # For backbone: P, C4'
        if residue_name in ["A", "C", "G" , "U"]:
            if atom_name in ["P", "C4'"]:
                return True

        return False

def convert_to_grain(output, structure):
    io = PDBIO()
    io.set_structure(structure)
    io.save(output,select=CoarseGrainSelect())
    pass


def main(input_file,output_file):
    mmcif_parser = PDB.MMCIFParser()
    pdb_parser = PDB.PDBParser(QUIET=True)
    structure = pdb_parser.get_structure(input_file, input_file)
    convert_to_grain(args.output, structure)



    # structure = structure_load(input_file)

    print(f"Plik zostal zapisany jako {output_file} w katalogu Struktury")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Przemiana z struktury pelnej na gruboziarnista")
    parser.add_argument('--input', type=str, required=True, help='Plik Wejsciowy.')
    parser.add_argument('--output', type=str, required=True, help='plik wyjsciowy.')
    args = parser.parse_args()
    main(args.input,args.output)