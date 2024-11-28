from Bio.PDB.PDBList import PDBList
from Bio import PDB
from Bio.PDB import PDBIO, Select, Superimposer
from Bio.PDB.Structure import Structure
import argparse

def convert_to_full(input_file, output_file):
    mmcif_parser = PDB.MMCIFParser()
    pdb_parser = PDB.PDBParser(QUIET=True)
    structure = pdb_parser.get_structure(input_file, input_file)
    template_adenine = pdb_parser.get_structure("adenine", 'templates/adenine.pdb')
    template_cytosine = pdb_parser.get_structure("cytosine", 'templates/cytosine.pdb')
    template_guanine = pdb_parser.get_structure("guanine", 'templates/guanine.pdb')
    template_uracil = pdb_parser.get_structure("uracil", 'templates/uracil.pdb')

    acgt_dict = {"A": template_adenine[0]["A"][26],
                 "C":template_cytosine[0]["A"][29],
                 "G":template_guanine[0]["A"][24],
                 "U":template_uracil[0]["A"][11]
                 }

    new_structure = Structure(output_file)
    for model in structure:
        new_model = PDB.Model.Model(model.id)
        new_structure.add(new_model)
        for chain in model:
            new_chain = PDB.Chain.Chain(chain.id)
            new_model.add(new_chain)
            for residue in chain:
                full_atom_res = acgt_dict[residue.get_resname().strip()].copy()
                full_atom_res.id = (" ",residue.id[1]," ")
                atom_list_grain = [atom for atom in residue]
                atom_list_template = [full_atom_res[atom.get_name().strip()] for atom in atom_list_grain]
                imposer = Superimposer()
                imposer.set_atoms(atom_list_grain,atom_list_template)
                for atom in full_atom_res:
                    atom.transform(imposer.rotran[0],imposer.rotran[1])

                new_chain.add(full_atom_res)

    io = PDB.PDBIO()
    io.set_structure(new_structure)
    io.save(output_file)

    pass


def main(input_file,output_file):
    convert_to_full(input_file,output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Przemiana z struktury pelnej na gruboziarnista")
    parser.add_argument('--input', type=str, required=True, help='Plik Wejsciowy.')
    parser.add_argument('--output', type=str, required=True, help='plik wyjsciowy.')
    args = parser.parse_args()
    main(args.input,args.output)

