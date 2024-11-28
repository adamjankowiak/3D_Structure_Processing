import ramachandraw.utils
import ramachandraw.parser

pdb_path = ramachandraw.utils.fetch_pdb("4YWO")
ramachandraw.utils.plot(pdb_path, show=1)