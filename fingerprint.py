import sys
from rdkit import Chem
from rdkit.Chem import AllChem

sdf = Chem.SDMolSupplier( sys.argv[1] )
mols = []
for mol in sdf:
    try:
        mols.append( mol )
    except:
        pass

fps = {}
for mol in mols:
    try:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol,2, nBits=1024)
        fp_bit = fp.ToBitString()
        string = ""
        for bit in fp_bit:
            string += bit+"\t"
        fps[ mol.GetProp("COMPOUND_ID") ] = string.rstrip("\t")
    except:
        pass
print len(fps)

f = open("bitmatrix.txt", "w")
for k, v in fps.items():
    f.write( "%s\t%s" % (k, v) )
    f.write( "\n" )

f.close()
