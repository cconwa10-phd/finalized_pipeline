# Given a molecule and the list of atoms in the molecule,
# generate an equivalent SMILES string where the atoms are in
# the same order as the list.
# Caution: this version is limited to 100 bonds.

from rdkit import Chem
from rdkit.Chem.rdchem import BondType
import random
from rdkit.Chem import Draw

def format_atom(atom):
    # A general purpose SMILES writer would check to see if the atom
    # can be written without the []s and if so, only show the symbol.
    # This version is easier, though incomplete (no isotope support)
    # The goal is to get the point across and not make the best SMILES
    smiles = '['
    if atom.GetIsAromatic():
        smiles += atom.GetSymbol().lower()
    else:
        smiles += atom.GetSymbol()

    hcount = atom.GetNumExplicitHs() + atom.GetNumImplicitHs()
    if hcount != 0:
        smiles += "H%s" % hcount

    charge = atom.GetFormalCharge()
    if charge != 0:
        # Represent "[Cl-]" and "[NH4+]" as "[Cl-1]" and "[NH4+1]"
        # That's a bit more cumbersome, but it's easy to code
        smiles += "%+d" % charge

    return smiles + "]"

_bond_symbols = {
    BondType.SINGLE: "",
    BondType.AROMATIC: "",  # No need to use ":" to force aromatic perception since
                            # an implicit bond between two aromatics is aromatic
    BondType.DOUBLE: "=",
    BondType.TRIPLE: "#"
    }

def format_bond(bond):
    bondtype = bond.GetBondType()
    # Handle the special case where the end atoms are aromatic but the
    # bond itself is not aromatic. c1ccccc1-c1ccccc1 . This preserves
    # aromatic information for some ring systems when parsed by tools
    # like OEChem which don't atomatically perceive aromaticity on input.
    if (bondtype is BondType.SINGLE and
        bond.GetBeginAtom().GetIsAromatic() and
        bond.GetEndAtom().GetIsAromatic()):
        return "-"
    return _bond_symbols[bondtype]

#_closure_symbols = list("123456789") + ["%"+"%d"%i for i in range(10, 100)] + ["0"]
def format_closure(i, atoms):
    _closure_symbols = list("123456789") + ["%" + "%d" % i for i in range(10, (len(atoms) + 1))] + ["0"]
    return _closure_symbols[i]

# Caution: this version is limited to 100 bonds.
def reordered_smiles_100(mol, atoms):
    # Make a mapping from bond.GetIdx() to a unique closure string
    #_closure_symbols = list("123456789") + ["%" + "%d" % i for i in range(10, len(atoms) - 1)] + ["0"]
    closure_map = {}
    for bond_i, bond in enumerate(mol.GetBonds()):
        closure_map[bond.GetIdx()] = format_closure(bond_i, atoms)

    smiles = ""
    for atom_i, atom in enumerate(atoms):
        if atom_i != 0:
            # Each atom is dot disconnected
            # (Don't put a "." before the first atom)
            smiles += "."
        # Save the atom, in the correct order
        smiles += format_atom(atom)

        # Make the bond connection through the closure
        for bond in atom.GetBonds():
            smiles += format_bond(bond)
            smiles += closure_map[bond.GetIdx()]

    return smiles

target_kay = "CCSC1=C(N)N(N=C1C#N)C1=C(Cl)C=C(C=C1Cl)C(F)(F)F"

mol = Chem.MolFromSmiles(target_kay)
atoms = list(mol.GetAtoms())

# Generate some random permutations and verify matching SMILES
list_smiles = []
list_mols = []
for i in range(10):
    random.shuffle(atoms)
    new_smiles = reordered_smiles_100(mol, atoms)
    new_canonical_smiles = Chem.CanonSmiles(new_smiles)
    list_smiles.append(new_canonical_smiles)
    new_mol = Chem.MolFromSmiles(new_smiles)
    list_mols.append(new_mol)
    print(new_canonical_smiles)
    #assert target_kay == new_canonical_smiles

# mol_show = Draw.MolsToGridImage(list_mols, molsPerRow=3, subImgSize = (500, 500), returnPNG = False)
# print(mol_show)

mol_show = Draw.MolsToGridImage(list_mols, molsPerRow=3, subImgSize = (500, 500), returnPNG = True)
mol_show.save('output.png')
