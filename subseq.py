from pymol import cmd

cmd.set("seq_view", 1)

oneLetter = {
	'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
	'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
	'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
	'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
}

def subseq():
	# Gaunam bendros sekos aa, pozicija, grandine
	aaList = {'aa': []}
	cmd.iterate("name ca", "aa.append([resn, resi, chain])", space=aaList)

	# Mums reikia vienraidzio aa kodo, o pymol grazina triraidi koda
	for aa in aaList['aa']: aa[0] = oneLetter[aa[0]]

	print aaList
	
	

cmd.extend('subseq', subseq)
