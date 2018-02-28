from pymol import cmd
import re

cmd.set("seq_view", 1)

oneLetter = {
	'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
	'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
	'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
	'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
}

def subseq(target):
	# Gaunam bendros sekos aa pozicijas, pacias aa, grandine
	aaPosList = {'aa': []}
	cmd.iterate("name ca", "aa.append([resn, resi, chain])", space=aaPosList)

	# Mums reikia vienraidzio aa kodo, o pymol grazina triraidi koda
	for aa in aaPosList['aa']: aa[0] = oneLetter[aa[0]]

	# Susikonstruojam vientisa aa grandine (panasiai kaip fasta)
	aaComplete = ''
	for aa in aaPosList['aa']: aaComplete += aa[0]

	# RegEx validavimas
	try:
		_target = re.compile(target, re.I)		
	except Exception as e:
		print "Error: for target(RegEx):", target, "- syntax error."
		return

	matchObj = _target.search(aaComplete)
	while matchObj:
		print matchObj.start(), matchObj.end()
		matchObj = _target.search(aaComplete, matchObj.start() + 1)
	

	
	
	

cmd.extend('subseq', subseq)
