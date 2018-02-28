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
	for aa in aaPosList['aa']: aaComplete += str(aa[0])

	# RegEx validavimas
	try:
		_target = re.compile(target, re.I)		
	except Exception as e:
		print "Error: for target(RegEx):" + target + "- syntax error."
		return

	matchObj = _target.search(aaComplete)

	# Tuscias select'as, kuri papildysim.
	selectName = 'subseq'
	cmd.select(selectName, 'none')

	while matchObj:
		start = str(aaPosList['aa'][matchObj.start()][1])
		end = str(aaPosList['aa'][matchObj.end() - 1][1])

		chainStart = str(aaPosList['aa'][matchObj.start()][2])
		chainEnd = str(aaPosList['aa'][matchObj.end() - 1][2])

		if chainStart == chainEnd:
			cmd.select(selectName, selectName +  " | ( i. " + start + "-" + end + " & c. " + chainStart + ")")
		
		matchObj = _target.search(aaComplete, matchObj.start() + 1)
	

	
	
	

cmd.extend('subseq', subseq)
