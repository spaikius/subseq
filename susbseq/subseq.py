from pymol import cmd
from ast import literal_eval

import re

cmd.set("seq_view", 1)

type_vals = ['re', 'gl']

params = {
	'type'  : type_vals[0],
	'target': None,
	'chains': None,
}

aaList = {'aa': []}
complete_chains = dict()

def subseq(**_kwargs):
	_kwargs.pop('_self', None)

	if bool(_kwargs):
		try:
			validate_args(_kwargs)
			init_data()
		except Exception as e:
			print('Error: ' + str(e))
			# print e
			return

# Function validates and sets parameters
def validate_args(_kwargs):
	for key, value in _kwargs.iteritems():
		# change to lower case
		(key, value) = (key.lower(), value.lower())

		if key in params:
			if key == 'type' and value in type_vals:
				params[key] = value
			elif key == 'target':
				params[key] = value
			elif key == 'chains':
				_list = literal_eval(re.sub('([A-Za-z])', r'"\1"', value))
				params[key] = [ i.upper() for i in _list ]
			else:
				raise Exception('Unknown value in ' + key + ' : ' + value)
		else:
			raise Exception('Unknown parameter: ' + key)

def init_data():
	cmd.iterate("name ca", "aa.append([resn, resi, chain, model])", space=aaList)
	
	# for aa in aaList['aa']: aa[0] = oneLetter[aa[0]]

	for name in cmd.get_names():
		complete_chains[name] = {}
		for chain_name in cmd.get_chains(name):
			complete_chains[name][chain_name] = ''

	for i in aaList['aa']:
		complete_chains[i[3]][i[2]] += str(i[0])

	print complete_chains


def subseq_re(_target):

	# RegEx validavimas
	try:
		target = re.compile(_target, re.I)		
	except Exception as e:
		print "Error: for target(RegEx): " + _target + " - syntax error."
		return


	matchObj = _target.search(aaComplete)
	if matchObj:
		# Tuscias select'as, kuri papildysim.
		selectName = 'subseq'
		cmd.select(selectName, 'none')
	else:
		print "No match found for given target: " + _target
		return

	# Iteruojam, kol nebus rastas match'as
	while matchObj:
		# Match'o pradzios ir pagaibos koordinates
		start = str(aaList['aa'][matchObj.start()][1])
		end = str(aaList['aa'][matchObj.end() - 1][1])

		# kokiai grandiniai priklauso match'o prima ir paskutine aa
		chainStart = str(aaList['aa'][matchObj.start()][2])
		chainEnd = str(aaList['aa'][matchObj.end() - 1][2])

		if chainStart == chainEnd:
			# Papildom select'a, o ne isnaujo perasom...
			cmd.select(selectName, selectName +  " | ( i. " + start + "-" + end + " & c. " + chainStart + ")")
		
		# kadangi re.search iesko tik pacio pirmo match'o, tai tesiant
		# pradedam nuo paskutinios rastos kooridnates
		# klausimas: ar leidziam match'ams overlap'int?
		matchObj = target.search(aaComplete, matchObj.start() + 1)
	

cmd.extend('subseq', subseq)

oneLetter = {
	'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
	'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
	'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
	'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
}

