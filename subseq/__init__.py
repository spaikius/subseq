import sys
import os

from pymol import cmd

sys.path.append(os.path.dirname(__file__))

import subseq_re
import subseq_local_alignment
import subseq_global_alignment

def __init__(self):
    pass

cmd.extend('subseq', subseq_re.subseq_re)
cmd.extend('subseq.local', subseq_local_alignment.subseq_local_alignment)
cmd.extend('subseq.global', subseq_global_alignment.subseq_global_alignment)
