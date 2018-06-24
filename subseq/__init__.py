import sys
import os

from pymol import cmd


def __init__(self):
    # Append path
    sys.path.append(os.path.dirname(__file__))
    
    # Module files
    import subseq_re
    import subseq_local_alignment
    import subseq_global_alignment


def subseq_re_dialog(app):
    pass


def subseq_local_alignment_dialog(app):
    pass


def subseq_global_alignment_dialog(app):
    pass


# PyMOL console extends
cmd.extend('subseq', subseq_re.subseq_re)
cmd.extend('subseq.local', subseq_local_alignment.subseq_local_alignment)
cmd.extend('subseq.global', subseq_global_alignment.subseq_global_alignment)
