from pymol import cmd

def subseq():

    aaList = {'aa': []}
    cmd.iterate("name ca", "aa.append([resi, chain, resn])", space=aaList)
    print aaList

cmd.extend('subseq', subseq)
