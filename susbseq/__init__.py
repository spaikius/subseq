from pymol import cmd

import os
import sys
import re

def __init__(self):
    path = os.path.dirname(__file__)
    sys.path.append(path)

    import ss_ARG
    import ss_DATA


def subseq(*_argv, **_kwargs):
    try:
        ss_ARG.set_parameters(_kwargs)
        param = ss_ARG.get_parameters()
        data = ss_DATA.get_data(param['models'], param['chains'])
        print(data)
    except Exception as e:
        print('Error: ' + str(e))

    ss_ARG.reset_parameters()
    ss_DATA.del_data()


cmd.extend('subseq', subseq)



# def subseq_re():
    # print('test')
    # _target = params['target']
    # _chains = params['chains']
    # print
    # _chains
    # RegExp validation
    # try:
    #     re_target = re.compile(_target, re.I)
    # except:
    #     raise Exception('Bad target(RegExp): ' + str(_target))

    # matchObj = _target.search(aaComplete)
    #  if matchObj:
    #   # Tuscias select'as, kuri papildysim.
    #   selectName = 'subseq'
    #   cmd.select(selectName, 'none')
    # else:
    #   print "No match found for given target: " + _target
    #   return

    # # Iteruojam, kol nebus rastas match'as
    # while matchObj:
    #   # Match'o pradzios ir pagaibos koordinates
    #   start = str(aaList['aa'][matchObj.start()][1])
    #   end = str(aaList['aa'][matchObj.end() - 1][1])

    #   # kokiai grandiniai priklauso match'o prima ir paskutine aa
    #   chainStart = str(aaList['aa'][matchObj.start()][2])
    #   chainEnd = str(aaList['aa'][matchObj.end() - 1][2])

    #   if chainStart == chainEnd:
    #       # Papildom select'a, o ne isnaujo perasom...
    #       cmd.select(selectName, selectName +  " | ( i. " + start + "-" + end + " & c. " + chainStart + ")")

    #   # kadangi re.search iesko tik pacio pirmo match'o, tai tesiant
    #   # pradedam nuo paskutinios rastos kooridnates
    #   # klausimas: ar leidziam match'ams overlap'int?
    #   matchObj = target.search(aaComplete, matchObj.start() + 1)

subseq_help = '''
@usage : subseq target=DIEVDLLKNGER, type=gl, chains=[A,C,D]
@@ Please note: arguments must be sepereted with comma (,) 
@params:
(required) target=(str)   : target sequence
(optional) model=[array]  : default all
(optional) type=(str)     : re or lo default re
(optional) chains=[array] : default all
'''