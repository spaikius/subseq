---
title: "pymol subseq plug-in"
author: "Rimvydas Noreika"
date: "March 8, 2018"
---

## subseq plug-in installation guide
### - using pymol GUI
```
Plugin -> Plugin Manager -> Install New Plugin -> Choose file... -> PATH/TO/subseq.zip -> open -> choose your installion dir -> ok
```
## User guide
### Usage
```
subseq target=<str>, method=<str>, submatrix=<str>, models=<list>
       , chains=<list>, gapcost=<float>, minscore=<float>

Exaple usage: subseq target=KTGT, method=la, chains=[A, B, T], submatrix=PATH/TO/MATRIX
Please note: each keyword parameter should be seperated with comma (,)
```
### Parameters
```
Parameters:
    help                            ; Prints usage manual

Keyword parameters:
    target=<str>        Required    ; Target sequence
                                      Examples:
                                       If method type is re:
                                        - target=KTGTAVU
                                        - target=(TATA.{3,5}ATG(.{3,4}){3,})
                                       If method type is la:
                                        - target=KTGAT


    method=<str>        Optional    ; Method search type.
                                      - 're' for Regular Expression
                                      - 'la' for local alignment. Smith-Waterman 
                                      Default value: 're'

    models=<list>       Optional    ; The list of models that will be used for target search.
                                      If the list is not provided, then search for target will be 
                                      performed in all available models
                                      Example: models=[5ara, 2cif, a4s2]

    chains=<list>       Optional    ; The list of chains that will be used for target search.
                                      If the list is not priveded, then search for target will be
                                      performed in all available model chains
                                      Example: chains=[A, T, X, Q]

    submatrix=<PATH>    Optional    ; Path to substitution matrix for local alignment
                                      Default substitution matrix: Blossum62

    gapcost=<float>     Optional    ; The linear gap cost for local alignment
                                      Default value: 10

    minscore=<float>    Optional    ; The minimum score in precentages for throwing off low score alignments
                                      in local alignment search.
                                      Default value: 51.00 ( 51% )
                                      Example: minscore=75.25
```