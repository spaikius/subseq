---
title: "pymol subseq plug-in"
author: "Rimvydas Noreika"
date: "March 8, 2018"
---

## subseq plug-in installation guide
Using PyMOL Plugin Manager, follow these steps: 
```
1) Press Plugin in the menu bar 
2) Press Plugin Manager
3) Press Install New Plugin 
4) Press Choose file... 
5) Navigate to subseq.zip 
6) Press open
7) Choose where to install plugin
8) Press ok
```
## subseq.py script guide
Type in PyMOL console
```
run PATH/TO/subseq.py
```

Or configurate your pymolrc file
On Windows: Start > Run and then paste
notepad "%HOMEDRIVE%%HOMEPATH%\pymolrc.pml"

On Unix/Linux-type system (including Mac OS X): Open a terminal and type
nano ~/.pymolrc

## Help
```
For RegExp search method type: help subseq
For local alignment search method type: help subseq.local
For global alignment search method type: help subseq.global
```
