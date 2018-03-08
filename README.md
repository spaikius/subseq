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
usage: 
```
subseq target=DIEVDLLKNGER, type=lo, chains=[A,C,D]
Please note: parameters must be sepereted with comma (,) 
```
params:
```
(required) | target=(str)   : target sequence
(optional) | model=[array]  : default all
(optional) | type=(str)     : re or lo default re
(optional) | chains=[array] : default all
```