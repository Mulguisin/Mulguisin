<div align="center">
   <!-- MGS logo -->
   <img src="./MGS_logo.png" width=300>
</div>

# Mulguisin

The Mulguisin(MGS) is a cluster/group finding algorithm in galaxy data. 
It can find cluster/group structure from the galaxy data and also provide topological informations.

It has been developed by InKyu Park to find jet structure for the LHC experiment.
Now, he modify the Mulguisin that can be applied to galaxy data.

## Install

The MGS can be install via pip install: 


```
pip install Mulguisin
```

### Source code

The latest version of MGS can be download with: 

```
git clone https://github.com/youngju20/Mulguisin.git
```

## Basic usage

```python
>>> import mgs_class
>>> import numpy as np
>>> data = np.random.random((10,2))
>>> boundaries = [0,1,0,1]
>>> Rcut = 0.1
>>> MGS = mgs_class.mulguisin(Rcut,data[:,0],data[:,1],boundaries=boundaries)
>>> Nmgs, imgs, clg, clm, cng = MGS.get_mgs()
>>> Nmgs
7
```

More examples are in the [test](test). We provide several jupyter notebooks to test MGS code.

