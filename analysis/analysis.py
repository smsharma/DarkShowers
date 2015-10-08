import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib import rc
import os, sys


class Data:
    def __init__(self, input, display=True, lumi=20*1000):
        self.data_ = pd.read_csv(input+".evt", sep=',',skipinitialspace=True, comment="#", index_col=('evt') ) 
        self.obj_ = pd.read_csv(input+".obj", sep=',',skipinitialspace=True, comment="#", index_col=('evt','type','n') ) 
        self.meta = pd.read_csv(input+".meta", sep=',',skipinitialspace=True, comment="#", index_col=None )
        self.w_=np.asarray(self.meta.at[0,'cxn']*lumi/self.meta.at[0,'nevt'])
        
        self.data = self.data_
        self.obj = self.obj_
        self.w =np.repeat(self.w_, len(self.data.index))
        
    def select(self, criterion):
        self.data = self.data_.query(criterion)
        evt_list = list(self.data.index)
        self.obj = self.obj_.loc[evt_list]
        self.w =np.repeat(self.w_, len(self.data.index))
        
    def get_obj(self, type, n=1):
        return self.obj.xs((type, n), level=('type','n'))

# a class to keep track of cuts
class Cuts:
    def __init__(self):
        self._cuts = []
        
    # pretend Cuts is a list
    def __len__(self): return len(self._cuts)
    def __getitem__(self, key): return self._cuts[key]
    def __setitem__(self, key, value): self._cuts[key] = value
    def __delitem__(self, key): del self._cuts[key]
    def append(self, *args, **kwargs): return self._cuts.append(*args, **kwargs)
    def extend(self, *args, **kwargs): return self._cuts.extend(*args, **kwargs)
    def count(self, *args, **kwargs): return self._cuts.count(*args, **kwargs)
    def index(self, *args, **kwargs): return self._cuts.index(*args, **kwargs)
    def insert(self, *args, **kwargs): return self._cuts.insert(*args, **kwargs)
    def remove(self, *args, **kwargs): return self._cuts.remove(*args, **kwargs)
    def pop(self, *args, **kwargs): return self._cuts.pop(*args, **kwargs)
    def reverse(self, *args, **kwargs): return self._cuts.reverse(*args, **kwargs)
    def sort(self, *args, **kwargs): return self._cuts.sort(*args, **kwargs)

    # display cuts in a useful way
    def __repr__(self):
        if len(self._cuts) == 0:
            return "No cuts have been defined"
        output = ["Requiring all of the following:"]
        for i, c in enumerate(self._cuts):
            output.append("    cuts[%d] \"%s\"" % (i, c))
        return os.linesep.join(output)

    # give a cut string to ROOT
    def __str__(self):
        if len(self._cuts) == 0:
            return ""
        else:
            return "(" + ") && (".join(self._cuts) + ")"

    # copy this set of cuts
    def copy(self):
        output = Cuts()
        output._cuts = self._cuts[:]
        return output

    # forget all cuts
    def clear(self): self._cuts = []

# class Plots:
#     def __init__(self):
#         self._cuts = []

#     def histogram_all(self):
#     def histogram_some(self):
#     def event_display(self):
#     def histogram_all(self):
#     def histogram_all(self):