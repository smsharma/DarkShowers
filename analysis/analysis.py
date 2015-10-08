import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib import rc

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

