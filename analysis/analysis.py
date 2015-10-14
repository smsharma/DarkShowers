import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib import rcParams
import os, sys


class Data:
    """
    Class to store data contents (event, object and run variables) and do simple cuts etc
    """
    def __init__(self, input, display=True, lumi=20*1000, load_data=True):
        if load_data:
            if os.path.exists(input+".evt"): self.data_ = pd.read_csv(input+".evt", sep=',',skipinitialspace=True, comment="#", index_col=('evt') ) 
            if os.path.exists(input+".obj"): 
                self.obj_ = pd.read_csv(input+".obj", sep=',',skipinitialspace=True, comment="#", index_col=('evt','type','n') ) 
            if os.path.exists(input+".meta"): self.meta = pd.read_csv(input+".meta", sep=',',skipinitialspace=True, comment="#", index_col=None )
        self.w_=np.asarray(self.meta.at[0,'cxn']*lumi/self.meta.at[0,'nevt'])
        
        self.data = self.data_
        self.obj = self.obj_
        self.w =np.repeat(self.w_, len(self.data.index))
        self.cuts = Cuts()
        self.lumi = lumi
        
    def select(self, cuts):
        query = ''
        for cut in cuts._cuts:
            if query == '':
                query += cut
            else: query += '& ' + cut 
        if query != '': self.data = self.data_.query(query)
        else: self.data = self.data_
        evt_list = list(self.data.index)
        # evt_list = [s + 1 for s in list(self.data.index)]
        self.obj = self.obj_.loc[evt_list]
        self.w =np.repeat(self.w_, len(self.data.index))
        self.cuts = cuts

    def eff(self):
        return len(self.data)/np.float(len(self.data_))

    def n_evts(self, L=1):
        pb = 1;
        ab = 10**-6*pb;
        return self.meta['cxn'][0]*pb*(1/ab)*self.eff()
        
    def get_obj(self, type, n=1):
        return self.obj.xs((type, n), level=('type','n'))

    def etaphi_to_xyz(self,pt, mass, eta, phi):
        e = sqrt( (pt*cosh(eta))**2 + mass**2 )
        px = pt*cos(phi)
        py = pt*sin(phi)
        pz = pt*sinh(eta)
        return e, px, py, pz

    def dR(self,dphi, deta=0):
        dphi = abs(dphi % (2*pi))
        dphi = dphi if dphi < pi else dphi-pi
        return sqrt(dphi**2 + deta**2)

class Cuts:
    """
    Class to keep track of cuts
    """
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

class Plots:
    """
    Class to plot histograms and event displays
    """
    def __init__(self):
        # rcParams for plots 
        rcParams['text.usetex']=False
        rcParams['xtick.labelsize'] = 16
        rcParams['ytick.labelsize'] = 16
        rcParams['axes.labelsize']=18
        rcParams['axes.titlesize']=18
        rcParams['font.family']='serif'
        rcParams['font.serif']='CMU Serif'

    def hist_all(self,sig,bkgs,bkg_labels,*args,**kwargs):
        n_vars = len(sig.data.columns.values)
        n_cols = 2
        n_rows = np.int(np.ceil(n_vars/n_cols))
        fig, axs = plt.subplots(n_rows,n_cols)
        rcParams['figure.figsize'] = 14., 5*n_rows
        
        n_var = 0

        for i in range(n_rows):
            for j in range(n_cols):
                if n_var >= n_vars: continue
                var = sig.data.columns.values[n_var]
                binning = np.linspace(np.min(sig.data[var]),np.max(sig.data[var]),25)
                for bkg,bkg_label in zip(bkgs,bkg_labels):
                    axs[i][j].hist([list(bkg.data[var])], bins=binning, stacked=True,color="green", weights=[bkg.w],
                             log=True,
                             alpha=0.4, lw=0, label=bkg_label,*args,**kwargs)
                axs[i][j].hist(list(sig.data[var]), bins=binning, color="red", weights=sig.w,
                         alpha=1.0, lw=1.5, label=r"${\rm Signal}$",histtype="step",linestyle="dashed",log=True,*args,**kwargs)
                axs[i][j].set_xlabel(var, labelpad=10)
                axs[i][j].legend(loc='best', fancybox=True, framealpha=0.5)
                n_var+=1

    def hist_some(self,vars, sig,bkgs,bkg_labels,*args,**kwargs):
        n_vars = len(vars)
        n_cols = 2
        n_rows = np.int(np.ceil(n_vars/n_cols))
        fig, axs = plt.subplots(n_rows,n_cols)
        rcParams['figure.figsize'] = 14., 5*n_rows
        
        n_var = 0

        for i in range(n_rows):
            for j in range(n_cols):
                if n_var >= n_vars: continue
                var = vars[n_var]
                binning = np.linspace(np.min(sig.data[var]),np.max(sig.data[var]),25)
                for bkg,bkg_label in zip(bkgs,bkg_labels):
                    axs[i][j].hist([list(bkg.data[var])], bins=binning, stacked=True,color="green", weights=[bkg.w],
                             log=True,
                             alpha=0.4, lw=0, label=bkg_label,*args,**kwargs)
                axs[i][j].hist(list(sig.data[var]), bins=binning, color="red", weights=sig.w,
                         alpha=1.0, lw=1.5, label=r"${\rm Signal}$",histtype="step",linestyle="dashed",log=True,*args,**kwargs)
                axs[i][j].set_xlabel(var, labelpad=10)
                axs[i][j].legend(loc='best', fancybox=True, framealpha=0.5)
                n_var+=1

    def hist_var(self,var,sig,bkgs,bkg_labels,*args,**kwargs):
        binning = np.linspace(np.min(sig.data[var]),np.max(sig.data[var]),25)
        for bkg,bkg_label in zip(bkgs,bkg_labels):
            plt.hist([list(bkg.data[var])], bins=binning, stacked=True,color="green", weights=[bkg.w], log=True, alpha=0.4, lw=0, label=bkg_label,*args,**kwargs)
        plt.hist(list(sig.data[var]), bins=binning, color="red", weights=sig.w, alpha=1.0, lw=1.5, label=r"${\rm Signal}$",histtype="step",linestyle="dashed",log=True,*args,**kwargs)
        plt.xlabel(var, labelpad=10)
        plt.legend(loc='best', fancybox=True, framealpha=0.5)

    def event_display(self, dat, n_events=16):
        uq=[]
        for i in dat.obj.xs(('vis'),level=('type')).index.values:
            uq.append(i[0])
        uq=np.unique(uq)

        plt.rcParams['figure.figsize'] = 14, 14 
        plt.rcParams['axes.labelsize'] = 16
        plt.rcParams['font.size'] = 16

        fig, axs = plt.subplots(int(np.sqrt(n_events)),int(np.sqrt(n_events)))

        i_events=0;
        for i in range(int(np.sqrt(n_events))):
            for j in range(int(np.sqrt(n_events))):
            
                vislist=dat.obj.xs((uq[i_events],'vis'),level=('evt','type'))
                invislist=dat.obj.xs((uq[i_events],'invis'),level=('evt','type'))
                muonlist=dat.obj.xs((uq[i_events],'muon'),level=('evt','type'))
                jetlist=dat.obj.xs((uq[i_events],'jet'),level=('evt','type'))

                axs[i][j].scatter(vislist['eta'], vislist['phi'], s=10*vislist['pt'], c=[1,0,0],alpha=0.5)
                axs[i][j].scatter(invislist['eta'], invislist['phi'], s=10*invislist['pt'], c=[0,0,0], alpha=0.5)
                axs[i][j].scatter(muonlist['eta'], muonlist['phi'], s=10*muonlist['pt'], c=[0,0,1],alpha=0.5, marker='$\mu$' )
                axs[i][j].scatter(jetlist['eta'], jetlist['phi'], s=10*jetlist['pt'], c=[0,1,0],alpha=0.5, marker='$j$' )

                axs[i][j].set_xlabel('$\eta$')
                axs[i][j].set_ylabel('$\phi$')

                axs[i][j].set_xlim(-5,  5)
                axs[i][j].set_ylim(0, 2*np.pi)
         
                i_events+=1;
                
        plt.tight_layout()
        plt.show()
