import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.animation as animation
import numpy as np
from numpy.random import normal as normal
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib
from astropy import constants as const
from astropy import units as u
from astropy.visualization import quantity_support
from astropy.constants import G
from astropy.io import ascii
from tqdm import tqdm
import pickle
import struct
from astropy.table import Table

class snapshot:
    """Each time snapshot in the nbody6 run.
    """

    def __init__(self,time,header,M,X,V,N):
        """__init__.

        Parameters
        ----------
        time :
            time of the snapshot in Myr
        header :
            header of the OUT3 file 
        M :
            Masses of each object in the snapshot 
        X :
            Pos of each object in the snapshot
        V :
            Vel of each object in the snapshot
        N :
            Name of each object in the snapshot
        """
        self.time=time
        self.o3_header=header
        self.out3=Table([M, X, V, N],names=('Mass', 'X', 'V', 'name'),meta={'name': 'OUT3'})

        
class nb6run:
    """NBody6 run
    """

    def __init__(self,PATH):
        """__init__.

        Parameters
        ----------
        PATH :
            Path of the Output files of Nbody6
        """
        self.PATH=PATH
        n,headers,Ms,Xs,Vs,Ns=self.read_out3()
        self.snaps=self.snapshots_out3(n,headers,Ms,Xs,Vs,Ns)
        n,time_arr,headers,E_ks,ks_bindata_snp,Npair_arr=self.read_out9()
        self.snapshots_out9(n,time_arr,headers,E_ks,ks_bindata_snp,Npair_arr)
        del ks_bindata_snp
        
        
    def snapshots_out3(self,n,headers,Ms,Xs,Vs,Ns):
        """Creates snapshots out of the OUT3 data.

        Parameters
        ----------
        n :
            number of snapshots
        headers :
            array of headers of OUT3 
        Ms :
            array of masses in each snapshot
        Xs :
            array of positions of each object in each snapshot
        Vs :
            array of velocities of each object in each snapshot
        Ns :
            array of names of each object in each snapshot
        """
        snap_array=[]
        for i in range(0,n):
            time=headers[i][3]
            ntot=headers[i][0]
            X_arr=Xs[i]
            X_arr=np.reshape(X_arr,(ntot,3))
            V_arr=Vs[i]
            V_arr=np.reshape(V_arr,(ntot,3))
            snap_array.append(snapshot(time,headers[i],Ms[i],X_arr,V_arr,Ns[i]))
        return snap_array
    
    def snapshots_out9(self,n,time_arr,headers,E_ks,ks_bindata_snp,Npair_arr):
        """snapshots_out9.

        Parameters
        ----------
        n :
            n
        time_arr :
            time_arr
        headers :
            headers
        E_ks :
            E_ks
        ks_bindata_snp :
            ks_bindata_snp
        Npair_arr :
            Npair_arr
        """
        for i in range(0,n):
            if time_arr[i]*self.snaps[i].o3_header[4][10]==self.snaps[i].time:
                self.snaps[i].o9_header=headers[i]
                self.snaps[i].E_k=E_ks[i]
                self.snaps[i].ks_bindata=ks_bindata_snp[i]
                self.snaps[i].NPAIRS=Npair_arr[i]
            else:
                raise Exception("time of snapshots do not match")

            
    def read_out9(self):
        """read_out9.
        """
        filename=self.PATH+'OUT9'
        with open(filename, "r") as f:
            #Information in bindat.f
            n=0
            headers=[]
            ks_bindata_snp=[]
            E_ks=[]
            Npair_arr=[]
            time_arr=[]
            while True:
                try:
                    h1=ascii.read(f.readline(),format='no_header',guess=False,fast_reader=True)
                except:
                    break #End of file
                NPAIRS=h1[0][0]
                TIME=h1[0][6]
                header=list(h1[0])
                h2=ascii.read(f.readline(),format='no_header',guess=False,fast_reader=True)
                E_k=list(h2[0])
                h3=ascii.read(f.readline(),format='no_header',guess=False,fast_reader=True)
                header=header+list(h3[0])
                temp=[]
                for t in range(0,NPAIRS):
                    temp.append(f.readline())
                ks_bindata=ascii.read(temp,format='no_header',guess=False,fast_reader=True,names=('E_b','ecc','E_cm','r','m1','m2','P','name1','name2','Stype1','Stype2','merge_type'))
                n=n+1
                headers.append(header)
                ks_bindata_snp.append(ks_bindata)
                E_ks.append(E_k)
                Npair_arr.append(NPAIRS)
                time_arr.append(TIME)
        return n,time_arr,headers,E_ks,ks_bindata_snp,Npair_arr
    
    
    def read_out3(self):
        """read_out3.
        """
        #Info in output.f
        filename=self.PATH+'OUT3'
        with open(filename, "rb") as f:
            n=0
            headers=[]
            Ms=[]
            Xs=[]
            Vs=[]
            Ns=[]
            while True:
                try:
                    bod1=struct.unpack('i',f.read(4))[0]
                except:
                    break #End of file
                NTOT,MODEL,NRUN,NK=struct.unpack('i'*4,f.read(4*4))
                bod2=struct.unpack('i',f.read(4))[0]
                if bod1==bod2:
                    pass
                else:
                    raise Exception("BODY SIZE ERROR")
                bod1=struct.unpack('i',f.read(4))[0]
                AS=struct.unpack('f'*NK,f.read(4*NK))
                time=AS[0]*AS[10]
                M=struct.unpack('f'*NTOT,f.read(4*NTOT))
                X=struct.unpack('f'*NTOT*3,f.read(4*NTOT*3))
                V=struct.unpack('f'*NTOT*3,f.read(4*NTOT*3))
                N=struct.unpack('i'*NTOT,f.read(4*NTOT))
                bod2=struct.unpack('i',f.read(4))[0]
                if bod1==bod2:
                    pass
                else:
                    raise Exception("BODY SIZE ERROR")
                n=n+1
                header=[NTOT,MODEL,NRUN,time,AS]
                headers.append(header)
                Ms.append(M)
                Xs.append(X)
                Vs.append(V)
                Ns.append(N)
            return n,headers,Ms,Xs,Vs,Ns

