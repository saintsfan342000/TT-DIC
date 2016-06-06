import numpy as n
import os
import glob
import pandas as pd
from scipy.interpolate import interp1d
os.getcwd()

expt = 12
frompath = r'F:\TensTors\TT2-12_FS32SS8\AllPts'
fromprefix = r'TT2-12_FS32SS8-Stage-0-' 
topath = r'F:\TT_New\TT2-12_FS32SS8\AramisExport_MissingRemoved'
toprefix = r'TT2-12_FS32SS8_'
BIN=True

LVpath = '{}/../../LVFiles/TT2-{:.0f}_LV.dat'.format(topath,expt)

last = 611

FMT='%.0f,%.0f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f'
headerline = '[0]Index_x [1]Index_y [2,3,4]Undef_X,Y,Z inches [5,6,7]Def_X,Y,Z inches [8,9,10,11]DefGrad (11 12 21 22)'

key = pd.read_excel('{}/../../TT-Summary.xlsx'.format(topath),sheetname='Summary',header=None,index_col=None,skiprows=1).values
R,th = key[ key[:,0] == expt, [3,4]].flatten()

STF = n.empty( (last+1,6) )

for i in range(last+1):
    
    d = n.genfromtxt('{}/{}{}.txt'.format(frompath,fromprefix,i),delimiter=',')
    
    if i != d[0,0]:
        print('Stage numbering inconsistent!')
        break
    else:
        STF[i,0] = d[0,0]
        STF[i,1] = d[0,1]
        
        if BIN == True:
            n.save('{}/{}{:.0f}'.format(topath,toprefix,STF[i,0]), d[ ~n.any( n.isnan(d),axis=1), :][1:])
        else:
            n.savetxt('{}/{}{:.0f}.dat'.format(topath,toprefix,STF[i,0]) , X=d[ ~n.any( n.isnan(d),axis=1), :][1:], fmt=FMT, header=headerline)
        print(i)
        
    
LV = n.loadtxt(LVpath,delimiter='\t')
STF[0,1] = LV[0,-1]
STF[:,[4,5]] = interp1d(LV[:,-1],LV[:,[2,0]],axis=0).__call__(STF[:,1])
STF[:,2] = STF[:,4]/(2*n.pi*R*th)/1000
STF[:,3] = STF[:,5]/(2*n.pi*R*R*th)/1000

n.savetxt('{}/STF.dat'.format(topath),X=STF,fmt='%.0f,%.6f,%.8f,%.8f,%.8f,%.8f',header='[0]Stg [1]Time [2]AxSts [3]ShSts [4]AxForce [5]Torque')