import numpy as n
vstack, hstack = n.vstack, n.hstack
linspace = n.linspace
array = n.array
nanmean, nanstd = n.nanmean, n.nanstd
abs = n.abs
import matplotlib.pyplot as p
from scipy.interpolate import griddata, interp1d
from scipy.linalg import sqrtm, eigvalsh
from scipy.spatial.distance import pdist
from pandas import read_excel, read_csv
import os

'''
Designed to be the complete TT analysis code.
Does the same analysis as Nico's code, but more:
-Techincal strains are calcualted as done by Nico, but also uses Dr. K's definitions.  
-Record the aramis x,y indices of the max point in every stage; written to max.dat
-Profiles pass thru the last stage max point; to accomodate this the stages are processed in reverse order
-Computes rotation and delta/L
-profStgs are determined from a plot of torque vs MTS rotation, and
taken as equally-spaced increments of rotation rather than stage number.
**
Requires a file ArSTF, which is the Aramis stage, force, time output.
This is used for the interpolation of the Labview data.
**
Still to add:  
--Radial contraction!
--Max point's strain path
--Epeq vs rot

'''

###### PRELIMINARY DATA to ENTER ######
last = 580
expt=21
tube=16
FS = 19
SS = 6

arampath = r'F:\TensTors\TT2-21_DC15_FS19_SS6\AllPts'
prefix = 'TT2-21_FS19SS6-Stage-0-'
BIN = False      # IF Bin is true, then it will load npy files from arampath and prefix
                # Otherwise, it will used read_excel to openup .txt files

savepath = r'F:\TT_New\TT2-21_FS19SS6'  #Folder will results are going
makecontours = True    #Make contours every 50 stages?
saveAram = True    # Whether to save the missing removed array to a .npy binary
saveprefix = 'TT2-23_FS19SS6_'  #Prefix for the npy files

key = read_excel('{}/../TT-Summary.xlsx'.format(savepath),sheetname='Summary',header=None,index_col=None,skiprows=1).values
alpha, Rm,thickness = key[ key[:,0] == expt, [1,3,4]].flatten()

os.chdir(savepath)
Size_Scott = 1/16 #[in]

###########################################################################
# Make STF file
###########################################################################
if not os.path.exists('STF.dat'):
    try:
        ArSTF = read_csv('{}/../ArSTF.dat'.format(arampath),sep=',',header=None,comment='#').values
    except:
        ArSTF = read_csv('F:\TensTors\TT2-23_DC16_FS19_SS6\ArSTF.dat',sep=',',header=None).values
    LV = read_csv('./../LVFiles/TT2-{:.0f}_LV.dat'.format(expt),sep='\t',header=None,comment='#').values
    STF = n.empty( (len(ArSTF[:,0]),6) )
    STF[:,0] = ArSTF[:,0]
    STF[:,1] = ArSTF[:,2]
    STF[0,1] = LV[0,-1]
    TRFD = interp1d(LV[:,-1],LV[:,[0,1,2,3]],axis=0).__call__(STF[:,1])
    #STF[:,[4,5]] = interp1d(LV[:,-1],LV[:,[2,0]],axis=0).__call__(STF[:,1])
    STF[:,[4,5]] = TRFD[:,[2,0]]
    STF[:,2] = STF[:,4]/(2*n.pi*Rm*thickness)/1000
    STF[:,3] = STF[:,5]/(2*n.pi*Rm*Rm*thickness)/1000
    n.savetxt('STF.dat',X=STF,fmt='%.0f, %.6f, %.6f, %.6f, %.6f, %.6f')


###########################################################################
# Make the contour in which we select the area we analyze
###########################################################################

if BIN == True:
    A = n.load('{}/{}{:.0f}.npy'.format(arampath,prefix,last))
else:
    A = read_csv('{}/{}{:.0f}.txt'.format(arampath,prefix,last),sep=',',na_values=' ',skiprows=3,header=None,index_col=None).values
    A = A[ ~n.any(n.isnan(A),axis=1), :]

    #[0]Index_x [1]Index_y [2,3,4]Undef_X,Y,Z inches [5,6,7]Def_X,Y,Z inches [8,9,10,11]DefGrad (11 12 21 22) *)
A = A[ n.abs(A[:,3])<0.2, :]        # Reduces the amount we'll have to plot

if not os.path.exists('box_limits.dat'):
    x=A[:,2]
    y=A[:,3]

    #Calculation for each facet of the Logarithmic cumulative plastic strain
    LEp = n.empty(len(x))
    for i in range(len(x)):
        F = A[i,8:].reshape(2,2)
        U = sqrtm( (F.T).dot(F) )   #Stretch Tensor
        eigU = eigvalsh(U)
        #logarithmic strain in the principal coordinate system
        LE = n.log(eigU)
        #Logarithmic Cumulative Plastic strain
        LEp[i] = ( 2/3 * ( LE[0]**2 + LE[1]**2 + (-LE[0]-LE[1])**2 ) )**0.5

    xspace=linspace(n.min(x),n.max(x),len( n.unique(A[:,0]) )/2 )
    yspace=linspace(n.min(y),n.max(y),len( n.unique(A[:,1]) )/2 )

    LEP = griddata( vstack((x,y)).T , LEp, (xspace[None,:],yspace[:,None]), method='linear')
    meanLEP = nanmean(LEP.flatten())
    stdLEP = nanstd(LEP.flatten())

    p.close('all')
    p.figure(1,facecolor='w',figsize=(16,12))
    p.contourf(xspace,yspace,LEP,linspace(meanLEP-0*stdLEP,meanLEP+3*stdLEP,256),extend='both')
    p.axis([n.min(x), n.max(x), n.min(y), n.max(y)])
    #p.axis('equal')
    p.grid(True,linestyle='-',linewidth=2)
    p.xlabel('Undeformed X (in)')
    p.ylabel('Undeformed Y (in)')

    crop = n.asarray( p.ginput(2) )
    Xmin = n.min( crop[:,0] )
    Xmax = n.max( crop[:,0] )
    Ymin = n.min( crop[:,1] )
    Ymax = n.max( crop[:,1] )
    
    headerline='Xmin, Xmax, Ymin, Ymax (all in inches!)'
    n.savetxt('box_limits.dat',X=array([Xmin, Xmax, Ymin, Ymax])[None,:],fmt='%.6f',delimiter=', ',header=headerline)
    
    box_x = array([Xmin , Xmax , Xmax , Xmin , Xmin])
    box_y = array([Ymin , Ymin , Ymax , Ymax , Ymin])

    p.plot(box_x,box_y,'w',linewidth=2.5)
    p.grid(False)
    titlestring = 'TT2-{:.0f}, $\\alpha$ = {}, FS{:.0f}SS{:.0f}'.format(expt,alpha,FS,SS)
    p.title(titlestring)
    p.draw()
    p.savefig('last_stage.png')
    p.close()
else:
    Xmin,Xmax,Ymin,Ymax = n.loadtxt('box_limits.dat',delimiter=',')
    box_x = array([Xmin , Xmax , Xmax , Xmin , Xmin])
    box_y = array([Ymin , Ymin , Ymax , Ymax , Ymin])

A = A[ (abs(A[:,2])<=0.2) & (abs(A[:,3])<=0.1) , :]
min_disp = n.min( pdist(A[:,[2,3,4]]) )
SS_th = (min_disp/thickness)
FS_th = FS / SS * (min_disp/thickness)

#Definition of the size of the averaging zone for grid method
size_av = Size_Scott - thickness * FS_th

#Stage Time Force data...Create the 10 stations at which we'll plot profiles
STF = n.loadtxt('STF.dat',delimiter=',')           
STF = STF[0:last+1,:]               #In case I exported past the last stage
LL = n.max(n.argmax(STF[:,2:4], axis=0))

if not os.path.exists('prof_stages.dat'):
    p.figure()
    p.plot(TRFD[:,1],TRFD[:,0])
    p.plot(TRFD[LL,1],TRFD[LL,0],'ro',ms=10)
    p.title('Click where the first\ncalc stage will be', fontsize=30)
    #p.axis([0,1.1*TRFD[-1,1],0,1.2*n.max(TRFD[:,0])])
    p.gca().set_ylim(top=1.2*n.max(TRFD[:,0]))
    yldrot = n.asarray(p.ginput(1)).flatten()[0]
    p.close()

    # Make evenly spaced increments
    for i in range(len(TRFD[:,1])):
        if TRFD[i,1] >= yldrot:
            break

    #Two before the limit load, seven after
    x = hstack( (n.linspace(yldrot,TRFD[LL,1],3)[0:2], n.linspace(TRFD[LL,1],TRFD[-1,1],8)) )
    profStg = n.empty(len(x))
    for i in range(len(x)):
        profStg[i] = n.where(TRFD[:,1]>=x[i])[0][0]
        
    headerline = 'Stages at which profiles were generated'
    n.savetxt('prof_stages.dat',X=profStg[None,:],fmt='%.0f',delimiter=', ',header=headerline)
else:
    profStg = n.loadtxt('prof_stages.dat',delimiter=',')

if not os.path.exists('disp-rot_limits.dat'):
    ## Now identify the upper and lower ranges for calculating delta and ph
    
    if BIN == True:
        A = n.load('{}/{}{:.0f}.npy'.format(arampath,prefix,0))
    else:
        A = read_csv('{}/{}{:.0f}.txt'.format(arampath,prefix,0),sep=',',na_values=' ',skiprows=3,header=None,index_col=None).values
        A = A[ ~n.any(n.isnan(A),axis=1), :]
    #[0]Index_x [1]Index_y [2,3,4]Undef_X,Y,Z inches [5,6,7]Def_X,Y,Z inches [8,9,10,11]DefGrad (11 12 21 22) *)
    p.figure(figsize=(16,12))
    p.title('Click the four points that bound the cusp of the thick edges',size=20)
    p.plot(A[:,3],A[:,4],'.',alpha=0.5,ms=4,markevery=2)
    p.gcf().gca().set_ylim([.98*1.826/2,1.02*1.9685/2])
    p.gcf().gca().set_xlim(right=-.1)
    rdlim = n.sort(n.asarray(p.ginput(2))[:,0])
    p.close()
    p.figure(figsize=(16,12))
    p.plot(A[:,3],A[:,4],'.',alpha=0.5,ms=4,markevery=2)
    p.gcf().gca().set_ylim([.98*1.826/2,1.02*1.9685/2])
    p.gcf().gca().set_xlim(left=.1)
    rdlim = n.sort( hstack( (rdlim, n.asarray(p.ginput(2))[:,0]) ) )
    p.close()
    headerline='Lower sxn Ymin, Lower sxn Ymax, Upper sxn Ymin, Upper sxn Ymax'
    n.savetxt('disp-rot_limits.dat',X=rdlim[None,:],fmt='%.6f',delimiter=', ',header=headerline)
else:
    rdlim = n.sort(n.loadtxt('disp-rot_limits.dat',delimiter=','))

#Initialize a few things before looping and calculating every stage
export_max=n.zeros( (last+1,13) )   #MaxPt data
export_mean=n.zeros( (last+1,12) )  #MeanPt data
export_stdv=n.zeros( (last+1,12) )  #Std Deviation data
export_Scott=n.zeros( (last+1,12) ) #DIC Macro data
prof_count=0                        #For indexing the profLEp list        
profLEp = [0]*len(profStg)          #List of proLEp arrays
L = n.zeros(last+1)                   #For delta/L
up_th = n.zeros(last+1)               #For phi
lo_th = n.zeros(last+1)               #For phi

#Cycle through the stages
z=1
for k in range(last,-1,-1):
    print(k)
    if k == 0:
        if BIN:
            A = n.load('{}/{}{:.0f}.npy'.format(arampath,prefix,k))
        else:
            A = read_csv('{}/{}{:.0f}.txt'.format(arampath,prefix,k),sep=',',na_values=' ',skiprows=3,header=None,index_col=None).values
            A = A[ ~n.any(n.isnan(A),axis=1), :]
            if saveAram:
                if not os.path.exists('./AramisBinary'):
                    os.mkdir('./AramisBinary')
                n.save('./AramisBinary/{}{:.0f}'.format(saveprefix,k),A)
        export_max[k] =  [STF[k,0], STF[k,1], STF[k,2], STF[k,3], 0,0,0,0,0,0,0,0,0]
        export_mean[k] = [STF[k,0], 0, STF[k,1], STF[k,2], STF[k,3], 0,0,0,0,0,0,0]
        export_stdv[k] = [STF[k,0], 0, STF[k,1], STF[k,2], STF[k,3], 0,0,0,0,0,0,0]
        export_Scott[k] = [STF[k,0], 0, STF[k,1], STF[k,2], STF[k,3], 0,0,0,0,0,0,0]
    else:
        LEp, NEx, NEy, NExy, gamma, xcoord, ycoord, aramX, aramY, NEx_alt, NEy_alt, gamma_alt = ( [] for _ in range(12) )
        scot_count=0                #Count the number of points in the Scott method
        if BIN:
            A = n.load('{}/{}{:.0f}.npy'.format(arampath,prefix,k))
        else:
            A = read_csv('{}/{}{:.0f}.txt'.format(arampath,prefix,k),sep=',',na_values=' ',skiprows=3,header=None,index_col=None).values
            A = A[ ~n.any(n.isnan(A),axis=1), :]
            if saveAram:
                if not os.path.exists('./AramisBinary'):
                    os.mkdir('./AramisBinary')
                n.save('./AramisBinary/{}{:.0f}'.format(saveprefix,k),A)
        #[0]Index_x [1]Index_y [2,3,4]Undef_X,Y,Z inches [5,6,7]Def_X,Y,Z inches [8,9,10,11]DefGrad (11 12 21 22) *)
        Ascan = A[ (A[:,2]>=Xmin) & (A[:,2]<=Xmax) & (A[:,3]>=Ymin) & (A[:,3]<=Ymax), :]
        Xind = Ascan[:,0]
        for j in n.unique(Xind):                  #<-- transpose because this must be a row vector for the for loop to incrememnt correctly!
            Yind = Ascan[ Ascan[:,0] == j, 1]
            if len(Yind) == (n.max(Yind) - n.min(Yind) +1):

                # Intialization of col_ variables using generator
                colLEp, colNEx, colNEy, colNExy, colG, colx, coly, colindx, colindy = ( n.zeros( len(Yind) ) for _ in range(9) )
                colNEx_alt, colNEy_alt, colG_alt = (n.zeros( len(Yind) ) for _ in range(3) )

                for m in range( len(Yind) ):
                    Acol = Ascan[ (Ascan[:,0] == j) & (Ascan[:,1]== Yind[m]), :].flatten()
                    F = Acol[8:].reshape(2,2)
                    U = sqrtm( (F.T).dot(F) )   
                    eigU = eigvalsh(U)
                    LE = n.log(eigU)
                    colLEp[m] = ( 2/3 * ( LE[0]**2 + LE[1]**2 + (-LE[0]-LE[1])**2 ) )**0.5
                    #R = F * U^(-0.5)  #Rotation tensor
                    #NE_calc = R.'*(U^0.5 - eye([2 2]))*R
                    NE_calc = U- n.eye(2)
                    colNEx[m] = NE_calc[0,0]
                    colNEy[m] = NE_calc[1,1]
                    colNExy[m] = NE_calc[0,1]
                    colG[m] = n.arctan(colNExy[m]/(1+colNEx[m])) + n.arctan(colNExy[m]/(1+colNEy[m]))
                    colNEx_alt[m] = F[0,0]-1
                    colNEy_alt[m] = F[1,1]-1
                    colG_alt[m]=n.arctan(F[0,1]/F[1,1]);
                    colx[m], coly[m] = Acol[[2,3]]
                    colindx[m], colindy[m] = Acol[[0,1]]
                
                locLEp = n.argmax(colLEp)                #Location of...
                LEp.append( colLEp[locLEp] )             #Max LEp in the current column
                NEx.append( colNEx[locLEp] )          
                NEy.append( colNEy[locLEp] )
                gamma.append( colG[locLEp] )      
                xcoord.append( colx[locLEp] )
                ycoord.append( coly[locLEp] )
                aramX.append( colindx[locLEp] )
                aramY.append( colindy[locLEp] )
                NEx_alt.append( colNEx_alt[locLEp] )
                NEy_alt.append( colNEy_alt[locLEp] )
                gamma_alt.append( colG_alt[locLEp] )

        LEp, NEx, NEy, gamma, xcoord, ycoord, aramX, aramY = map(array,[LEp, NEx, NEy, gamma, xcoord, ycoord, aramX, aramY])    #Convert lists to arrays
        NEx_alt, NEy_alt, gamma_alt = map(array,[NEx_alt, NEy_alt, gamma_alt])
        ratio = NEy / gamma
        ratioAvg = nanmean(ratio)
        ratioSDEV = nanstd(ratio)
        passed = (ratio >= ratioAvg - 0.5 * ratioSDEV) & (ratio <= ratioAvg + 0.5 * ratioSDEV)
        LEp=LEp[passed]
        NEx=NEx[passed]
        NEy=NEy[passed]
        gamma=gamma[passed]
        NEx_alt=NEx_alt[passed]
        NEy_alt=NEy_alt[passed]
        gamma_alt=gamma_alt[passed]        
        xcoord=xcoord[passed]
        ycoord=ycoord[passed]
        aramX = aramX[passed]
        aramY = aramY[passed]
                
        locmax = n.argmax( LEp )
        if k == last:
            aramXmaxlast, aramYmaxlast = aramX[locmax], aramY[locmax]   #Save for making strain profiles
            MaxTen = n.flipud(array( vstack( (LEp,aramX,aramY,xcoord,ycoord) ) ).T[n.argsort(LEp),:][-10:,:])
                        
        
        export_max[k] =  [ STF[k,0], STF[k,1], STF[k,2], STF[k,3], NEx[locmax], NEy[locmax], abs(gamma[locmax]), 
                                NEx_alt[locmax], NEy_alt[locmax], abs(gamma_alt[locmax]), LEp[locmax], aramX[locmax], aramY[locmax] ]
        export_mean[k] = [ STF[k,0],  STF[k,1], sum(passed), STF[k,2], STF[k,3], nanmean(NEx), nanmean(NEy), abs(nanmean(gamma)), 
                                nanmean(NEx_alt), nanmean(NEy_alt), abs(nanmean(gamma_alt)), nanmean(LEp)]
        export_stdv[k] = [ STF[k,0], STF[k,1], sum(passed), STF[k,2], STF[k,3], nanstd(NEx), nanstd(NEy), abs( nanstd(gamma) ), 
                                nanstd(NEx_alt), nanstd(NEy_alt), abs(nanstd(gamma_alt)), nanstd(LEp)]

        ######################################################
        #Average using Scot-size box around the max point only
        ######################################################
        rgn = (A[:,2] <= (xcoord[locmax]+size_av/2)) & (A[:,2] >= (xcoord[locmax]-size_av/2)) & (A[:,3] <= (ycoord[locmax]+size_av/2)) &  (A[:,3] >= (ycoord[locmax]-size_av/2))
        Acol = A[ rgn , :]
        
        colLEp, colNEx, colNEy, colNExy, colG, colx, coly  = ( n.zeros( len(Acol[:,0]) ) for _ in range(7) )
        colNEx_alt, colNEy_alt, colG_alt = (n.zeros( len(Acol[:,0]) ) for _ in range(3) )
        
        for m in range(len(Acol[:,0])):
            F = Acol[m,8:].reshape(2,2)
            U = sqrtm( (F.T).dot(F) ) 
            eigU = eigvalsh(U)
            LE = n.log(eigU)
            colLEp[m] = ( 2/3 * ( LE[0]**2 + LE[1]**2 + (-LE[0]-LE[1])**2 ) )**0.5
            NE_calc = U - n.eye(2)             #Technical strain calculation in the stretching coordinate system
            colNEx[m] = NE_calc[0,0]
            colNEy[m] = NE_calc[1,1]
            colNExy[m] = NE_calc[0,1]
            colG[m] = n.arctan(colNExy[m]/(1+colNEx[m])) + n.arctan(colNExy[m]/(1+colNEy[m]))
            colNEx_alt[m] = F[0,0]-1
            colNEy_alt[m] = F[1,1]-1
            colG_alt[m]=n.arctan(F[0,1]/F[1,1]);
            colx[m], coly[m] = Acol[m,[2,3]]
        OurSize_Scott = FS_th * thickness + .5*( (max(colx) - min(colx)) +( max(coly) - min(coly)) )
        export_Scott[k] = [STF[k,0], STF[k,1], nanmean(OurSize_Scott), STF[k,2], STF[k,3], nanmean(colNEx), nanmean(colNEy), abs(nanmean(colG)), 
                            nanmean(colNEx_alt), nanmean(colNEy_alt), abs(nanmean(colG_alt)), nanmean(colLEp)]                                
                                
        ######################################################
        # Profile of ep for different data_force_points
        ######################################################
        if k in profStg:
            Acol = A[ A[:,0] == aramXmaxlast, :]
            ystn = n.empty( (len(Acol[:,0]),2) )
            ystn[:,0] = Acol[:,6]/thickness
            for j in range(len(Acol[:,0])):
                F = Acol[j,8:].reshape(2,2)
                U = sqrtm( (F.T).dot(F) ) 
                eigU = eigvalsh(U)
                LE = n.log(eigU)
                ystn[j,1] = ( 2/3 * ( LE[0]**2 + LE[1]**2 + (-LE[0]-LE[1])**2 ) )**0.5
            profLEp[-(1+prof_count)] = ystn #Want the profiles put in ascendng order
            prof_count=prof_count+1        
        ######################################################
        # Strain contours
        ######################################################        
        if (k in profStg) and (makecontours):
            Acol = A[ abs(A[:,3])<=0.2, :]
            x=Acol[:,2] # [in]                                #Coordinates of all points
            y=Acol[:,3]
            LEp = n.empty(len(x))
            for m in range(len(x)):
                F = Acol[m,8:].reshape(2,2)
                U = sqrtm( (F.T).dot(F) ) 
                eigU = eigvalsh(U)
                LE = n.log(eigU)
                LEp[m] = ( 2/3 * ( LE[0]**2 + LE[1]**2 + (-LE[0]-LE[1])**2 ) )**0.5
            xspace=linspace(n.min(x),n.max(x),len( n.unique(Acol[:,0]) )/2 )
            yspace=linspace(n.min(y),n.max(y),len( n.unique(Acol[:,1]) )/2 )
            LEP = griddata( vstack((x,y)).T , LEp, (xspace[None,:],yspace[:,None]), method='linear')
            meanLEP = nanmean(LEP.flatten())
            stdLEP = nanstd(LEP.flatten())
            p.figure(facecolor='w',figsize=(16,12))
            p.contourf(xspace,yspace,LEP,linspace(meanLEP-0*stdLEP,meanLEP+3*stdLEP,256),extend='both')
            #p.axis('equal')
            p.xlabel('Undeformed X (in)')
            p.ylabel('Undeformed Y (in)')
            p.axis([n.min(x), n.max(x), n.min(y), n.max(y)])
            if k == last:
                z=1
                for m in range(len(MaxTen[:,0])):
                    p.plot(MaxTen[m,3],MaxTen[m,4],marker='${}$'.format(z),ms=12,mfc='w',mec='w')
                    z+=1
            else:
                p.plot(xcoord,ycoord,'o',mec='w',mfc='w',ms=6)
            p.plot(box_x,box_y,'w',linewidth=2.5)
            p.xlabel('Undeformed X (in)')
            p.ylabel('Undeformed Y (in) (axial)')
            if k == LL:
                p.title('$\mathsf{{e^{{p}}_{{e}}}}$: Station {} (Limit Load)'.format(n.where(profStg==k)[0][0]+1))
            else:
                p.title('$\mathsf{{e^{{p}}_{{e}}}}$: Station {}'.format(n.where(profStg==k)[0][0]+1))
            cbar = p.colorbar()
            cbar.set_label(r'     $\mathsf{e_{e}^{p}}$',rotation=0,fontsize=25)
            if not os.path.exists('./Contours'):
                os.mkdir('./Contours')
            p.savefig('./Contours/Stg{:.0f}'.format(k))
            p.close()
        ###############################################
        # Calculate delta and phi
        ###############################################

    lower = A[ (A[:,3] >= rdlim[0]) & (A[:,3] <= rdlim[1]) & (n.abs(A[:,2]) <= 0.3), :]
    upper = A[ (A[:,3] >= rdlim[2]) & (A[:,3] <= rdlim[3]) & (n.abs(A[:,2]) <= 0.3), :]
    L[k] = n.mean(upper[:,6]) - n.mean(lower[:,6])
    up_th[k] = n.mean( n.arctan2(upper[:,7], upper[:,5]) - n.arctan2(upper[:,4], upper[:,2]) )*180/n.pi
    lo_th[k] = n.mean( n.arctan2(lower[:,7], lower[:,5]) - n.arctan2(lower[:,4], lower[:,2]) )*180/n.pi


#########
# End of Loop
#########

# delta/L and phi
dL = (L - L[0]) / L[0]
dth = -(lo_th - up_th)

# Save data files!
#########
#Max.dat
headerline='[0]Stage [1]Time [2]AxSts [3]ShSts [4]NEx [5]NEy [6]Gamma [7]F11-1 [8]F22-1 [9]atan(F12/F22) [10]epeq [11]AramX [12]AramY'
n.savetxt('max.dat', X=export_max, fmt='%.0f, %.2f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.0f, %.0f',header=headerline)
#mean.dat
headerline='[0]Stage [1]Time [2]NumPtsPassed [3]AxSts [4]ShSts [5]NEx [6]NEy [7]Gamma [8]F11-1 [9]F22-1 [10]atan(F12/F22) [11]epeq'
n.savetxt('mean.dat', X=export_mean, fmt='%.0f, %.2f, %.0f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f',header=headerline)
#std.dat
headerline='[0]Stage [1]Time [2]NumPtsPassed [3]AxSts [4]ShSts [5]NEx [6]NEy [7]Gamma [8]F11-1 [9]F22-1 [10]atan(F12/F22) [11]epeq'
n.savetxt('std.dat', X=export_stdv, fmt='%.0f, %.2f, %.0f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f',header=headerline)
#like_Scott
headerline='[0]Stage [1]Time [2]SizeAveragingZone(in) [3]AxSts [4]ShSts [5]NEx [6]NEy [7]Gamma [8]F11-1 [9]F22-1 [10]atan(F12/F22) [11]epeq'
n.savetxt('like_Scott.dat', X=export_Scott, fmt='%.0f, %.2f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f',header=headerline)
#disp-rot
headerline='[0]Stage [1]Time [2]AxSts [3]ShSts [4]Delta/L [5]Phi. Lg = {:.6f} inch'.format(L[0])
n.savetxt('disp-rot.dat', X=hstack( (STF[:,[0,1,2,3]],dL[:,None],dth[:,None]) ), fmt='%.0f, %.2f, %.6f, %.6f, %.6f, %.6f',header=headerline)
# Save the Max10
headerline = '[0]Epeq [1]AramXIndex [2]AramYIndex [3]UndefXCoord [4]UndefYCoord'
n.savetxt('Max10.dat',X=MaxTen,fmt='%.6f, %.0f, %.0f, %.6f, %.6f',header=headerline)

# Save profiles (takes some work)
maxprolen = 0
for i in range(len(profStg)):
    prolen = len(profLEp[i][:,0])
    if prolen > maxprolen:
        maxprolen = prolen
allprofs = n.zeros( (1+maxprolen,len(profStg)*2) )
allprofs[allprofs == 0] = n.nan
for i in range(len(profStg)):
    allprofs[0,2*i] = profStg[i]
    allprofs[1:len(profLEp[i][:,0])+1,[2*i,2*i+1]] = profLEp[i]
# Save the Max10
headerline = 'First Row: Stage number. 2nd Row Begin y/t vs epeq'
n.savetxt('StrainProfiles.dat',X=allprofs,fmt='%.6f',delimiter=',') 

os.system('python F:\TT_New\AA_PyFiles\TT2_DIC_Analysis_Plots.py {} {} {} {} {}'.format(int(expt), alpha, int(FS), int(SS), savepath))
