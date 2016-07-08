'''
Really designed for plotting all experiments.  
The color pallette is different, and radial contraction is split into two figures:
One at station 2 and one at LL
'''


import numpy as n
L = n.genfromtxt
abs=n.abs
import os
from sys import argv
from pandas import read_excel
import figfun as ff
import os
import matplotlib.pyplot as p
p.close('all')
from cycler import cycler
 
colors = n.array([["#771155", "#AA4488", "#CC99BB"],
             ["#114477", "#4477AA", "#77AADD"],  
             ["#117777", "#44AAAA", "#77CCCC"],
             ["#117744", "#44AA77", "#88CCAA"],
             ["#777711", "#AAAA44", "#DDDD77"],
             ["#774411", "#AA7744", "#DDAA77"],
             ["#771122", "#AA4455", "#DD7788"]])

#colors = colors[::2].ravel()
#colors = n.hstack((colors,colors))
#colors = colors[[1,2,3,5,6,0]][:,[0,2]].ravel()
colors = colors[[2,3,-1,0]].ravel() 
#skip = 5           
#colors = n.hstack((colors.ravel()[::skip],colors.ravel()[::skip],colors.ravel()[::skip],colors.ravel()[::skip],colors.ravel()[::skip]))
 
expts = n.array([35,20,24,17,34,32,8,9,30,16,15,18])
alpha = ...
FS, SS = 19,6

savefigs = False

key = read_excel('../TT-Summary.xlsx',sheetname='Summary',header=None,index_col=None,skiprows=1).values

#expt, alpha, tube no, thickness
if (len(expts) >= 1) and (alpha == ...) :
    expinfo = key[ n.in1d(key[:,0], expts), :][:,[0,1,4,5]]
    savepath = '../ComparisonFigs/PaperSet'.format(alpha)
elif (type(alpha) in [int,float]) and (len(expts)==0):
    expinfo = key[ key[:,1]==alpha,: ][:,[0,1,4,5]]
    savepath = '../ComparisonFigs/Alpha-{}'.format(alpha)
    if n.isnan(alpha):
        expinfo = key[ n.isnan(key[:,1]),: ][:,[0,1,4,5]]
        savepath = '../ComparisonFigs/Alpha-Inf'.format(alpha)
else:
    raise ValueError('expts must be empty array OR alpha must be ellipsis')
    

if (len(expinfo.shape)==2) and (expinfo.shape[0]==1):
    expinfo = expinfo.ravel()
    expts = n.array([expinfo[0]]).astype(int)
    limloads = n.empty( (len(expts),4) ) #Force, torque, disp, rot
    stat2 = n.empty_like(limloads)
    expinfo = expinfo[None,:]
else:
    # A couple of bad expts we want to exclude
    expinfo = expinfo[ ~n.in1d(expinfo[:,0],[10,26]), :]
    order = n.flipud(n.argsort(expinfo[:,1]))
    expinfo = expinfo[order,:]
    expts = expinfo[:,0].astype(int)
    limloads = n.empty( (len(expts),4) ) #Force, torque, disp, rot
    stat2 = n.empty_like(limloads)
    
    
if (savefigs == True) and not (os.path.exists(savepath)):
    os.mkdir(savepath)

for G in range(len(expts)):
    
    relpath  = '../TT2-{}_FS{}SS{}'.format(expts[G],FS,SS)
            
    #########
    #Max.dat
    #   '[0]Stage [1]Time [2]AxSts [3]ShSts [4]NEx [5]NEy [6]Gamma [7]F11-1 [8]F22-1 [9]atan(F12/F22) [10]epeq [11]AramX [12]AramY'
    #mean.dat
    #   [0]Stage [1]Time [2]NumPtsPassed [3]AxSts [4]ShSts [5]NEx [6]NEy [7]Gamma [8]F11-1 [9]F22-1 [10]atan(F12/F22) [11]epeq'
    #like_Scott
    #   [0]Stage [1]Time [2]SizeAveragingZone(in) [3]AxSts [4]ShSts [5]NEx [6]NEy [7]Gamma [8]F11-1 [9]F22-1 [10]atan(F12/F22) [11]epeq'
    # profStgs
    #   'Stages at which profiles were generated'
    # profUr
    #   'First Row Stage number.  Second row begin data.\n[0]Ycoord [1:] Stage Ur/Ro'
    #MaxPt.dat
    #   [0]Stage [1]Time [2]AxSts [3]ShSts [4]NEx [5]NEy [6]Gamma [7]F11-1 [8]F22-1 [9]atan(F12/F22) [10]epeq'

    STF = L('{}/STF.dat'.format(relpath),delimiter=',')
    profStg = L('{}/prof_stages.dat'.format(relpath),delimiter=',').astype(int)
    dmax = L('{}/max.dat'.format(relpath),delimiter=',')
    dmaxPt = L('{}/MaxPt.dat'.format(relpath),delimiter=',')
    dmean = L('{}/mean.dat'.format(relpath),delimiter=',')
    dscot = L('{}/like_Scott.dat'.format(relpath),delimiter=',')
    DR = L('{}/disp-rot.dat'.format(relpath),delimiter=',')
    #'[0]Stage [1]Time [2]AxSts [3]ShSts [4]Delta/L [5]Phi. Lg = {:.6f} inch'
    limloads[G] = DR[int(profStg[2]),2:]
    stat2[G] = DR[int(profStg[1]),2:]
    profLEp = L('{}/StrainProfiles.dat'.format(relpath),delimiter=',')[1:]
    profUr = L('{}/RadialContraction.dat'.format(relpath),delimiter=',')[1:]
    profUr[:,1:] -= n.nanmean( profUr[profUr[:,0]>=.35,:][:,1:], axis=0 )
    
    ##################################################
    # Figure 1 - AxSts-Delta and ShearSts-Rot
    ##################################################
    if G == 0:
        p.style.use('mysty-sub')
        p.rcParams['axes.prop_cycle'] = cycler('color',colors)
        fig1 =  p.figure(1,facecolor='w',figsize=(8,12))
        ax11 = fig1.add_subplot(2,1,1)
        ax12 = fig1.add_subplot(2,1,2)
    
    masterline, = ax11.plot(DR[:,4],DR[:,2],label = '{:.0f}; {}; {:.0f}'.format(expts[G],expinfo[G,1],expinfo[G,3]))
    mastercolor = masterline.get_color()
    masterlabel = masterline.get_label()
    ax12.plot(DR[:,5],DR[:,3],label = '{:.0f}; {}; {:.0f}'.format(expts[G],expinfo[G,1],expinfo[G,3]), color=mastercolor)
    
    if G == len(expts)-1:
        for J in range(len(expts)):
            m11 = ax11.plot(limloads[J,2],limloads[J,0],'^',mec='r',mfc='r',ms=6)[0]
            m21 = ax12.plot(limloads[J,3],limloads[J,1],'^',mec='r',mfc='r',ms=6)[0]
            m12 = ax11.plot(stat2[J,2],stat2[J,0],'o',mec='r',mfc='r',ms=6)[0]
            m22 = ax12.plot(stat2[J,3],stat2[J,1],'o',mec='r',mfc='r',ms=6)[0]
        
        ax11.set_title('Nominal Response',size=18)
        ax11.set_xlabel('$\\delta/\\mathsf{L}$')
        ax11.set_ylabel('$\\Sigma$\n$(\\mathsf{ksi})$')
        ax11.set_ylim([0,1.2*n.max(limloads[:,0])])
        ax11.set_xlim(left=0)
        l1 = ax11.legend([m12,m11],["Station 2", "LL"],loc='upper right',numpoints=1,fontsize=10,frameon=False)
        p.setp(l1.get_texts(),color='r')
        l2 = ax11.legend(loc='lower right',fontsize=10,title='Expt; $\\alpha$; Tube')
        p.setp(l2.get_title(),fontsize=10)
        ax11.add_artist(l1)
        ax12.set_xlabel('$\\phi^\\circ$')
        ax12.set_ylabel('$\\mathcal{T}$\n$(\\mathsf{ksi})$')
        ax12.set_ylim([0,1.2*n.max(limloads[:,1])])
        ax12.set_xlim(left=0)
        l1 = ax12.legend([m22,m21],["Station 2", "LL"],loc='upper right',numpoints=1,fontsize=10,frameon=False)
        p.setp(l1.get_texts(),color='r')
        l2 = ax12.legend(loc='lower right',fontsize=10,title='Expt; $\\alpha$; Tube')
        p.setp(l2.get_title(),fontsize=10)
        ax12.add_artist(l1)
        
        p.sca(ax11)
        ff.myax(fig1,ff.ksi2Mpa,'$\\Sigma$\n$(\\mathsf{MPa})$')
        p.sca(ax12)
        ff.myax(fig1,ff.ksi2Mpa,'$\\mathcal{T}$\n$(\\mathsf{MPa})$')

        if savefigs:
            #fig1.savefig('{}/1 - Sts-Delta-Rot.png'.format(savepath),dpi=125)
            #fig1.savefig('{}/1 - Sts-Delta-Rot.pdf'.format(savepath),dpi=125)
            p.close()

    ##################################################
    # Figure 2 - Epsilon-Gamma
    ##################################################
    if G == 0:
        p.style.use('mysty-sub')
        p.rcParams['axes.prop_cycle'] = cycler('color',colors)
        fig2 = p.figure(2,facecolor='w',figsize=(8,12) )
        ax21 = fig2.add_subplot(2,1,1)
        ax22 = fig2.add_subplot(2,1,2)
    
    ax21.plot(abs(dmean[:,7]),dmean[:,6],'o',ms=4,mfc=mastercolor,mec='none',label=masterlabel)
    ax21.plot(abs(dmax[-1,6]),dmax[-1,5],'s',ms=8,mfc=mastercolor,mec='none')
    ax22.plot(abs(dmean[:,10]),dmean[:,9],'o',ms=4,mfc=mastercolor,mec='none',label=masterlabel)
    ax22.plot(abs(dmax[-1,9]),dmax[-1,8],'s',mfc=mastercolor,mec='none',ms=8)
    
    if G == len(expts)-1:
        ax21.set_title('Mean strain path; Aramis User manual Definitions',fontsize=14)
        ax21.set_xlabel('$\gamma$')
        ax21.set_ylabel('$\\epsilon_y$')
        ax21.set_ylim(bottom=0)
        ax21.set_xlim(left=0)
        l2 = ax21.legend(loc='center left',bbox_to_anchor=(1.01,.5),fontsize=10,numpoints=1,handletextpad=.1,title='Expt; $\\alpha$; Tube')
        p.setp(l2.get_title(),fontsize=10)
        #ax21.legend(loc='lower left')
        ax22.set_title('Mean strain path; Haltom 2013 Definitions',fontsize=14)
        ax22.set_xlabel('$\\gamma = atan(\\mathsf{F}_{\\mathsf{01}}/\\mathsf{F}_{\\mathsf{11}}$)')
        ax22.set_ylabel('$\\epsilon_{\\mathsf{y}}$\n$\\mathsf{F}_{\\mathsf{11}}-\\mathsf{1}$')
        ax22.set_ylim(bottom=0)
        ax22.set_xlim(left=-.015)
        l2 = ax22.legend(loc='center left',bbox_to_anchor=(1.01,.5),fontsize=10,numpoints=1,handletextpad=.1,title='Expt; $\\mathregular{\\alpha}}$; Tube')
        p.setp(l2.get_title(),fontsize=10)
        #ax22.legend(loc='lower left')
        
        p.sca(ax21)
        ff.myax(fig2)
        p.sca(ax22)
        ff.myax(fig2)
        
        if savefigs:
            #fig2.savefig('{}/2 - StrainPath.png'.format(savepath),dpi=125,bbox_inches='tight')
            #fig2.savefig('{}/2 - StrainPath.pdf'.format(savepath),dpi=125,bbox_inches='tight')
            p.close()

    ##################################################
    # Figure 3 - Radial contraction at Station 2
    ##################################################
    if G == 0:
        p.style.use('mysty')
        p.rcParams['axes.prop_cycle'] = cycler('color',colors)
        fig3 = p.figure(3,facecolor='w',figsize=(12,6) )
        ax3 = fig3.add_axes([.12,.12,.8,.78])
    
    # Station 2
    ax3.plot(2*profUr[:,0]/0.62,profUr[:,1+1],'-',lw=1.5,color=mastercolor,label=masterlabel)

    if G == len(expts) - 1:
        
        ax3.set_title('Radial Contraction at Station 2 and LL',size=18)
        ax3.set_xlabel('$\\mathsf{2}\\mathsf{y}_{\\mathsf{o}}/\\mathsf{L}_{\\mathsf{g}}$')
        ax3.set_ylabel('$\\frac{\\mathsf{u}_{\\mathsf{r}}}{\\mathsf{R}_{\\mathsf{o}}}$')
        ax3.set_xlim([-1,1])
        ax3.grid(True,which='both',axis='y',linestyle='--',alpha=0.5)
        l3 = ax3.legend(loc='lower right',title='Expt; $\\alpha$; Tube',frameon=True)
        p.setp(l3.get_title(),fontsize=12)
        p.setp(l3.get_frame(),lw=0)
        ax3.yaxis.set_ticks(n.arange(-.024,0+.004,.004))
        ax3.set_ylim([-.024,.002])
        
        p.sca(ax3)
        ff.myax(fig3,TW=.002,HW=.3,HL=.05,OH=.2)
        
        if savefigs:
            fig3.savefig('{}/3 - RadialContraction.png'.format(savepath),dpi=125,bbox_inches='tight')
            fig3.savefig('{}/3 - RadialContraction.pdf'.format(savepath),dpi=125,bbox_inches='tight')
            pass
    ##################################################
    # Figure 5 - Radial contraction at LL
    ##################################################
    if G == 0:
        p.style.use('mysty')
        p.rcParams['axes.prop_cycle'] = cycler('color',colors)
        fig5 = p.figure(5,facecolor='w',figsize=(12,6) )
        ax5 = fig5.add_axes([.12,.12,.8,.78])
    
    # Station 2
    ax5.plot(2*profUr[:,0]/0.62,profUr[:,1+2],'-',lw=1.5,color=mastercolor,label=masterlabel)
   
    if G == len(expts) - 1:
        
        ax5.set_title('Radial Contraction at Station 2 and LL',size=18)
        ax5.set_xlabel('$\\mathsf{2}\\mathsf{y}_{\\mathsf{o}}/\\mathsf{L}_{\\mathsf{g}}$')
        ax5.set_ylabel('$\\frac{\\mathsf{u}_{\\mathsf{r}}}{\\mathsf{R}_{\\mathsf{o}}}$')
        ax5.set_xlim([-1,1])
        ax5.grid(True,which='both',axis='y',linestyle='--',alpha=0.5)
        l3 = ax5.legend(loc='lower right',title='Expt; $\\alpha$; Tube',frameon=True)
        p.setp(l3.get_title(),fontsize=12)
        p.setp(l3.get_frame(),lw=0)
        ax5.yaxis.set_ticks(n.arange(-.024,0+.004,.004))
        ax5.set_ylim([-.024,.002])
        
        p.sca(ax5)
        ff.myax(fig5,TW=.002,HW=.3,HL=.05,OH=.2)
        
        if savefigs:
            fig5.savefig('{}/3 - RadialContraction_LL.png'.format(savepath),dpi=125,bbox_inches='tight')
            fig5.savefig('{}/3 - RadialContraction_LL.pdf'.format(savepath),dpi=125,bbox_inches='tight')
            pass
            
    ##################################################
    # Figure 4 - Strain profile
    ##################################################
    if G == 0:
        p.style.use('mysty')
        p.rcParams['axes.prop_cycle'] = cycler('color',colors)
        fig4 = p.figure(4,facecolor='w',figsize=(12,6) )
        ax4 = fig4.add_axes([.12,.12,.8,.78])
    
    ax4.plot(profLEp[:,-3]/expinfo[G,2],profLEp[:,-1],color=mastercolor,label=masterlabel)
    if G == len(expts) - 1:
        
        ax4.set_title('$\\mathsf{e}^{\\mathsf{p}}_{\\mathsf{e}}$ at Failure',size=18)
        ax4.set_xlabel('y$_{\\mathsf{o}}$/t$_{\\mathsf{o}}$')
        ax4.set_ylabel('$\\mathsf{e}^{\\mathsf{p}}_{\\mathsf{e}}$')
        ax4.set_xlim([-8,8])
        l2 = ax4.legend(loc='upper right',title='Expt; $\\alpha$; Tube')
        p.setp(l2.get_title(),fontsize=12)
        p.sca(ax4)
        ff.myax(fig4,TW=.0025,HW=.3,HL=.05,OH=.2)
        
        if savefigs:
            fig4.savefig('{}/4 - Strain Profile.png'.format(savepath),dpi=125,bbox_inches='tight')
            fig4.savefig('{}/4 - Strain Profile.pdf'.format(savepath),dpi=125,bbox_inches='tight')
            p.close()   
            pass
            
    ##################################################
    # Figure 6 - Equiv Stn at LL 
    ##################################################
    if G == 0:
        p.style.use('mysty')
        p.rcParams['axes.prop_cycle'] = cycler('color',colors)
        fig6 = p.figure(6)
        ax6 = fig6.gca()
    
    sig, tau = dmax[profStg[2]][2:4]
    triax = (sig/2)/n.sqrt(.75*(sig**2+4*tau**2))
    LL, = ax6.plot(triax, dmax[profStg[2],10],'gs',mec='g',label='LL')
    sig, tau = dmax[profStg[-1]][2:4]
    triax = (sig/2)/n.sqrt(.75*(sig**2+4*tau**2))
    END, = ax6.plot(triax, dmax[profStg[-1],10],'rs',mec='none',label='Fail.')
    tx = ax6.text(triax,dmax[profStg[-1],10],'{}\n'.format(expinfo[G,1]),color=mastercolor,ha='center',size=10)
    if expts[G] == 8:
        # This experiment is a mess!
        LL.remove()
        END.remove()
        tx.remove()
    if G == len(expts) - 1:
        ax6.axis(xmin=0)
        ax6.set_title('$\\mathsf{e}^{\\mathsf{p}}_{\\mathsf{e}}$ at Limit Load',size=18)
        ax6.set_xlabel('$\\sigma_{\\mathsf{m}}/\\sigma_{\\mathsf{e}}$')
        ax6.set_ylabel('$\\mathsf{e}^{\\mathsf{p}}_{\\mathsf{e}}$')
        leg6 = ax6.legend([LL,END],['LL','Fail'],loc='upper right')
        [tx.set_color(leg6.get_lines()[k].get_color()) for k,tx in enumerate(leg6.get_texts())]
        ff.myax(fig6)
        
        if savefigs:
            fig6.savefig('{}/6 - eeq at LL.png'.format(savepath),dpi=125,bbox_inches='tight')
            fig6.savefig('{}/6 - eeq at LL.pdf'.format(savepath),dpi=125,bbox_inches='tight')
            p.close()   
            pass
            
if not savefigs:
    p.show('all')        
else:
    p.close('all')
