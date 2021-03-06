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

# Specify expts or alpha!

try:
    argv = argv[1:]
    if len(argv) == 1:
        # we only have alpha!
        alpha = float(argv[0])
        expts = n.array([])
    elif len(argv)>1:
        alpha = ...
        expts = n.array(argv).astype(int)
    else:
        raise
except:
    expts = n.array([])
    alpha = n.nan

FS, SS = 19, 6
savefigs = False

key = read_excel('../TT-Summary.xlsx',sheetname='Summary',header=None,index_col=None,skiprows=1).values

#expt, alpha, tube no, thickness, true alpha, ecc
if (len(expts) >= 1) and (alpha == ...) :
    expinfo = key[ n.in1d(key[:,0], expts), :][:,[0,1,4,5,2,6]]
    savepath = '../ComparisonFigs/PaperSet'.format(alpha)
elif (type(alpha) in [int,float]) and (len(expts)==0):
    expinfo = key[ key[:,1]==alpha,: ][:,[0,1,4,5,2,6]]
    savepath = '../ComparisonFigs/Alpha-{}'.format(alpha)
    print(savepath)
    if n.isnan(alpha):
        expinfo = key[ n.isnan(key[:,1]),: ][:,[0,1,4,5,2,6]]
        savepath = '../ComparisonFigs/Alpha-Inf'.format(alpha)
else:
    raise ValueError('expts must be empty array OR alpha must be ellipsis')
    

if (len(expinfo.shape)==2) and (expinfo.shape[0]==1):
    # Specifying Experiments
    expinfo = expinfo.ravel()
    expts = n.array([expinfo[0]]).astype(int)
    limloads = n.empty( (len(expts),4) ) #Force, torque, disp, rot
    stat2 = n.empty_like(limloads)
    expinfo = expinfo[None,:]
else:
    expinfo = expinfo[ ~n.in1d(expinfo[:,0],[10,26,36]), :]    # A couple of bad expts we want to exclude
    order = n.argsort(expinfo[:,1])
    expinfo = expinfo[order,:]
    expts = expinfo[:,0].astype(int)
    limloads = n.empty( (len(expts),4) ) #Force, torque, disp, rot
    stat2 = n.empty_like(limloads)
    savepath = '../ComparisonFigs/Alpha-1.5_NoPlugExpt'
    
legtitle = 'Expt; $\\alpha$; Tube; $\\Xi$'
    
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
    profStg = L('{}/prof_stages.dat'.format(relpath),delimiter=',')
    dmax = L('{}/max.dat'.format(relpath),delimiter=',')
    dmaxPt = L('{}/MaxPt.dat'.format(relpath),delimiter=',')
    dmean = L('{}/mean.dat'.format(relpath),delimiter=',')
    dscot = L('{}/like_Scott.dat'.format(relpath),delimiter=',')
    DR = L('{}/disp-rot.dat'.format(relpath),delimiter=',')
    #'[0]Stage [1]Time [2]AxSts [3]ShSts [4]Delta/L [5]Phi. Lg = {:.6f} inch'
    limloads[G] = DR[int(profStg[2]),2:]
    stat2[G] = DR[int(profStg[1]),2:]
    profLEp = L('{}/StrainProfiles.dat'.format(relpath),delimiter=',')
    profUr = L('{}/RadialContraction.dat'.format(relpath),delimiter=',')[1:]
    
    ##################################################
    # Figure 1 - AxSts-Delta and ShearSts-Rot
    ##################################################
    if G == 0:
        p.style.use('mysty-sub')
        fig1 =  p.figure(1,facecolor='w',figsize=(8,12))
        ax11 = fig1.add_subplot(2,1,1)
        ax12 = fig1.add_subplot(2,1,2)
    
    LINE, = ax11.plot(DR[:,4],DR[:,2],label = '{:.0f}; {:.2f}; {:.0f}; {:.1f}'.format(expts[G],expinfo[G,4],expinfo[G,3],expinfo[G,5]))
    ax12.plot(DR[:,5],DR[:,3],label=LINE.get_label())
    
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
        l2 = ax11.legend(loc='lower right',fontsize=10,title=legtitle)
        p.setp(l2.get_title(),fontsize=10)
        ax11.add_artist(l1)
        ax12.set_xlabel('$\\phi^\\circ$')
        ax12.set_ylabel('$\\mathcal{T}$\n$(\\mathsf{ksi})$')
        ax12.set_ylim([0,1.2*n.max(limloads[:,1])])
        ax12.set_xlim(left=0)
        l1 = ax12.legend([m22,m21],["Station 2", "LL"],loc='upper right',numpoints=1,fontsize=10,frameon=False)
        p.setp(l1.get_texts(),color='r')
        l2 = ax12.legend(loc='lower right',fontsize=10,title=legtitle)
        p.setp(l2.get_title(),fontsize=10)
        ax12.add_artist(l1)
        
        p.sca(ax11)
        ff.myax(fig1,ff.ksi2Mpa,'$\\Sigma$\n$(\\mathsf{MPa})$')
        p.sca(ax12)
        ff.myax(fig1,ff.ksi2Mpa,'$\\mathcal{T}$\n$(\\mathsf{MPa})$')

        if savefigs:
            fig1.savefig('{}/1 - Sts-Delta-Rot.png'.format(savepath),dpi=125)
            #p.close()

    ##################################################
    # Figure 2 - Epsilon-Gamma
    ##################################################
    if G == 0:
        p.style.use('mysty-sub')
        fig2 = p.figure(2,facecolor='w',figsize=(8,12) )
        ax21 = fig2.add_subplot(2,1,1)
        ax22 = fig2.add_subplot(2,1,2)
    
    ax21.plot(abs(dmean[:,7]),dmean[:,6],'o',ms=4,mfc=LINE.get_color(),mec=LINE.get_color(),label=LINE.get_label())
    ax21.plot(abs(dmax[-1,6]),dmax[-1,5],'s',ms=8,mfc=LINE.get_color(),mec=LINE.get_color())
    ax22.plot(abs(dmean[:,10]),dmean[:,9],'o',ms=4,mfc=LINE.get_color(),mec=LINE.get_color(),label=LINE.get_label())
    ax22.plot(abs(dmax[-1,9]),dmax[-1,8],'s',mfc=LINE.get_color(),mec=LINE.get_color(),ms=8)
    
    if G == len(expts)-1:
        ax21.set_title('Mean strain path; Aramis User manual Definitions',fontsize=14)
        ax21.set_xlabel('$\gamma$')
        ax21.set_ylabel('$\\epsilon_y$')
        ax21.set_ylim(bottom=0)
        ax21.set_xlim(left=0)
        l2 = ax21.legend(loc='center left',bbox_to_anchor=(1.01,.5),fontsize=10,numpoints=1,handletextpad=.1,title=legtitle)
        p.setp(l2.get_title(),fontsize=10)
        ax22.set_title('Mean strain path; Haltom 2013 Definitions',fontsize=14)
        ax22.set_xlabel('$\\gamma = atan(\\mathsf{F}_{\\mathsf{01}}/\\mathsf{F}_{\\mathsf{11}}$)')
        ax22.set_ylabel('$\\epsilon_{\\mathsf{y}}$\n$\\mathsf{F}_{\\mathsf{11}}-\\mathsf{1}$')
        ax22.set_ylim(bottom=0)
        ax22.set_xlim(left=0)
        l2 = ax22.legend(loc='center left',bbox_to_anchor=(1.01,.5),fontsize=10,numpoints=1,handletextpad=.1,title=legtitle)
        p.setp(l2.get_title(),fontsize=10)
        
        p.sca(ax21)
        ff.myax(fig2)
        p.sca(ax22)
        ff.myax(fig2)
        
        if savefigs:
            fig2.savefig('{}/2 - StrainPath.png'.format(savepath),dpi=125,bbox_inches='tight')
            #p.close()

    ##################################################
    # Figure 3 - Radial contraction at the LL and at stage prior to
    ##################################################
    if G == 0:
        p.style.use('mysty')
        fig3 = p.figure(3,facecolor='w',figsize=(12,6) )
        ax3 = fig3.add_axes([.12,.12,.8,.78])
        plt3 = []
    
    # Station 2
    mark1, = ax3.plot(2*profUr[:,0]/0.62,profUr[:,1+1],'--',lw=1.5,color=LINE.get_color())
    # Limit Load
    mark2, = ax3.plot(2*profUr[:,0]/0.62,profUr[:,2+1],lw=1.5,color=LINE.get_color(),label=LINE.get_label())
    
    if G == len(expts) - 1:
        
        ax3.set_title('Radial Contraction at Station 2 and LL',size=18)
        ax3.set_xlabel('$\\mathsf{2}\\mathsf{y}_{\\mathsf{o}}/\\mathsf{L}_{\\mathsf{g}}$')
        ax3.set_ylabel('$\\frac{\\mathsf{u}_{\\mathsf{r}}}{\\mathsf{R}_{\\mathsf{o}}}$')
        ax3.set_xlim([-1,1])
        ax3.grid(True,which='both',axis='y',linestyle='--',alpha=0.5)
        l2 = ax3.legend([mark1, mark2] ,['Station 2','Limit Load'],loc='lower left')
        p.setp(l2.get_lines(),color='k')
        l3 = ax3.legend(loc='lower right',title=legtitle,frameon=True)
        p.setp(l3.get_title(),fontsize=12)
        p.setp(l3.get_frame(),lw=0)
        ax3.add_artist(l2)
        yl = ax3.get_ylim()
        ax3.yaxis.set_ticks(n.arange(-.024,0+.004,.004))
        ax3.set_ylim(yl)
        
        p.sca(ax3)
        ff.myax(fig3,TW=.002,HW=.3,HL=.05,OH=.2)
        
        if savefigs:
            fig3.savefig('{}/3 - RadialContraction.png'.format(savepath),dpi=125,bbox_inches='tight')
            #p.close()   

    ##################################################
    # Figure 4 - Strain profile
    ##################################################
    if G == 0:
        p.style.use('mysty')
        fig4 = p.figure(4,facecolor='w',figsize=(12,6) )
        ax4 = fig4.add_axes([.12,.12,.8,.78])
    
    ax4.plot(profLEp[:,-3]/expinfo[G,2],profLEp[:,-1],color=LINE.get_color(),label=LINE.get_label())
    
    if G == len(expts) - 1:
        
        ax4.set_title('$\\mathsf{e}^{\\mathsf{p}}_{\\mathsf{e}}$ at Failure',size=18)
        ax4.set_xlabel('y$_{\\mathsf{o}}$/t$_{\\mathsf{o}}$')
        ax4.set_ylabel('$\\mathsf{e}^{\\mathsf{p}}_{\\mathsf{e}}$')
        ax4.set_xlim([-8,8])
        l2 = ax4.legend(loc='upper right',title=legtitle)
        p.setp(l2.get_title(),fontsize=12)
        
        p.sca(ax4)
        ff.myax(fig4,TW=.0025,HW=.3,HL=.05,OH=.2)
        
        if savefigs:
            fig4.savefig('{}/4 - Strain Profile.png'.format(savepath),dpi=125,bbox_inches='tight')
            #p.close()   
            
if not savefigs:
    p.show('all')        
else:
    p.close('all')
