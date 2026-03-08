import numpy as np
from astropy.io import fits,ascii
from astropy.table import Table,vstack
from dlnpyutils import utils as dln
import matplotlib.pyplot as plt
from dlnpyutils import utils as dln,plotting as pl
from scipy.stats import binned_statistic
from scipy.optimize import curve_fit

from scipy.special import gamma
def c4(n):
    return np.sqrt(2/(n-1)) * gamma(n/2) / gamma((n-1)/2)
corr = np.array([c4(n) for n in np.arange(20)])
corr[:2] = 1.0

def precision():


    elem = ['c','n','o','na','mg','al','si','s','k','ca','ti','v','cr','mn','co','ni','ce']
    tags = [e+'_fe' for e in elem]
    tags = ['teff','logg','fe_h']+tags

    tab=Table.read('allStarLiter-dr17-synspec_rev1.fits.gz')
    for c in tab.colnames:tab[c].name=c.lower()
    for t in tags: tab[t+'_sigma']=np.nan

    index=dln.create_index(tab['apogee_id'])
    gd, = np.where(index['num']>1)
    out = tab[index['lo'][gd]]
    for i in range(len(gd)):
        ind = index['index'][index['lo'][gd[i]]:index['hi'][gd[i]]+1]
        nind = len(ind)
        print(i+1,nind,tab['apogee_id'][ind[0]])
        tab1 = tab[ind]
        for j in range(len(tags)):
            vals = tab1[tags[j]]
            gg, = np.where(np.isfinite(vals) & (np.abs(vals) < 1e6))
            if len(gg)<2: continue
            vals = vals[gg]
            mn = np.mean(vals)
            std = np.std(vals)
            sig = np.sqrt(np.sum((vals-mn)**2)/(nind-1))
            err = sig/corr[len(gg)]
            tab1[tags[j]+'_sigma'] = err
        out[i] = tab1[0]
        
    out.write('dr17_precision_repeats.fits',overwrite=True)

def abunderrfunc(x,*pars):
    return np.sqrt(pars[0]**2+(pars[1]/x)**2)

def abunderrfunc2(x,*pars):
    return np.sqrt(pars[0]**2+(pars[1]/x))

def plots():
    tab = Table.read('dr17_precision_repeats.fits')

    elem = ['c','n','o','na','mg','al','si','s','k','ca','ti','v','cr','mn','co','ni','ce']
    tags = [e+'_fe' for e in elem]
    tags = ['teff','logg','fe_h']+tags

    gd, = np.where((tab['teff']<5000) & (tab['logg']<2.0))
    tab2 = tab[gd]

    fig = plt.figure(1,figsize=(8.5,5))
    
    snrbins = np.arange(20,500,10)
    for i in range(len(tags)):
        print(i+1,tags[i])
        name = tags[i].upper()
        if len(name)>1:
            name = tags[i][:1].upper()+tags[i][1:]
        if name[-3:]=='_fe': name=name[:-3]
        if name=='Logg': name='logg'
        vals = tab2[tags[i]+'_sigma'].data.data
        gg, = np.where(np.isfinite(vals))
        yr = [2e-3,1]
        if tags[i]=='teff': yr=[0.1,100]
        o=pl.scatter(tab2['snr'][gg],vals[gg],tab2['logg'][gg],size=5,ylog=True,
                     xr=[20,200],yr=yr,xtitle='S/N',ytitle=tags[i]+' sigma',title=name+' Uncertainty',charsize=18)
        res,xedge,_=binned_statistic(tab2['snr'][gg],vals[gg],bins=snrbins,statistic=np.nanmedian)
        plt.plot(xedge[:-1],res,c='r',label='Median binned')

        initpars = [0.01,0.02]
        pars1,pcov1 = curve_fit(abunderrfunc2,tab2['snr'][gg],vals[gg],p0=initpars,maxfev=5000)
        gdres, = np.where(np.isfinite(res))
        pars2,pcov2 = curve_fit(abunderrfunc2,xedge[gdres],res[gdres],p0=pars1,maxfev=5000)
        print(pars2)
        model1 = abunderrfunc2(snrbins,*pars1)
        model2 = abunderrfunc2(snrbins,*pars2)
        #plt.plot(snrbins,model1,c='purple')
        plt.plot(snrbins,model2,c='blue',label='Model fit')
        err50 = abunderrfunc2(50,*pars2)
        err70 = abunderrfunc2(70,*pars2)
        if tags[i]=='teff':
            plt.fill([10,120,120,10],[0.1,0.1,10,10],color='white')
            plt.annotate('S/N=50 $\sigma$ = {:4.1f} K'.format(err50),xy=[0.03,0.90],
                         xycoords='axes fraction',ha='left',fontsize=18)
            plt.annotate('S/N=70 $\sigma$ = {:4.1f} K'.format(err70),xy=[0.03,0.81],
                         xycoords='axes fraction',ha='left',fontsize=18)
        else:
            plt.fill([10,120,120,10],[0.3,0.3,0.9,0.9],color='white')
            plt.annotate('S/N=50 $\sigma$ = {:5.3f} dex'.format(err50),xy=[0.03,0.90],
                         xycoords='axes fraction',ha='left',fontsize=18)
            plt.annotate('S/N=70 $\sigma$ = {:5.3f} dex'.format(err70),xy=[0.03,0.81],
                         xycoords='axes fraction',ha='left',fontsize=18)
        plt.tick_params(axis='both', which='major', length=10)   # longer major ticks
        plt.tick_params(axis='both', which='minor', length=5)    # optional minor ticks
        plt.legend(fontsize=15,loc='upper right')
        plt.savefig('plots/dr17_precision_repeats_'+tags[i]+'.png',bbox_inches='tight')

def mkhtml():
        

    # Make HTML
    green = '#32CD32'
    lgreen = '#CCFF00'
    orange = '#FFBF00'
    yellow = '#FFFF00'
    red = '#FF0000'
    white = '#FFFFFF'
    ascolor = white

    elem = ['c','n','o','na','mg','al','si','s','k','ca','ti','v','cr','mn','co','ni','ce']
    tags = [e+'_fe' for e in elem]
    tags = tags+['teff','logg','fe_h']
        
    lines = []
    lines += ['<html>']
    lines += ['<head>']
    lines += ['<title>']
    lines += ['</title>']
    lines += ['</head>']
    lines += ['<body>']
    lines += ['<h1>APOGEE Abundance Uncertainties from Repeat Observations</h1>']
    lines += ['logg<2.0, Teff<5000, DR17 6128 stars<br>']
    lines += ['Click on figure to enlarge']
    lines += ['<table border=1>']
    lines += ['<tr><th>Number</th>']
    lines += ['<th>Element</th><th>DR17 Abund</th><th>Empirical Scatter</th>']
    lines += ['<th>Uncertainty from repeats</th><tr>']
    for i in range(len(tags)):
        name = tags[i].upper()
        if len(name)>1:
            name = tags[i][:1].upper()+tags[i][1:]
        if name[-3:]=='_fe': name=name[:-3]
        if name=='Logg': name='logg'
        print(i+1,tags[i],name)
        lines += ['<tr>']
        lines += ['<td style="text-align:center">{:}</td>'.format(i+1)]
        lines += ['<td style="text-align:center">{:}</td>'.format(name)]
        width = 500
        lines += ['<td><a href="plots/dr17_{:s}_feh_hist2d.png" target="_blank"><img src="plots/dr17_{:s}_feh_hist2d.png" width={:d}></a></td>'.format(name.lower(),name.lower(),width)]
        lines += ['<td><a href="plots/dr17_{:s}_scatter_snr.png" target="_blank"><img src="plots/dr17_{:s}_scatter_snr.png" width={:d}></a></td>'.format(name.lower(),name.lower(),width)]
        lines += ['<td><a href="plots/dr17_precision_repeats_{:s}.png" target="_blank"><img src="plots/dr17_precision_repeats_{:s}.png" width={:d}></a></td>'.format(tags[i],tags[i],width)]
        lines += ['</tr>']
    lines += ['</table>']
    lines += ['</body>']
    lines += ['</html>']
    dln.writelines('dr17_precision_repeats.html',lines)

        
def finalresults():
    """ Make final results plot """

    tab = ascii.read('precision_results.txt',format='commented_header')

    fig = plt.figure(1,figsize=(15,5))

    plt.rcParams['font.family'] = 'Arial'
    
    for i in range(len(tab)):
        name = tab['Element'][i]
        group = tab['group'][i]
        if i==0:
            label1 = 'Scatter'
            label2 = 'Pipeline'
            label3 = 'Repeats'
        else:
            label1 = None
            label2 = None
            label3 = None
        plt.scatter([i+1-0.25],tab['scatter_snr50'][i],s=10,c='k',label=label1)
        plt.scatter([i+1+0.25],tab['scatter_snr70'][i],s=20,c='k')
        plt.plot([i+1-0.25,i+1+0.25],[tab['scatter_snr50'][i],tab['scatter_snr70'][i]],color='k')
        plt.scatter([i+1-0.25],tab['pipeline_snr50'][i],s=10,c='r',label=label2)
        plt.scatter([i+1+0.25],tab['pipeline_snr70'][i],s=20,c='r')
        plt.plot([i+1-0.25,i+1+0.25],[tab['pipeline_snr50'][i],tab['pipeline_snr70'][i]],color='r')
        plt.scatter([i+1-0.25],tab['repeat_snr50'][i],s=10,c='green',label=label3)
        plt.scatter([i+1+0.25],tab['repeat_snr70'][i],s=20,c='green')
        plt.plot([i+1-0.25,i+1+0.25],[tab['repeat_snr50'][i],tab['repeat_snr70'][i]],color='green')
        plt.annotate(name,xy=[i+0.5,-0.05],xycoords='data',ha='center',fontsize=18)
        if i in [1,7,11,15,16]:
            plt.axvline(i+1.5,c='k')
        else:
            plt.axvline(i+1.5,c='k',linestyle='dotted',linewidth=1)
    plt.axvline(0.5,c='k')
    plt.axhline(0.1,c='gray',linewidth=0.3)
    plt.annotate('0.10 dex',xy=[1.01,0.49],xycoords='axes fraction',ha='left',fontsize=13,clip_on=False)
    plt.axhline(0.05,c='gray',linewidth=0.3)
    plt.annotate('0.05 dex',xy=[1.01,0.34],xycoords='axes fraction',ha='left',fontsize=13,clip_on=False)
    plt.axhline(0.02,c='gray',linewidth=0.3)
    plt.annotate('0.02 dex',xy=[1.01,0.14],xycoords='axes fraction',ha='left',fontsize=13,clip_on=False)
    plt.annotate('C+N',xy=[0.06,1.03],xycoords='axes fraction',ha='center',fontsize=16)
    plt.annotate(r'$\alpha$-elements',xy=[0.29,1.03],xycoords='axes fraction',ha='center',fontsize=16)
    plt.annotate('Odd-Z',xy=[0.59,1.03],xycoords='axes fraction',ha='center',fontsize=16)
    plt.annotate('Fe-peak',xy=[0.83,1.03],xycoords='axes fraction',ha='center',fontsize=16)
    plt.annotate('n-capture',xy=[0.97,1.03],xycoords='axes fraction',ha='center',fontsize=16)
    plt.xticks(np.arange(len(tab))+1,tab['Element'])
    ax = fig.axes[0]
    labels = ax.get_xticklabels()
    for i in range(len(labels)):
        if tab['priority'][i]==1:
            labels[i].set_fontweight('bold')
            labels[i].set_fontsize(20)
        elif tab['priority'][i]==2:
            labels[i].set_fontweight('medium')
            labels[i].set_fontsize(17)
        elif tab['priority'][i]==3:
            labels[i].set_fontweight('light')
            labels[i].set_fontsize(12)
    plt.xlim(-0.5,17.5)
    plt.ylim(1e-2,1)
    plt.yscale('log')
    plt.ylabel('Scatter/Precision (dex)')
    plt.legend(fontsize=12,framealpha=1)
    plt.scatter([0.75-1],[0.25],s=10,color='gray')
    plt.scatter([1.25-1],[0.20],s=20,color='gray')
    plt.plot([0.75-1,1.25-1],[0.25,0.20],color='gray')
    plt.annotate('S/N=50',xy=[0.9-1,0.27],xycoords='data',ha='center',fontsize=8)
    plt.annotate('S/N=70',xy=[1.15-1,0.17],xycoords='data',ha='center',fontsize=8)
    plt.savefig('plots/precision_results.png',bbox_inches='tight')
