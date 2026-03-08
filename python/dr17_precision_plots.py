import os
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from dlnpyutils import utils as dln,plotting as pl
from scipy.stats import binned_statistic


def precision_plots():
    """ APOGEE DR17 abundance precision plots. """
    
    #tab = Table.read('/Users/nidever/as5/hge/abundances/allStarLiter-dr17-synspec_rev1_matched.fits.gz')
    #tab = tab.filled()
    ##gd,=np.where( (np.abs(tab['glat'])<10) )
    ##gd,=np.where( (np.abs(tab['glat'])<10) & (tab['starflag']==0) & (tab['logg']<2.0))
    #gd, = np.where( (np.abs(tab['glat'])<10) & (tab['starflag']==0) & (tab['teff']<5000) & (tab['logg']<1.5))
    #tab = tab[gd]
    #tab.write('/Users/nidever/as5/hge/abundances/allStarLiter-dr17-synspec_rev1_matched_rgb.fits')
    tab = Table.read('/Users/nidever/as5/hge/abundances/allStarLiter-dr17-synspec_rev1_matched_rgb.fits')
    tab = tab.filled()

    # Low-alpha sequence
    lowalpha, = np.where((tab['o_fe'] < 0.2) & (tab['o_fe'] > 0.0))
    
    elem = ['c','n','o','na','mg','al','si','s','k','ca','ti','v','cr','mn','co','ni','ce']
    
    
    #plt.figure(1,figsize=(10,7))
    plt.figure(1,figsize=(8.5,5))
    
    xr = [-2.5,1.0]
    yr = [-0.5,0.75]
    fehbins = np.arange(-0.8,0.6,0.05)
    snrbins = np.arange(20,220,10)
    
    # Element plots
    for i in range(len(elem)):
        el = elem[i]
        name = el.upper()
        if len(el)==2:
            name = el[0].upper()+el[1]
        print(i+1,name)
            
        # Abundance vs. feh plot, hist2d
        abund1 = tab[el+'_fe']
        gd1, = np.where((np.abs(abund1) < 10) & np.isfinite(abund1))
        o=pl.hist2d(tab['fe_h'][gd1],abund1[gd1],log=True,xr=xr,yr=yr,xtitle='[Fe/H]',
                    ytitle='['+name+'/Fe]',title='APOGEE DR17 - '+name,charsize=20)

        res,xedge,_ = binned_statistic(tab['fe_h'][gd1],abund1[gd1],bins=fehbins,statistic=np.nanmedian)
        sm1 = dln.interp(xedge[:-1],res,tab['fe_h'][gd1])
        diff1 = abund1[gd1]-sm1
        sig1 = dln.mad(diff1)
        plt.plot(xedge[:-1],res,c='r')
        plt.annotate('Scatter = {:.3f} dex'.format(sig1),xy=[0.05,0.92],xycoords='axes fraction',ha='left',fontsize=18)
        plt.ylim(yr)
        plt.savefig('plots/dr17_'+el+'_feh_hist2d.png',bbox_inches='tight')

        # Abundance vs. feh plot, scatter color-coded by Teff
        abund1 = tab[el+'_fe']
        gd1, = np.where((np.abs(abund1) < 10) & np.isfinite(abund1))
        o=pl.scatter(tab['fe_h'][gd1],abund1[gd1],tab['teff'][gd1],size=5,xr=xr,yr=yr,xtitle='[Fe/H]',
                     ytitle='['+name+'/Fe]',title='APOGEE DR17 - '+name+' (Teff)',charsize=20)
        res,xedge,_ = binned_statistic(tab['fe_h'][gd1],abund1[gd1],bins=fehbins,statistic=np.nanmedian)
        sm1 = dln.interp(xedge[:-1],res,tab['fe_h'][gd1])
        diff1 = abund1[gd1]-sm1
        sig1 = dln.mad(diff1)
        plt.plot(xedge[:-1],res,c='r')
        plt.annotate('Scatter = {:.3f} dex'.format(sig1),xy=[0.05,0.90],xycoords='axes fraction',ha='left',fontsize=18)
        plt.ylim(yr)
        plt.savefig('plots/dr17_'+el+'_feh_teff.png',bbox_inches='tight')
        
        # Scatter vs. S/N plot

        # scatter in low-alpha sequence vs. S/N
        if el not in ['ti','mn','ce']:
            low, = np.where((tab['o_fe'] < 0.2) & (tab['o_fe'] > 0.0) & (np.abs(abund1) < 10) & np.isfinite(abund1))
            res,xedge,_ = binned_statistic(tab['fe_h'][low],abund1[low],bins=fehbins,statistic=np.nanmedian)
            sm1 = dln.interp(xedge[:-1],res,tab['fe_h'][low])
            diff1 = abund1[low]-sm1
            sig1 = dln.mad(diff1)
            res2,xedge2,_ = binned_statistic(tab['snr'][low],diff1,bins=snrbins,statistic=dln.mad)
            label = r'low-$\alpha$ Residual Scatter'
        else:
            res,xedge,_ = binned_statistic(tab['fe_h'][gd1],abund1[gd1],bins=fehbins,statistic=np.nanmedian)
            sm1 = dln.interp(xedge[:-1],res,tab['fe_h'][gd1])
            diff1 = abund1[gd1]-sm1
            sig1 = dln.mad(diff1)
            res2,xedge2,_ = binned_statistic(tab['snr'][gd1],diff1,bins=snrbins,statistic=dln.mad)
            label = r'Residual Scatter'
        pl.plot(snrbins[:-1],res2,c='r',label=label,xtitle='S/N',ytitle='Scatter',title='',charsize=18)
        # pipeline uncertainties
        res3,xedge3,_ = binned_statistic(tab['snr'][gd1],np.log10(tab[el+'_fe_err'][gd1]),bins=snrbins,statistic='median')
        plt.plot(snrbins[:-1],10**res3,c='k',label='pipeline uncertainties')

        # Values at S/N=50 and 70
        scat50 = res2[np.where(snrbins==50)][0]
        scat70 = res2[np.where(snrbins==70)][0]
        err50 = 10**res3[np.where(snrbins==50)][0]
        err70 = 10**res3[np.where(snrbins==70)][0]
        plt.annotate('S/N=50 $\sigma$ = {:5.3f}'.format(err50),xy=[0.03,0.90],xycoords='axes fraction',
                     ha='left',fontsize=18)
        plt.annotate('{:5.3f}'.format(scat50),color='r',xy=[0.39,0.90],xycoords='axes fraction',ha='left',fontsize=18)
        plt.annotate('S/N=70 $\sigma$ = {:5.3f}'.format(err70),xy=[0.03,0.81],xycoords='axes fraction',
                     ha='left',fontsize=18)
        plt.annotate('{:5.3f}'.format(scat70),color='r',xy=[0.39,0.81],xycoords='axes fraction',ha='left',fontsize=18)
        plt.ylim(2e-3,1)
        plt.yscale('log')
        plt.title(name+' Scatter')
        plt.tick_params(axis='both', which='major', length=10)   # longer major ticks
        plt.tick_params(axis='both', which='minor', length=5)    # optional minor ticks
        plt.legend(fontsize=15)
        plt.savefig('plots/dr17_'+el+'_scatter_snr.png',bbox_inches='tight')
        

        
    #import pdb; pdb.set_trace()


    # Make HTML
    green = '#32CD32'
    lgreen = '#CCFF00'
    orange = '#FFBF00'
    yellow = '#FFFF00'
    red = '#FF0000'
    white = '#FFFFFF'
    ascolor = white
    
    lines = []
    lines += ['<html>']
    lines += ['<head>']
    lines += ['<title>']
    lines += ['</title>']
    lines += ['</head>']
    lines += ['<body>']
    lines += ['<h1>APOGEE Abundances</h1>']
    lines += ['-10 &lt; b &lt; +10 deg, logg<1.5, Teff<5000, STARFLAG=0, DR17 25279 stars<br>']
    lines += ['Click on figure to enlarge']
    lines += ['<table border=1>']
    lines += ['<tr><th>Number</th>']
    lines += ['<th>Element</th><th>DR17 Abund</th><th>DR17 Abund (Teff)</th>']
    lines += ['<th>Error</th><tr>']
    for i in range(len(elem)):
        el = elem[i]
        name = el.upper()
        if len(el)==2:
            name = el[0].upper()+el[1]
        print(i+1,el)
        lines += ['<tr>']
        lines += ['<td style="text-align:center">{:}</td>'.format(i+1)]
        lines += ['<td style="text-align:center">{:}</td>'.format(name)]
        width = 500
        lines += ['<td><a href="plots/dr17_{:s}_feh_hist2d.png" target="_blank"><img src="plots/dr17_{:s}_feh_hist2d.png" width={:d}></a></td>'.format(el,el,width)]
        lines += ['<td><a href="plots/dr17_{:s}_feh_teff.png" target="_blank"><img src="plots/dr17_{:s}_feh_teff.png" width={:d}></a></td>'.format(el,el,width)]
        lines += ['<td><a href="plots/dr17_{:s}_scatter_snr.png" target="_blank"><img src="plots/dr17_{:s}_scatter_snr.png" width={:d}></a></td>'.format(el,el,width)]
        lines += ['</tr>']
    lines += ['</table>']
    lines += ['</body>']
    lines += ['</html>']
    dln.writelines('dr17_precision.html',lines)
