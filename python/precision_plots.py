import os
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from dlnpyutils import utils as dln,plotting as pl
from scipy.stats import binned_statistic


def precision_plots():
    """ Compare APOGEE DR17 abundance precision to SDSS-V Astra ASPCAP values. """
    
    #astraAllStarASPCAP-0.8.fits.gz
    #In [8]: len(tab)
    #Out[8]: 2466006

    #gd,=np.where((np.abs(tab['b'])<5) & ((tab['l']<60) | (tab['l']>300)))

    #In [12]: len(gd)
    #Out[12]: 583909

    #o=pl.hist2d(tab['fe_h'][gd],tab['alpha_m_atm'][gd],log=True,xr=[-3,1],yr=[-1,1])

    #o=pl.hist2d(tab['fe_h'][gd],tab['alpha_m_atm'][gd],log=True,xr=[-2.5,1],yr=[-0.75,1.0])

    #o=pl.hist2d(tab['fe_h'][gd],tab['coarse_alpha_m_atm'][gd],log=True,xr=[-2.5,1],yr=[-0.75,1.0])
    #looks pretty horrible

    #o=pl.hist2d(tab['teff'][gd],tab['logg'][gd],log=True)

    #gd,=np.where((np.abs(tab['b'])<5) & ((tab['l']<60) | (tab['l']>300)) &
    #             (tab['teff']<6300) & (tab['logg']<3.2))


    #In [27]: len(gd)
    #Out[27]: 513735

    #o=pl.hist2d(tab['fe_h'][gd],tab['mg_h'][gd]-tab['fe_h'][gd],log=True,xr=[-2.5,1],yr=[-0.75,1.0])

    #dr17=Table.read('/Users/nidever/apogee2/spectro/allStarLiter-dr17-synspec_rev1.fits.gz')
    #for c in dr17.colnames:dr17[c].name=c.lower()

    # xmatch with dr17 to make sure 

    # tab['sdss4_apogee_id']=tab['sdss4_apogee_id'].data.data
    # _,ind1,ind2=np.intersect1d(dr17['apogee_id'],tab['sdss4_apogee_id'][gd],return_indices=True)

    #In [50]: len(ind1)
    #Out[50]: 40218

    # dr17b = dr17[ind1]
    # tab2 = tab[gd][ind2]
    # tab2['mg_fe'] = tab2['mg_h']-tab2['fe_h']


    # tab2.write('/Users/nidever/as5/hge/abundances/astraAllStarASPCAP-0.8_mwmdisk_dr17matched.fits')
    # dr17b.write('/Users/nidever/as5/hge/abundances/allStarLiter-dr17-synspec_rev1_mwmdisk_sdss5astramatched.fits')
    
    # -0.054374

    # o=pl.hist2d(dr17b['fe_h'],dr17b['mg_fe'],log=True,xr=[-2.5,1.0],yr=[-0.5,0.5])

    # o=pl.hist2d(tab2['fe_h'],tab2['mg_fe']+0.0544,log=True,xr=[-2.5,1.0],yr=[-0.5,0.5])


    # o=pl.hist2d(dr17b['fe_h'],dr17b['c_fe'],log=True,xr=[-2.5,1.0],yr=[-0.5,0.5])

    # o=pl.hist2d(tab2['fe_h'],tab2['c_h']-tab2['fe_h'],log=True,xr=[-2.5,1.0],yr=[-0.5,0.5])
    # this is horrible!


    # o=pl.hist2d(dr17b['fe_h'],dr17b['n_fe'],log=True,xr=[-2.5,1.0],yr=[-0.5,0.75])

    # o=pl.hist2d(tab2['fe_h'],tab2['n_h']-tab2['fe_h'],log=True,xr=[-2.5,1.0],yr=[-0.5,0.75])
    # this actually looks tighter

    # o=pl.hist2d(dr17b['fe_h'],dr17b['ni_fe'],log=True,xr=[-2.5,1.0],yr=[-0.5,0.75])
    # o=pl.hist2d(tab2['fe_h'],tab2['ni_h']-tab2['fe_h'],log=True,xr=[-2.5,1.0],yr=[-0.5,0.75])
    # more outliers

    # o=pl.hist2d(dr17b['fe_h'],dr17b['si_fe'],log=True,xr=[-2.5,1.0],yr=[-0.5,0.75])
    # o=pl.hist2d(tab2['fe_h'],tab2['si_h']-tab2['fe_h'],log=True,xr=[-2.5,1.0],yr=[-0.5,0.75])
    # slightly worse systematics

    # o=pl.hist2d(dr17b['fe_h'],dr17b['co_fe'],log=True,xr=[-2.5,1.0],yr=[-0.5,0.75])
    # o=pl.hist2d(tab2['fe_h'],tab2['co_h']-tab2['fe_h'],log=True,xr=[-2.5,1.0],yr=[-0.5,0.75])
    # somewhat worse

    astra = Table.read('/Users/nidever/as5/hge/abundances/astraAllStarASPCAP-0.8_mwmdisk_dr17matched.fits')
    astra = astra.filled()
    dr17 = Table.read('/Users/nidever/as5/hge/abundances/allStarLiter-dr17-synspec_rev1_mwmdisk_sdss5astramatched.fits')
    dr17 = dr17.filled()
    

    #astra 21
    #elem = ['al_h', 'ca_h', 'ce_h', 'c_h', 'co_h', 'cr_h', 'cu_h', 'fe_h', 'k_h', 'mg_h',
    #        'mn_h', 'na_h', 'nd_h', 'ni_h', 'n_h', 'o_h', 'p_h', 'si_h', 's_h', 'ti_h', 'v_h']


    # dr17 elements
    #elem2 = ['c_fe', 'ci_fe', 'n_fe', 'o_fe', 'na_fe', 'mg_fe', 'al_fe', 'si_fe', 'p_fe', 's_fe',
    #         'k_fe', 'ca_fe', 'ti_fe', 'tiii_fe', 'v_fe', 'cr_fe', 'mn_fe', 'co_fe', 'ni_fe', 'cu_fe', 'ce_fe']


    # ci ?
    # tiii ?
    # fe  ?

    #elem = ['al', 'ca', 'ce', 'c', 'co', 'cr', 'cu', 'k', 'mg',
    #        'mn', 'na', 'nd', 'ni', 'n', 'o', 'p', 'si', 's', 'ti', 'v']

    #elem = ['c','n','o','na','mg','al','si','p','s','k','ca','ti','v','cr','mn','co','ni','cu','ce','nd']
    elem = ['c','n','o','na','mg','al','si','s','k','ca','ti','v','cr','mn','co','ni','ce']
    
    # 20 species: C, C I, N, O, Na, Mg, Al, Si, S, K, Ca, Ti, Ti II, V, Cr, Mn, Fe, Co, Ni, and Ce.

    plt.figure(1,figsize=(10,7))
    
    xr = [-2.5,1.0]
    yr = [-0.5,0.75]
    fehbins = np.arange(-0.8,0.6,0.05)
    snrbins = np.arange(20,100,10)
    
    # Element plots
    for i in range(len(elem)):
        el = elem[i]
        name = el.upper()
        if len(el)==2:
            name = el[0].upper()+el[1]
            
        # DR17 plot
        abund1 = dr17[el+'_fe']
        gd1, = np.where((np.abs(abund1) < 10) & np.isfinite(abund1))
        o=pl.hist2d(dr17['fe_h'][gd1],abund1[gd1],log=True,xr=xr,yr=yr,xtitle='[Fe/H]',
                    ytitle='['+name+'/Fe]',title='APOGEE DR17 - '+name,charsize=17)

        res,xedge,_ = binned_statistic(dr17['fe_h'][gd1],abund1[gd1],bins=fehbins,statistic=np.nanmedian)
        sm1 = dln.interp(xedge[:-1],res,dr17['fe_h'][gd1])
        diff1 = abund1[gd1]-sm1
        sig1 = dln.mad(diff1)
        plt.plot(xedge[:-1],res,c='r')
        plt.annotate('Scatter = {:.3f} dex'.format(sig1),xy=[0.05,0.05],xycoords='axes fraction',ha='left',fontsize=16)
        plt.ylim(yr)
        plt.savefig('plots/dr17_'+el+'_feh.png',bbox_inches='tight')
        
        
        # Astra ASPCAP plot
        abund2 = astra[el+'_h']-astra['fe_h']
        gd2, = np.where((np.abs(abund2) < 10) & np.isfinite(abund2))
        med = np.nanmedian(abund2[gd2]-dr17[el+'_fe'][gd2])
        calabund = abund2-med
        o=pl.hist2d(astra['fe_h'][gd2],calabund[gd2],log=True,xr=[-2.5,1.0],yr=[-0.5,0.75],xtitle='[Fe/H]',
                    ytitle='['+name+'/Fe]',title='Astra ASPCAP - '+name,charsize=17)
        res2,xedge2,_ = binned_statistic(astra['fe_h'][gd2],calabund[gd2],bins=fehbins,statistic=np.nanmedian)
        sm2 = dln.interp(xedge2[:-1],res2,astra['fe_h'][gd2])
        diff2 = calabund[gd2]-sm2
        sig2 = dln.mad(diff2)
        plt.plot(xedge2[:-1],res2,c='r')
        plt.annotate('Scatter = {:.3f} dex'.format(sig2),xy=[0.05,0.05],xycoords='axes fraction',ha='left',fontsize=16)
        plt.ylim(yr)
        plt.annotate('{:.3f} dex offset'.format(med),xy=[0.05,0.93],xycoords='axes fraction',ha='left',fontsize=16)
        plt.savefig('plots/astra_'+el+'_feh.png',bbox_inches='tight')
        
        # Scatter vs. S/N plot
        
        res,xedge,_ = binned_statistic(dr17['snr'][gd1],diff1,bins=snrbins,statistic=dln.mad)
        res2,xedge2,_ = binned_statistic(astra['snr'][gd2],diff2,bins=snrbins,statistic=dln.mad)
        pl.plot(xedge[:-1],res,c='r',label='DR17',xtitle='S/N',ytitle='Scatter',title='')
        plt.plot(xedge2[:-1],res2,c='blue',label='Astra')
        plt.legend()
        plt.savefig('plots/scatter_'+el+'_feh.png',bbox_inches='tight')
        
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
    lines += ['-5 &lt; b &lt; +5 deg, 300 &lt; l &lt; 60, DR17-Astra xmatched, 40218 stars<br>']
    lines += ['Click on figure to enlarge']
    lines += ['<table border=1>']
    lines += ['<tr><th>Number</th>']
    lines += ['<th>Element</th><th>DR17</th>']
    lines += ['<th>Astra</th><th>Scatter</th><tr>']
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
        lines += ['<td><a href="plots/dr17_{:s}_feh.png" target="_blank"><img src="plots/dr17_{:s}_feh.png" width={:d}></a></td>'.format(el,el,width)]
        lines += ['<td><a href="plots/astra_{:s}_feh.png" target="_blank"><img src="plots/astra_{:s}_feh.png" width={:d}></a></td>'.format(el,el,width)]
        lines += ['<td><a href="plots/scatter_{:s}_feh.png" target="_blank"><img src="plots/scatter_{:s}_feh.png" width={:d}></a></td>'.format(el,el,width)]
        lines += ['</tr>']
    lines += ['</table>']
    lines += ['</body>']
    lines += ['</html>']
    dln.writelines('precision.html',lines)
