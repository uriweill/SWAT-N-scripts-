import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress 
import string
from itertools import product
import matplotlib.patches as patches
import matplotlib

# set font of figures
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"

# set days to use (remove day 0) and to generations - given 10 generations per day
# number of generations between day zero and one is 9 generations because of serial dilution

data_points = [1,2,3,4]
generations = [0.,9.,19.,29., 39.]

generations = generations[-len(data_points):]
generations = [x-generations[0] for x in generations]
    
letters = list(string.ascii_uppercase[0:8])
numbers = [str(i+1).zfill(2) for i in range(12)]

wells = ["".join(x) for x in product(letters, numbers)]

# parse plate information and pepare dataframe for analysis
plate = pd.DataFrame.from_csv("plates.csv").unstack().reset_index()
plate.columns = ['replicate', 'ORF', 'location']
plate['plate'] = plate['location'].apply(lambda x: x[:2].replace('0',''))
plate['well'] = plate['location'].apply(lambda x: x[2:])
FACS = []
# iterate over all csv files and read the FACS data in the form of fitness
for i in range(3):
    df = pd.read_csv("plate%i_mcherry_ratio.csv"%(i+1), header=None).T
    #remove samples in which the initial mcherry/gfp ratio is too off.
    df[df[0] > 5 ] = np.nan
    df[df[0] < 0.2 ] = np.nan    
    # normalize cell counts to t=0
    df_N = df.div(df.loc[:,data_points[0]], axis=0)
    df_N = df_N.apply(np.log2)
    
    df_N['well'] = wells
    df_N['plate'] = '%i'%(i+1)
    
    FACS.append(df_N)    

# combine numeric data with data abaut each well (strains, localization. etc)
FACS = pd.concat(FACS)
DF = FACS.merge(plate)
localization = pd.read_csv("localization.csv")
DF = localization.merge(DF)

# pefrom linear regression
reg = DF[data_points].apply(func=lambda x: linregress(generations, x), axis=1)
reg = pd.DataFrame(reg)
reg[['slope', 'intercept', 'r-value', 'p-value', 'stderr']] = reg[0].apply(pd.Series)
reg.drop(0, 1, inplace=True)
reg['d_mu[%]'] = (reg['slope']*100)
DF = DF.merge(reg, left_index=True, right_index=True)
DF.groupby(['ORF', 'Gene', 'N', 'C']).mean().to_csv("results_df.csv")
groups = DF.groupby('ORF')

# set global params
variability = 0.5
rsq = 0.5
replicates = 2
noise_range = 0.015 

# plot fintess graph for each gene
genes = []
for n, (i, ORF) in enumerate(groups):
    gene = ORF.Gene.values[0][0]+ORF.Gene.values[0][1:].lower()
    mu = ORF.mean()['d_mu[%]']
    Nloc = ORF.N.values[0]
    Cloc = ORF.C.values[0]
    y_mean = ORF[data_points].mean(axis=0)
    y_std = ORF[data_points].std(axis=0)
    cv = np.abs(y_std / y_mean).max()

    genes.append(gene)
    plt.figure(0, figsize=(6,6))
    ax1 = plt.subplot2grid((3,3), (0,2), colspan=1)
    ax2 = plt.subplot2grid((3,3), (1,2), colspan=1)
    ax = plt.subplot2grid((3,3), (0, 0), rowspan=2, colspan=2)
    
    f = lambda x: noise_range*(x-generations[0])
    x0 = np.arange(generations[0],generations[-1]+1,0.01)
    ax.fill_between(x0,f(x0),-f(x0), facecolor='0.9', zorder=0, edgecolor='0.7', linestyle=':', linewidth=0.5)

    a, b = ORF.mean()[['slope', 'intercept']]
    if abs(a) > noise_range:
        c='k'
    else:
        c='0.7'
    ax.annotate('$ \Delta \mu  = %.1f$%s' %(mu, '%'), (17,0), 
                ha='center', va='center', size=13)
        
    ax.plot(x0, a*x0+b, c=c)
    ax.errorbar(generations, y_mean, yerr=y_std, fmt='o',
                c=c, mec=c, mfc='w', zorder=2,  capsize=3, lw=1.5, ms=3)
    
    ax.set_yticks([-1,0,1])
    ax.set_yticklabels(['0.5', '1', '2'])
    ax.set_xlabel('Generations')
    ax.set_xlim(generations[0],generations[-1]+1)
    ax.set_ylim(-1.5,1.5)
    
    fontsize=13
    ax.tick_params(direction='out', labelsize=fontsize)
    ax.set_xlabel(ax.get_xlabel(), size=fontsize)
    ax.set_ylabel(ax.get_ylabel(), size=fontsize)
    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()
    
    ax.set_ylabel(r'$\frac{C_{tag}}{N_{tag}}$', rotation=0, va='center', size=15)
    
    ax.add_patch(patches.Rectangle((1, -1.3),3.3,0.25, fill=True, ec='k', fc='#00ff00', zorder=10))
    ax.add_patch(patches.Rectangle((4.7, -1.3),7,0.25, fill=True, ec='k', fc='w', zorder=10))
    ax.annotate(gene, (8.2,-1.2), size=11, zorder=10, ha='center', va='center')
    ax.annotate('GFP', (2.65,-1.2), size=11, zorder=10, ha='center', va='center')
        
    ax.add_patch(patches.Rectangle((8.32, 1.1),3.3,0.25, fill=True, ec='k', 
                                   fc='#00ff00', zorder=10))
    ax.add_patch(patches.Rectangle((1, 1.1),7,0.25, fill=True, ec='k', fc='w', zorder=10))
    ax.annotate(gene, (4.8,1.2), size=11, zorder=10, ha='center', va='center')
    ax.annotate('GFP', (10,1.2), size=11, zorder=10, ha='center', va='center')

    try:
        figure_gene_name = gene[0]+gene[1:].lower()
        Ctag = plt.imread('pictures/%s-GFP.tif'%figure_gene_name)
        Ntag = plt.imread('pictures/GFP-%s.tif'%figure_gene_name)
        ax1.imshow(Ctag)
        ax2.imshow(Ntag)
            
    except IOError:
        print("missing image for " + gene)
        
    
    cweight='normal';nweight='normal'
    if mu > noise_range*100:
        cweight = 'bold' 
    elif mu < -noise_range*100:
        nweight = 'bold'
    ax1.set_title("C:"+Cloc, weight=cweight)
    ax2.set_title("N:"+Nloc, weight=nweight)
    
    ax1.set_xticks([])
    ax1.set_yticks([])

    ax2.set_xticks([])
    ax2.set_yticks([])
    plt.savefig('results/%s.svg'%gene)
    plt.close()

# export results
df = DF[DF['Gene'].isin(genes)]
df = df.groupby(['Gene', 'N', 'C']).mean()
df.to_csv("final_genes.csv")

dmu = pd.read_csv("results_df.csv")
uncorrupted = pd.read_csv('uncorrupted.csv')
dmu = dmu.merge(uncorrupted, how='inner')
dmu = dmu['d_mu[%]'].dropna()
plt.figure()
ax = plt.axes()
plt.hist(dmu.dropna(), bins=np.arange(-15,15,1.5))
ax.axvline(1.5, c='0.7', ls=':')
ax.axvline(-1.5, c='0.7', ls=':')
ax.axvspan(-1.5, 1.5, color='0.7', alpha=0.4)
ax.set_xticks(np.arange(-10,11,2))
ax.set_xlim(-9,9)
plt.savefig('results/dmu histogram.svg')

