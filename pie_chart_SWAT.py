# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 09:53:45 2017

@author: dan
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
###   FILES   ###
location_mapping_file = 'loc_code.csv'
loc_code = pd.Series.from_csv(location_mapping_file).to_dict()

data = 'SWAT FG Library updated annotations V2.csv'

### COLUMN NAMES (OF DIFFERENT TAGGINGS)   ###
N_G = "NOP1pr-GFP localization new Generic"
C   = "Known C' LocalizationWithambiguous"
N   = "Native-GFP localization new Seamless"
Cherry = 'TeF2-mCherry localization new cherry'

tags = [N_G, C, N, Cherry]

ignore = ['ambiguous', 'below threshold', 'missing']
def tidy_split(df, column, sep=',', keep=False):
    """
    Split the values of a column and expand so the new DataFrame has one split
    value per row. Filters rows where the column is missing.

    Params
    ------
    df : pandas.DataFrame
        dataframe with the column to split and expand
    column : str
        the column to split and expand
    sep : str
        the string used to split the column's values
    keep : bool
        whether to retain the presplit value as it's own row

    Returns
    -------
    pandas.DataFrame
        Returns a dataframe with the same columns as `df`.
    """
    indexes = list()
    new_values = list()
    df = df.dropna(subset=[column])
    for i, presplit in enumerate(df[column].astype(str)):
        values = presplit.split(sep)
        if keep and len(values) > 1:
            indexes.append(i)
            new_values.append(presplit)
        for value in values:
            indexes.append(i)
            new_values.append(value)
    new_df = df.iloc[indexes, :].copy()
    new_df[column] = new_values
    return new_df

def replace_letter_2_word(series):
    out = []
    for i, letters in series.iteritems():
        L = list(letters)
        try:
            locs = ",".join([loc_code[l.strip()] for l in L])
            out.append(locs)
        except:
            out.append(letters)
            print "unable to replace %s in row %i (%s)" %(letters,i+2, series.name)
    return out
    
def map_annotation(location_mapping_file, data):
    
    df = pd.read_csv(data, sep='\t')
    df.replace(np.isnan, '0', inplace=True)
#    df = df[[N_G, N, Cherry]].str.strip()
#    f = lambda x: ','.join(loc_code[list(x)].values)
#    df[[N_G, N, Cherry]] = df[[N_G, N, Cherry]].apply(f)
    df[N_G] = replace_letter_2_word(df[N_G])
    df[N] = replace_letter_2_word(df[N])
    df[Cherry] = replace_letter_2_word(df[Cherry])
        
    df.replace({C:{'0':'missing'}}, inplace=True)
    df.to_csv('data_with_mapping.csv')
    
    return df

def compare_localizations(tagX, tagY,
                          ignore=ignore):
    
    # convert localization code letters to full name localizations
    df = pd.DataFrame.from_csv('data_with_mapping.csv')    
    df[df.Gene=='None'] = np.nan
    df.dropna(subset=['Gene'], inplace=True)
    df.replace('bud', 'bud neck', inplace=True)
    # if a specific gene is associated with more than a single cellular location
    # split this gene into several rows according to the number of locations    
    X_loc = tidy_split(pd.DataFrame(df[tagX]), tagX)
    Y_loc = tidy_split(pd.DataFrame(df[tagY]), tagY)
        
    # remove 'ignore' locations - this can be editted - see function params    
    X_loc = X_loc[~X_loc.isin(ignore)].dropna()
    Y_loc = Y_loc[~Y_loc.isin(ignore)].dropna()
    
    Xonly = set(X_loc.index).difference(set(Y_loc.index))
    Yonly = set(Y_loc.index).difference(set(X_loc.index))
    
    Xloc = X_loc.copy()
    Xloc['ORF'] = X_loc.index
    
    Yloc = Y_loc.copy()
    Yloc['ORF'] = Y_loc.index
    
    same = set(Xloc.merge(Yloc, left_on=['ORF',tagX], right_on=['ORF', tagY], 
                          how='inner').ORF.values)
    
    ALL = set(X_loc.index).union(set(Y_loc.index))
    diff = ALL.difference(same.union(Xonly.union(Yonly)))
    
    return Xonly, Yonly, same, diff

def draw_pie_chart(tagX=C, tagY=N, ignore=ignore,
                   labels = ['same', 'N only', 'different', 'C only'],
                    colors=['','','','']):
    
    Xonly, Yonly, same, diff = compare_localizations(tagX, tagY, ignore=ignore)
    sizes = map(len, [same, Yonly, diff, Xonly])
    explode = (0, 0, 0, 0)  # only "explode" the 2nd slice (i.e. 'Hogs')
    
    fig1, ax1 = plt.subplots()
    total = sum(sizes)
    ax1.pie(sizes, explode=explode, labels=labels, 
            autopct=lambda(p): '{:.0f}'.format(p * total / 100),
            shadow=True, startangle=90, colors=colors)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    s = tagX + " vs. " + tagY
    ax1.set_title(s)
    plt.tight_layout()
    plt.savefig("pie[%s].png"%s)
    plt.show()

    return sizes

df = map_annotation(location_mapping_file, data)

for tagX, tagY in combinations(tags, 2):
        
    Xonly, Yonly, same, diff = compare_localizations(tagX, tagY)
    sizes = draw_pie_chart(tagX, tagY, 
                   labels=['Same', 'Y only', 'different', 'X only'],
                   colors =['#f3e6ff','#dab3ff', '#c180ff', '#a94dff'])
    
    s = tagX + " vs. " + tagY
    df[s] = ''
    df[s].loc[diff] = 'different'
    df[s].loc[same] = 'same'
    df[s].loc[Yonly] = 'Y only'
    df[s].loc[Xonly] = 'X only'

df.to_csv("out.csv")

