# /usr/bin/env/python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 08:32:01 2020

@author: cluzetb
"""


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from postes.explore_metadata import find_metadata
import seaborn as sns


def find_common_massif_obsbase(obsPath, metadata = '/home/cluzetb/snowtools_git/DATA/METADATA.xml', massifs = [12], droptypes = [], pivotval = 'HTN cm'):
    df = pd.read_csv(obsPath, delimiter = ';')

    # A. PREPARE the Dataframe
    # 2. to_datetime
    df['Date UTC'] = pd.to_datetime(df['Date UTC'], format = '%Y-%m-%d-%H-%M')

    # drop climatological posts
    for dt in droptypes:
        df = df[df.type != dt]
    # 1. pivot on posts
    df = df.pivot(index = 'Date UTC', columns = 'NUMPOST', values = pivotval)

    # /!\ CRUCIAL POINT : sort the columns (NUMPOINT) so that the points are in the same ordre in the observations and in the forcing.
    print('sooooorting df')
    df = df.sort_index(axis=1)
    # 3. keep only last years, 6:00 obs and winter
    df = df.drop([ d for d in df.index if ((d.hour != 6) or (7 <= d.month <= 9))])
    df = df.drop([ d for d in df.index if ((d.year < 2013) or (d.year > 2018))])

    # 4. drop columns with less than 40% points
    df = df.loc[:, df.isnull().mean() <= .95]

    # 4. convert to meters -_-:
    if pivotval == 'HTN cm':
        df = df / 100.
        df.rename(columns={'HTN cm': 'HTN m'})

    # sort the columnsf

    # B. select in the DF

    # 1'. select posts from massif
    listp, listpnb, listpname = find_metadata(filepath=metadata, massifs = massifs)

    listcols_in_sel = [l for l in listpnb if l in df.columns]

    df = df[listcols_in_sel]

    # Actualize the list of posts
    listp = [p for (p, nb) in zip(listp, listpnb) if nb in df.columns]
    listpname = [na  for (nb, na) in zip(listpnb, listpname) if nb in df.columns]
    listpnb = [nb for nb in listpnb if nb in df.columns]
    print(listpnb)
    print('common postes between metadata and csv :', len(listpnb))
    # print(listpname)

    return df, listp, listpnb, listpname
