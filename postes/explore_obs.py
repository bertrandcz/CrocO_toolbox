# /usr/bin/env/python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 08:32:01 2020

@author: cluzetb
"""


import datetime
import os

import numpy as np
import pandas as pd
from postes.explore_metadata import find_metadata


def find_common_massif_obsbase(obsPath, metadata = os.environ['SNOWTOOLS_CEN'] + '/DATA/METADATA.xml', massifs = [12], droptypes = [], pivotval = 'HTN cm'):
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
    listpname = [na for (nb, na) in zip(listpnb, listpname) if nb in df.columns]
    listpnb = [nb for nb in listpnb if nb in df.columns]
    print(listpnb)
    print('common postes between metadata and csv :', len(listpnb))
    # print(listpname)

    return df, listp, listpnb, listpname


def detect_windy_dates(obsTs, oper, year, sda, sdf):
    """
    BC 29/05/20
    detect wind ablation from observation (df or unpickled obsTs) and oper HS timeseries pickle and mask the observations after it.
    The vague formatting for observations comes from the fact that Preparator needs a df to produce daily observations
    while the wind detection was conceived to work with unpickled obsTs and thatI don't have enough time to harmonize

    An event is suspect if the rate of HS ablation between two consecutive daily observations (separated by less than 4 days):
      _
    _|  exceeds sda [m/day] (about .12 m/day)
     |_ is superior to sdf (about 5) times the HS ablation rate in the model during the same period.

    """

    # iterate over posts
    perc_ok = np.zeros((obsTs['DEP'].shape[1])) - 1
    suspect_beg_dates = dict()
    suspect_end_dates = dict()

    for ip in range(obsTs['DEP'].shape[1]):
        suspect_beg_dates[ip] = []
        suspect_end_dates[ip] = []
        # for ip in range(startp, endp):
        obsHtn = np.ma.masked_invalid(obsTs['DEP'][:, ip])
        time_obs = np.ma.masked_array(
            obsTs['time'], mask=obsHtn.mask).compressed()
        time_mod = oper['time']
        time_obs_mod, iobs, imod = np.intersect1d(time_obs, time_mod,
                                                  return_indices=True,
                                                  assume_unique=True)
        obsHtn = obsHtn.compressed()[iobs]
        modHtn = oper['DEP'][imod, ip, 0]
        if len(time_obs_mod) > 0:
            for i, t in enumerate(time_obs_mod[0:-1]):
                suspect = False
                t_n = time_obs_mod[i + 1]
                dt = (t_n - t).days
                if dt < 4.:
                    D_dep_tol = sda * dt
                    D_dep_obs = obsHtn[i] - obsHtn[i + 1]
                    if len(modHtn) > 0:
                        D_dep_mod = modHtn[i] - modHtn[i + 1]
                        modOk = True
                    else:
                        modOk = False
                    if D_dep_obs > D_dep_tol:
                        if modOk:
                            if D_dep_obs > sdf * D_dep_mod:
                                suspect = True
                        else:
                            suspect = True
                if suspect is True:
                    suspect_beg_dates[ip].append(t)
                    suspect_end_dates[ip].append(t_n)
            startD = datetime.datetime(year, 10, 1, 0, 0, 0)
            endD = datetime.datetime(year + 1, 7, 1, 0, 0, 0)
            if len(time_obs_mod) > 11:
                if len(suspect_beg_dates[ip]) > 0:
                    # one suspect date has been detected
                    perc_ok[ip] = (suspect_beg_dates[ip][0] - startD) / (endD - startD) * 100.
                else:
                    # no detection : 100% ok.
                    perc_ok[ip] = 100.
            else:
                # reject :  not enough obs this year.
                perc_ok[ip] = np.nan
        else:
            # reject : no obs this year.
            perc_ok[ip] = np.nan

    return suspect_beg_dates, suspect_end_dates, perc_ok
