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
from postes.explore_metadata import find_metadata, find_duplicates,\
    find_name_station
import matplotlib.pyplot as plt


def find_common_massif_obsbase(obsPath, metadata = os.environ['SNOWTOOLS_CEN'] + '/DATA/METADATA.xml',
                               massifs = [12], droptypes = [], pivotval = 'HTN cm',
                               startY=2013, endY = 2018,
                               percent_valid = 5, percent_zeros = 90,
                               dropCorsica = False, remove_duplicate_posts = False, plotClean = False):
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
    df = df.drop([ d for d in df.index if ((d.year < startY) or (d.year > endY))])

    # 4. drop columns with less than percent% points valid
    df = df.loc[:, df.isnull().mean() <= (1 - 0.01 * percent_valid)]

    # 4. convert to meters -_-:
    if pivotval == 'HTN cm':
        df = df / 100.
        df.rename(columns={'HTN cm': 'HTN m'})
    # BC august 2020 drop columns with too many zeros among valid values
    # (eliminate exotic posts with only a few non-null values.
    df = df.loc[:, (df < 0.001).mean() / (1 - df.isnull().mean()) <= 0.01 * percent_zeros]

    # BC august 2020 drop corsica
    # NB: Posts from Italy, starting in '204' are also removed (e.g. PRERICHARD wich is not trustworthy.)
    not_corsica = [c for c in df.columns if not str(c).startswith('20')]
    if dropCorsica is True:
        df = df[not_corsica]
    # sort the columns

    # 1'. select posts from massif
    listp, listpnb, listpname = find_metadata(filepath=metadata, massifs = massifs)

    listcols_in_sel = [l for l in listpnb if l in df.columns]

    df = df[listcols_in_sel]

    # remove the duplicated posts.
    # Actualize the list of posts
    listp = [p for (p, nb) in zip(listp, listpnb) if nb in df.columns]
    listpname = [na for (nb, na) in zip(listpnb, listpname) if nb in df.columns]
    listpnb = [nb for nb in listpnb if nb in df.columns]

    if remove_duplicate_posts:
        _, numbers, _ = find_duplicates(listp, listpnb, listpname)
        print(numbers)
        dups_to_remove = duplicates_to_remove(numbers)
        print(dups_to_remove, len(dups_to_remove))
        print(' posts to remove: ', find_name_station(dups_to_remove))
        not_to_remove = [c for c in df.columns if c not in dups_to_remove]
        print( 'df shape before removing the dups', df.shape)
        df = df[not_to_remove]
        print( 'df shape after removing the dups', df.shape)

    # Actualize the list of posts
    listp = [p for (p, nb) in zip(listp, listpnb) if nb in df.columns]
    listpname = [na for (nb, na) in zip(listpnb, listpname) if nb in df.columns]
    listpnb = [nb for nb in listpnb if nb in df.columns]

    # clean the data
    df = clean_obs(df, plot = plotClean)
    print(listpnb)
    print('common postes between metadata and csv :', len(listpnb))

    return df, listp, listpnb, listpname


def duplicates_to_remove(numbers):
    """
    mask duplicate posts. General rule: reject the nivological post (00X). many exceptions and particular cases
    """
    # BC 25/08/20 : copy from mask_duplicate_posts.
    remlist = set()
    # all those post must be deleted (thereby transgressing the general rule), and their duplicates kept.
    exceptions_remove_list = [66146002,  # FONT NEGRE-EDFNIVO
                              9200602,   # MONT D'OLMES
                              5085001,   # MONTGENEVRE
                              6073400,   # Isola
                              6073405,   # ISOLA
                              31042012,  # LUCHON
                              ]
    # all the posts marked as duplicated of the following posts must be kept.
    exceptions_keep_list = [5157001,     # ST VERAN
                            38472002,    # COL DE PORTE-EDFNIVO
                            74056005,    # LE TOUR
                            73132003,    # COL-DES-SAISIES
                            74085406,    # LES COMTAMINES MONTJOI
                            74136003,    # GRAND BORNAND
                            73206001,    # PRALOGNAN LA VANOISE
                            65283001,    # LOUDERVIELLE
                            65138400,    # CAUTERETS_LYS_0
                            ]
    for _, nums in numbers.items():
        # exceptions
        excs = set(nums) & set(exceptions_remove_list)
        kexcs = set(nums) & set(exceptions_keep_list)
        if len(excs) > 0:
            for e in excs:
                remlist.add(e)
        elif len(kexcs) > 0:
            pass
        else:
            # remove nivological posts only if a corresponding clim post is in its list of dups.
            climposts = [str(ggg)[-3] == '0' for ggg in nums]
            print(nums, climposts)
            if any(climposts):
                for ic, c in enumerate(climposts):
                    if c is False:
                        remlist.add(nums[ic])
    return list(remlist)


def clean_obs(df, plot = False):
    """
    BC remove false zeros in EDFNIVO timeseries and correct spurious errors including grass at the beg/end of season.

    """
    mows = mower_edf()
    for iiii, c in enumerate(df.columns):

        zeros = np.where(df[c] < 0.01, True, False)
        fifties = np.where(df[c] > 0.40, True, False)

        # prevent from false detection on first date.
        zeros[0] = False
        zeros[1] = False
        pos_suspect = []

        # detect drops
        pos_suspect_drop = [i for i, _ in enumerate(zeros)
                            if (zeros[i] and fifties[i - 1] and fifties[i - 2]) and (df.index[i] - df.index[i - 2] < pd.Timedelta(3, unit = 'd'))
                            ]
        # detect long drops
        pos_suspect_drop_continued = pos_suspect_drop[:]  # soft copy
        for i, z in enumerate(zeros):
            if z and i - 1 in pos_suspect_drop_continued and not (df.index[i].month == 10 and df.index[i].day == 1):
                pos_suspect_drop_continued.append(i)

        # detect peaks from zero
        pos_suspect_peak_from_zero = [i for i, _ in enumerate(zeros)
                                      if (zeros[i] and fifties[i - 1] and zeros[i - 2]) and (df.index[i] - df.index[i - 2] < pd.Timedelta(3, unit = 'd'))]

        # detect sudden relative peaks
        pos_suspect_rel_peak = [i for i, _ in enumerate(zeros) if (df[c][i] - df[c][i - 1]) > 0.6 and (df[c][i] - df[c][i + 1]) > 0.6]

        # mow the lawn
        pos_mow = []
        if c in mows.keys():
            print('mowing', c, iiii)
            pos_mow = [i for i, dt in enumerate(df.index) if any([pd.Timestamp(mow[0]) < dt < pd.Timestamp(mow[1]) for mow in mows[c]])]
        pos_suspect = list(set(pos_suspect_rel_peak + pos_suspect_peak_from_zero + pos_suspect_drop_continued + pos_suspect_drop + pos_mow))

        df[c][pos_suspect] = np.nan

        if plot:
            if len(pos_suspect) > 0  or c in mows.keys():
                # plot
                plt.figure()
                plt.plot(df.index, df[c], ls = 'None', marker = '.')
                for p in pos_suspect:
                    plt.axvline(x=df.index[p], c = 'r', lw = 3, zorder = -1)

                title = str(c) + '_' + find_name_station(c).replace(" ", "_")
                plt.title(title)
                plt.gcf().autofmt_xdate()
                plt.grid(True)

                plt.savefig('/home/cluzetb/Bureau/postes/mauvaises_donnees/suspect_' + title + '.png')
                plt.close()

    return df


def mower_edf():
    """
    dubious data to remove, mainly grass at the beginning/end of the season on edf stations. left and right: strict limits.
    Also treat the inconsistency between the teo stations at la clusaz.
    """
    d = {
        # LA CLUSAZ
        74080403: [('2015-12-31', '2016-08-01')],  # same data as 74080400 in early 2016: no sense !
        # ROUGNOUS-VERDON-EDFNIV
        4006006: [('2016-07-31', '2017-12-11'), ('2018-06-04', '2018-10-30'), ('2019-05-22', '2020-01-01')],
        # PASSAUR-EDFNIVO
        4193008: [('2016-07-31', '2016-10-14'), ('2017-05-09', '2017-10-21'), ('2018-05-10', '2018-10-29'), ('2019-05-01', '2020-01-01')],
        # IZOARD-EDFNIVO
        5027002: [('2016-07-31', '2016-10-14'), ('2017-06-01', '2017-11-06'), ('2018-08-01', '2018-11-01'), ('2019-06-04', '2020-01-01')],
        # LAC_NOIR-EDFNIVO
        5063003: [('2016-07-31', '2016-10-14'), ('2017-05-19', '2017-11-06'), ('2018-05-30', '2018-10-28'), ('2019-06-07', '2020-01-01')],
        # CHARDONNET-EDFNIVO
        5093003: [('2016-07-31', '2016-11-05'), ('2017-06-12', '2017-11-06'), ('2018-06-14', '2018-10-29'), ('2019-06-15', '2020-01-01')],
        # ROUGNOUS_PRAPIC-EDFNIV
        5096003: [('2016-07-31', '2016-10-14'), ('2017-06-09', '2017-11-06'), ('2018-06-15', '2018-10-28'), ('2019-06-17', '2020-01-01')],
        # CEZANNE-EDFNIVO
        5101003: [('2016-07-31', '2016-10-13'), ('2017-05-25', '2017-11-07'), ('2018-06-09', '2018-10-28'), ('2019-05-24', '2020-01-01')],
        # LES_MARROUS-EDFNIVO
        5157003: [('2016-07-31', '2016-10-13'), ('2017-06-06', '2017-11-06'), ('2018-06-13', '2018-10-29'), ('2019-05-07', '2019-05-19'),
                  ('2019-06-14', '2020-01-01')],  # bug in May 2019
        # SANGUINIERE-EDFNIVO
        6056005: [('2016-07-31', '2017-11-06'), ('2018-05-05', '2018-10-30'), ('2019-05-05', '2020-01-01')],
        # ALBEILLE-EDFNIVO
        9162001: [('2018-07-31', '2018-10-29'), ('2019-06-05', '2020-01-01')],
        # LES SONGES-EDFNIVO
        9220004: [('2016-07-31', '2020-01-01')],  # very bad data.
        # PRA_LONG-EDFNIVO
        31123005: [('2016-07-31', '2017-11-05'), ('2018-08-01', '2018-10-27'), ('2019-05-01', '2020-01-01')],
        # AGNELIN-EDFNIVO
        38163006: [('2019-07-01', '2020-01-01')],  # rather clean data
        # COL_DE_PORTE-EDFNIVO
        38472002: [('2019-05-01', '2020-01-01')],
        # MIGOUELOU
        65032004: [('2016-07-31', '2016-10-14'), ('2017-05-26', '2017-11-13'), ('2018-06-06', '2018-12-01'), ('2019-03-06', '2020-01-01')],  # clean except for March 2019 on.
        # TROUMOUSE-EDFNIVO
        65192006: [('2016-07-31', '2016-11-07'), ('2018-08-01', '2018-11-13'), ('2019-01-23', '2019-01-24'), ('2019-05-28', '2020-01-01')],
        # SPIGEOLES-EDFNIVO
        65282003: [('2016-07-31', '2016-11-05'), ('2017-06-22', '2017-11-06'), ('2018-08-01', '2018-11-11'), ('2019-06-25', '2020-01-01')],
        # BARRADA-EDFNIVO
        65295004: [('2017-11-01', '2018-10-29'), ('2019-06-21', '2020-01-01')],  # strange 2.5m in late 2018 season
        # LES_DOUGNES-EDFNIVO
        66005004: [('2016-07-31', '2016-10-11'), ('2017-05-08', '2017-11-06'), ('2018-05-29', '2018-10-28'), ('2019-05-11', '2020-01-01')],
        # GAOUGETA-EDFNIVO
        66081002: [('2016-07-31', '2016-11-06'), ('2017-05-23', '2017-11-06'), ('2018-06-15', '2018-11-01'), ('2019-06-05', '2020-01-01')],
        # CORMET_DE_ROSELEND-EDF
        73034008: [('2016-07-31', '2017-01-14'), ('2017-05-16', '2017-11-06'), ('2018-06-03', '2018-12-06'), ('2019-06-08', '2020-01-01')],  # strange start of the season
        # PLAN_SETI-EDFNIVO
        73040006: [('2017-06-18', '2017-10-23'), ('2018-06-29', '2018-10-28'), ('2019-06-28', '2020-01-01')],
        # PETIT_MONT_CENIS-EDFNI
        73056002: [('2017-05-29', '2017-10-04'), ('2018-06-04', '2018-10-28'), ('2019-05-27', '2020-01-01')],
        # BISSORTE-EDFNIVO
        73194005: [('2017-05-27', '2017-11-06'), ('2018-06-05', '2018-10-28'), ('2019-06-07', '2020-01-01')],
        # NOTRE_DAME_D'AOUT-EDFN
        73206003: [('2016-07-31', '2017-10-21'), ('2018-06-22', '2018-07-31'), ('2019-06-28', '2020-01-01')],  # very strange behaviour during 2016-2017
        # PETITE_GOUILLE_EDFNIVO
        73232004: [('2016-07-31', '2016-11-05'), ('2017-06-05', '2017-11-12'), ('2018-08-01', '2018-10-30'), ('2019-06-28', '2020-01-01')],
        # SOUS_LES_BARMES-EDFNIV
        73304008: [('2016-07-31', '2016-10-14'), ('2017-06-05', '2017-11-07'), ('2018-06-16', '2018-10-28'), ('2019-06-15', '2020-01-01')],
        # JUCLAR-EDFNIVO
        99130015: [('2016-07-31', '2016-11-07'), ('2017-04-19', '2018-10-31'), ('2019-01-10', '2019-02-04'), ('2019-04-29', '2020-01-01')],  # seems noisy in 2017-2018 and in the middle of last season
    }
    return d


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
