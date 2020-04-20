# -*- coding: utf-8 -*-
'''
Created on 6 févr. 2019

@author: cluzetb
'''
from SemiDistributed import PrepAbs
import os
from utilcrampon import dictErrors, niceName, niceLabel
from utilcrampon import setSubsetclasses, cm2inch
from utilpp import set_itimes

from matplotlib.colors import Normalize
from matplotlib.offsetbox import AnnotationBbox, OffsetImage

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np


if 'beaufix' not in os.uname()[1]:
    import seaborn as sns


class Pie(object):
    '''
    class to plot a sdObject into a pie
    '''

    def __init__(self, sdObj, focusCl = None, focusColors = None, maketitle = False, score = None):
        '''
        Constructor

        @param:
        focusColors : greylevel color for the focus dots

        '''

        if not os.path.exists('../pie/'):
            os.mkdir('../pie/')
        self.sdObj = sdObj
        if not self.sdObj.isloaded:
            self.sdObj.load()  # -> it is not always loaded (ex. of synth obs)
        self.sel = np.unique(self.sdObj.pgd.elev)
        self.ssl = sdObj.options.classesS
        self.sas = 'all'
        self.rmax = float(np.max(self.sdObj.pgd.elev)) + 300.

        self.dictlims = {'B1': [-0.2, 0.2], 'B2': [-0.2, 0.2], 'B3': [-0.2, 0.2], 'B4': [-0.2, 0.2],
                         'B5': [-0.2, 0.2], 'B6': [-0.2, 0.2], 'B7': [-0.2, 0.2], 'DEP': [-1., 1.]}
        if focusCl is not None:
            if type(focusCl) is int:
                self.focusCl = [focusCl]
            else:
                self.focusCl = focusCl
            if focusColors is not None:
                self.focusColors = focusColors
            else:
                self.focusColors = [0.] * len(self.focusCl)
            self.cmap_f = plt.cm.get_cmap(name = "gist_gray")
            self.cmap_f.set_bad('orange', 1.)
        self.maketitle = maketitle
        self.score = score

    def plot(self, ax = None, savefig = False, cmap = 'viridis', clims = None, title = None, colortitle = None):
        N = 8.
        width = 2. * np.pi / N

        if clims is not None:
            lims = clims
        else:
            lims = [0., 1.]
        self.title = title
        # if 0 and 20 slopes, 2 axes per var : 1 bar for flat, 1 pie for sloped
        for var in self.sdObj.options.ppvars:

            if ax is None:

                fig = plt.figure(figsize = (cm2inch(12, 3.5 * 12 / (1.5 + 2. * len(self.ssl) - 1 ) )))
                gs = gridspec.GridSpec(1, len(self.ssl) + 1,
                                       width_ratios=[0.8, 4, 4, 0.2] if len(self.ssl) == 3 else [0.8, 4, 0.2],
                                       wspace = 0.32,
                                       bottom = 0.15,
                                       top = 0.85,
                                       left = 0.08,
                                       right = 0.85 if len(self.ssl) == 2 else 0.92)

            # itération sur les slopes de plot
            for i, slope in enumerate(self.ssl):
                # avec mask, on masque toutes les classes qui ne sont pas dans la slope qu'on plote.
                fcl, mask = setSubsetclasses(self.sdObj.pgd, self.sel, self.sas, slope, )
                radii = np.ma.masked_array(self.rmax - self.sdObj.pgd.elev, mask = np.invert(mask)).compressed()
                theta = np.ma.masked_array(self.sdObj.pgd.aspect * np.pi / 180., mask = np.invert(mask)).compressed()
                # ici, on masque les nan dans la donnée à plotter
                ggg = np.ma.masked_invalid(self.sdObj.data[var])
                data = ggg.filled()
                validcolors = np.ma.masked_array(data, mask = np.invert(mask)).compressed()
                validcolors = np.ma.masked_where(validcolors == 1e+20, validcolors)
                # set corresponding fill_value
                cmap = plt.cm.get_cmap(cmap)
                cmap.set_bad(color = '0.5')
                # pie or imshow depending on the slope
                if slope == '0':
                    sax = plt.subplot(gs[0])
                    gg = sax.imshow(np.flipud(np.expand_dims(validcolors, 1)),
                                    interpolation  ='None', cmap = cmap, zorder = 0, vmin = lims[0], vmax = lims[1])
                    sax.autoscale(enable=False)
                    sax.set_yticks((self.sel - 600.) / 300.)
                    sax.set_yticklabels(reversed(list(map(str, list(map(int, self.sel))))))
                    sax.set_xticks([])
                    sax.set_xticklabels([])
                    if hasattr(self, 'focusCl'):
                        # BC 01/07/19 : set the following to false to print quantiles.
                        leg2ok = True
                        for iCl, cl in enumerate(self.focusCl):
                            if self.sdObj.pgd.aspect[cl] < -0.5:
                                if self.focusColors[iCl] > -0.5:
                                    _ = sax.scatter(1. / len(self.sel + 1.) - 1. / 20.,
                                                    (self.rmax - self.sdObj.pgd.elev[cl] - 300) / 300.,
                                                    c = self.focusColors[iCl], cmap = self.cmap_f,
                                                    vmin = 0., vmax = 1.,
                                                    zorder = 10, s=9)
                                else:
                                    _ = sax.scatter(1. / len(self.sel + 1.) - 1. / 20., (self.rmax - self.sdObj.pgd.elev[cl] - 300) / 300.,
                                                    c = 'r', zorder = 10, s=15)
                                if not leg2ok:
                                    xy = [0, (self.rmax - 600 - 300) / 300.]
                                    xybox = np.array([455., -31.])
                                    arr = np.ma.masked_invalid([[np.nan, 0., 1. / 3., 2. / 3., 1.]])
                                    im = OffsetImage(arr, zoom = 27, cmap = self.cmap_f, interpolation='None')
                                    im.image.axes = ax
                                    ab = AnnotationBbox(im, xy,
                                                        xybox=xybox,
                                                        xycoords='data',
                                                        boxcoords="offset points",
                                                        pad=0.03,
                                                        # arrowprops=dict(arrowstyle="->")
                                                        )
                                    sax.add_artist(ab)
                                    txtpos = [12.9, 11.8]
                                    sax.text(txtpos[0], txtpos[1],
                                             'NS Q1 Q2 Q3 Q4', fontsize = 12
                                             )
                                    leg2ok = True
                    sax.text(0, (self.rmax - 150 ) / 300., 'flat', horizontalalignment = 'center', fontsize = 15)
                else:  # 20 (or maybe 40)
                    sax = plt.subplot(gs[:, i], projection = 'polar')
                    sax.set_theta_zero_location('N')
                    sax.set_theta_direction(-1)
                    sax.set_yticklabels([])
                    bars = sax.bar(theta, radii, width = width, bottom=0.0, zorder = 0)

                    norm = Normalize(vmin = lims[0], vmax = lims[1])

                    for ib, b in enumerate(bars):
                        if np.ma.getmaskarray(validcolors)[ib]:
                            b.set_facecolor(np.squeeze(cmap(np.ma.masked_array([1], mask = [True]))))
                        else:
                            b.set_facecolor(cmap(norm(validcolors[ib])))
                    # if the constructor is provided a focus Cl, print x in it
                    if hasattr(self, 'focusCl'):
                        legOk = False
                        for iCl, cl in enumerate(self.focusCl):
                            if cl in fcl:
                                if self.focusColors[iCl] > -0.5:
                                    sax.scatter(self.sdObj.pgd.aspect[cl] * np.pi / 180.,
                                                self.rmax - self.sdObj.pgd.elev[cl] - 150.,
                                                c = str(self.focusColors[iCl]), cmap = self.cmap_f,
                                                vmin = 0., vmax = 1., zorder = 10,
                                                s=9, label = 'observed classes')
                                else:
                                    sax.scatter(self.sdObj.pgd.aspect[cl] * np.pi / 180., self.rmax - self.sdObj.pgd.elev[cl] - 150.,
                                                c = 'r', zorder = 10, s=15, label = 'observed classes')
                                if not legOk:
                                    # sax.legend(loc = (0.65, -0.13), scatterpoints = 1)
                                    legOk = True
                    # set the locations of the angular gridlines
                    _, _ = sax.set_thetagrids( np.arange(22.5, 360, 45) )
                    sax.set_xticklabels([])
                    sax.set_rticks(np.unique(radii))
                    # sax.grid(False)
                    [sax.text(th * np.pi / 180, 3700, tt,
                              horizontalalignment='center', verticalalignment='center',
                              fontsize = 8,)
                     for (th, tt) in zip(range(0, 360, 45), ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'])]
                    sax.text(-np.pi, 4450. if len(self.ssl) == 2 else 5600., '$' + slope + '^o$' + ' slope', horizontalalignment = 'center', fontsize = 15)

                    # lims = self.dictlims[var]
                    sax.set_ylim([np.min(radii) - 300., np.max(radii)])

                # colorbar
                # cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
                cax = plt.subplot(gs[i + 1])
                step = (lims[1] - lims[0]) / 10.
                cb = plt.colorbar(gg, cax=cax, ticks=np.arange(lims[0], lims[1] + step, step))
                cb.set_ticklabels(['{:.1f}'.format(x) for x in np.arange(lims[0], lims[1] + step, step)])
                # fig.subplots_adjust(left=0.05, right=0.85)
                plt.title(colortitle if colortitle is not None else niceLabel(var), fontsize = 10)
            if self.title == '':
                plt.suptitle(niceLabel(var, self.score, printunits=False), fontsize = 15)
            elif self.title is not None and self.maketitle is not False:
                plt.suptitle(self.title, fontsize = 15)
            if savefig:
                print('im here before saving', os.getcwd())
                if hasattr(self, 'focusCl'):
                    if len(self.focusCl) == 0:
                        fig.savefig('../pie/' + var + 'pie_' + self.sdObj.ptinom + '_' + '_' + str(int(self.sdObj.pgd.elev[self.focusCl])) +
                                    '_' + str(int(self.sdObj.pgd.aspect[self.focusCl])) +
                                    '_'      + str(int(np.arctan(self.sdObj.pgd.slope[self.focusCl]) * 180. / np.pi))  + '_' + self.sdObj.date + '.png')
                    else:
                        fig.savefig('../pie/' + var + 'pie_' + self.sdObj.ptinom + '.png', dpi=200)

                else:
                    fig.savefig('../pie/' + var + 'pie_' + self.sdObj.ptinom + '.png')
                # plt.close()
                plt.show()
            else:
                plt.show()


def prepare_spaghetti_Marie(run, var, cl, ax=None, color=None, density=False, label=None, score=None, title=True, legEl={}):

    if hasattr(run, 'ensProClim'):

        itimes, itimesC = set_itimes(run, clim=True, fromOl=True)
        itimesAn = set_itimes(run)
        xtimes = run.ensProOl['time'][itimes]
        clim = run.ensProClim[var][itimesC, cl, :]
    else:
        clim = None
        itimes = set_itimes(run, fromOl=True)
        itimesAn = set_itimes(run)
        xtimes = run.ensProOl['time'][itimes]

    ol = run.ensProOl[var][itimes, cl, :]
    truth = run.truth[var][itimes, :][:, cl]
    an = run.ensProAn[var][itimesAn, cl, :]
    # if assim with less member than openloop.
    if np.shape(an)[1] < np.shape(ol)[1]:
        ol = ol[:, 0:np.shape(an)[1]]
    legEl = spaghettiMarie(run, xtimes, ol, truth, an, cl, var,
                           clim=clim, ax=ax, color=color, density=density,
                           label=label, score=score, title=title, legEl=legEl)
    return legEl


def spaghettiMarie(run, xtimes, ol, truth, an, cl, var, clim=None,
                   ax=None, color=None, density=False, label=None, score=None, title=True, legEl={}):
    if ax is None:
        _, ax = plt.subplots(figsize=(14, 5))

    # plt.ylim([0., 700.])
    for tt in run.conf.assimdates:
        ax.axvline(x=tt, ymin=0, ymax=500, color='k', lw=0.5, linestyle='--', )
    if density is False:
        if score is not None:
            from utilpp import RMSE, spread
            anexp = np.expand_dims(an, axis=1)
            olexp = np.expand_dims(ol, axis = 1)
            truthexp = np.expand_dims(truth, axis=1)
            if ('sigma' in score) or ('spread' in score):
                ax.plot(xtimes, spread(
                    anexp), color='royalblue' if color is None else color, lw=2, alpha=0.7)
                ax.plot(xtimes, spread(
                    olexp), color='0.7', lw=2, alpha=0.7)
            else:
                ax.plot(xtimes, RMSE(anexp, truthexp), color='royalblue' if color is None else color, lw=2,
                        alpha=0.7)
                ax.plot(xtimes, RMSE(olexp, truthexp), color='0.7', lw=2,
                        alpha=0.7)

                # print(RMSE(anexp, truthexp))
            ax.set_ylabel(niceLabel(var, score=score))
        else:
            olp = ax.plot(xtimes, ol, color='0.9', alpha=0.5,
                          zorder=-1, label='openloop')
            if 'openloop' not in legEl.keys():
                legEl['openloop'] = ol
                print(np.shape(ol))
            ax.plot(xtimes, np.median(ol, axis=1), color='k')
            if clim is not None:
                ax.plot(xtimes, clim, color='brown', zorder=-1)
            ax.plot(
                xtimes, an, color='royalblue' if color is None else color, alpha=0.5)
            ax.plot(xtimes, np.mean(an, axis=1),
                    color='b' if color is None else color, lw=2, label=label)
            truthp = ax.plot(xtimes, truth, color='firebrick',
                             lw=2, label='truth')
            if 'truth' not in legEl.keys():
                legEl['truth'] = truthp
            ax.set_ylabel(niceLabel(var))
    else:
        olp = plot_density(xtimes, ol, ax, color='0.7',
                           zorder=-1, isOl=True, label='openloop')
        if 'openloop' not in legEl.keys():
            legEl['openloop'] = olp

        anp = plot_density(
            xtimes, an, ax, color='royalblue' if color is None else color, label=label)
        if label not in legEl.keys():
            legEl[label] = anp
        truthp = ax.plot(xtimes, truth, color='firebrick', lw=2, label='truth')
        if 'truth' not in legEl.keys():
            legEl['truth'] = truthp[0]
        ax.set_ylabel(niceLabel(var))
    if title is True:
        ax.set_title(niceName(run.pgd, cl))
    plt.gcf().autofmt_xdate(bottom=0.2)
    if not os.path.exists('../pie'):
        os.mkdir('../pie')
    plt.savefig('../pie/spaghetti_Marie_' + var + '_' + str(cl) + '.png')

    return legEl


def plot_density(xtimes, ens, ax, color='royablue', zorder=1, label=None, isOl=False):
    q = dict()
    if isOl:
        factalpha = 1
    else:
        factalpha = 1.
    for perc in range(0, 120, 20):
        q[perc] = np.percentile(ens, perc, axis=1)
    # ax.fill_between(xtimes,q[0],q[20], facecolor = color, alpha = factalpha *0.5, zorder = zorder,)
    # ax.fill_between(xtimes,q[80],q[100], facecolor = color, alpha = factalpha *0.5, zorder = zorder,)
    # ax.fill_between(xtimes,q[20],q[40], facecolor = color, alpha = factalpha *0.7, zorder = zorder,label = label, )
    # ax.fill_between(xtimes,q[60],q[80], facecolor = color, alpha = factalpha *0.7, zorder = zorder, )
    # ax.fill_between(xtimes,q[40],q[60], facecolor = color, alpha = factalpha *1, zorder = zorder,)
    plot = ax.fill_between(xtimes, q[0], q[100], facecolor=color,
                           alpha=factalpha * 0.7, zorder=zorder, label=label, )
    return plot


def spaghettispace(var, pb, pa, obs=None, savePath = None):
    plt.figure()
    gg = plt.plot(np.ma.masked_where(pb.stack[var] == 1e+20, pb.stack[var]), color = 'b', alpha = 0.5, label = 'bg')
    plt.setp(gg[1:], label="_")
    gg = plt.plot(np.ma.masked_where(pa.stack[var] == 1e+20, pa.stack[var]), color = 'g', alpha = 0.5, label = 'an')
    plt.setp(gg[1:], label="_")
    plt.ylabel(var)
    plt.xlabel('classes')
    if obs is not None:
        _ = plt.plot(np.ma.masked_where(obs.data[var] == 1e+20, obs.data[var]), color = 'm', linewidth =3, label = 'obs')
        plt.legend()
    if savePath is None:
        plt.show()
    else:
        plt.savefig( '../pie/' + savePath)
        plt.close()
    return 0


def cpdfs(options, pb, pa, obsReal, cl, alpha = None):
    '''
    pb : ensBg at a given date
    pb : ensAn at a given date
    obsReal : obs    ''
    alpha : inflation factor vector for this date (only for rlocal)
    '''
    sbplot = 100 + len(options.ppvars) * 10
    iv = 0
    plt.figure()
    if alpha is not None:
        errors = dictErrors()
    for var in options.ppvars:
        iv += 1
        plt.subplot(sbplot + iv)
        Xb = pb.stack[var][cl, :]
        Xb[np.isnan(Xb)] = 0.2
        Xb = np.sort(Xb)
        Fb = np.arange(len(Xb)) / float(len(Xb))
        Xa = pa.stack[var][cl, :]
        Xa[np.isnan(Xa)] = 0.2
        Xa = np.sort(Xa)
        Fa = np.arange(len(Xa)) / float(len(Xa))

        Xo = np.linspace(np.nanmin(np.array([Xa, Xb])), np.nanmax(np.array([Xa, Xb])), len(Xa))
        Xo = np.sort(np.append(Xo, obsReal.data[var][cl]))
        Fo = np.where(Xo < obsReal.data[var][cl], 0., 1.)
        lol = np.interp(Xo, Xb, Fb)
        integrand = (Fo - lol)**2
        if alpha is not None:
            plt.plot(Xo, 0.2 * np.exp(-0.5 * (Xo - obsReal.data[var][cl])**2 * alpha[cl] / errors[var]), color = 'r', label = 'like')
        plt.plot(Xb, Fb, color = 'b', linewidth = 2, label = 'bg')
        plt.plot(Xa, Fa, color = 'g', alpha = 0.5, linewidth = 2, label = 'an')
        plt.plot(Xo, Fo, color = 'm', linewidth = 2, label = 'obs')
        plt.plot(Xo, integrand, color = 'k')
        plt.legend()
        plt.title(var + ' CDFs in class ' + str(cl))
        from scipy import integrate
        CRPS_verif = integrate.simps(integrand, Xo)
    return CRPS_verif


def snsscatter_2class(options, pb, pa, obsReal, var, cl1, cl2, focusCl, savefig = False):

    p = sns.JointGrid(x = pb.stack[var][cl1, :], y = pb.stack[var][cl2, :])
    p = p.plot_joint(plt.scatter, s = 40, alpha = 0.5)
    if cl1 in focusCl and cl2 not in focusCl:
        plt.axvline(obsReal.data[var][cl1], color = 'm', linewidth = 3, alpha = 0.5)
        plt.scatter(obsReal.data[var][cl1], obsReal.data[var][cl2], color = 'm', s = 100, marker = '*')
    elif cl2 in focusCl and cl1 not in focusCl:
        plt.axhline(obsReal.data[var][cl2], color = 'm', linewidth = 3, alpha = 0.5)
        plt.scatter(obsReal.data[var][cl1], obsReal.data[var][cl2], color = 'm', s = 100, marker = '*')
    elif cl1 in focusCl and cl2 in focusCl:
        plt.scatter(obsReal.data[var][cl1], obsReal.data[var][cl2], color = 'm', s = 100, marker = '*')
    p.ax_marg_x.hist(pb.stack[var][cl1, :], alpha = 0.5, color ='b')
    p.ax_marg_y.hist(pb.stack[var][cl2, :], orientation = 'horizontal', alpha = 0.5)
    p.ax_marg_x.hist(pa.stack[var][cl1, :], alpha = 0.5, color ='g')
    p.ax_marg_y.hist(pa.stack[var][cl2, :], orientation = 'horizontal', alpha = 0.5, color = 'g')
    plt.scatter(pa.stack[var][cl1, :], pa.stack[var][cl2, :], color = 'g', alpha = 0.5, s=40)
    plt.xlabel('class ' + str(cl1))
    plt.ylabel('class ' + str(cl2))
    plt.suptitle(var)

    if savefig:
        plt.savefig('pie/snsscatter_' + var + '_' + str(cl1) + '_' + str(cl2) + '_' + str(options.dates) + '.png')
    plt.close()


def snsscatter_2vars(options, pb, pa, obsReal, var1, var2, cl, savefig = False):
    b1 = pb.stack[var1][cl, :]
    b2 = pb.stack[var2][cl, :]
    o1 = obsReal.data[var1][cl]
    o2 = obsReal.data[var2][cl]
    a1 = pa.stack[var1][cl, :]
    a2 = pa.stack[var2][cl, :]
    p = sns.JointGrid(x = b1, y = b2)
    p = p.plot_joint(plt.scatter, s = 40, alpha = 0.5)
    plt.scatter(o1, o2, color = 'm', s = 100, marker = '*')
    p.ax_marg_x.hist(b1, alpha = 0.5, color ='b')
    p.ax_marg_y.hist(np.ma.masked_invalid(b2).compressed(), orientation = 'horizontal', alpha = 0.5)
    p.ax_marg_x.hist(np.ma.masked_invalid(a1).compressed(), alpha = 0.5, color ='g')
    p.ax_marg_y.hist(np.ma.masked_invalid(a2).compressed(), orientation = 'horizontal', alpha = 0.5, color = 'g')
    plt.scatter(a1, a2, color = 'g', alpha = 0.5, s=40)
    plt.xlabel('class ' + str(cl) + ' ' + str(var1))
    plt.ylabel('class ' + str(cl) + ' ' + str(var2))
    plt.suptitle(cl)

    if savefig:
        plt.savefig('pie/snsscatter_' + var1 + '_' + var2 + '_' + str(cl) + '_' + str(options.dates) + '.png')
        plt.close()


def plot_pie_fromstack(run, kind, date, mb, listvar, ptinom='crpsa', focusCl=None, focusColors=None):
    sdObj = PrepAbs(date, run.options, ptinom=ptinom,
                    pgdPath=run.xpiddir + '/crampon/' + run.options.saverep + '/PGD.nc')
    data = dict()
    sdObj.options.ppvars = listvar
    for var in sdObj.options.ppvars:
        print(var)
        if kind == 'bg':
            data[var] = run.ensBg[date].stack[var][:, mb]
        elif kind == 'an':
            data[var] = run.ensAn[date].stack[var][:, mb]
    sdObj.data = data
    sdObj.isloaded = True
    sdObj.pgd = run.pgd
    sdObj.options.classesS = ['0', '20', '40']
    return sdObj, Pie(sdObj, focusCl=focusCl, focusColors=focusColors)


def plot_part_cesar_from_run(run, cl, varName, date, savefig = False, kindobs = 'Real'):
    if kindobs == 'Real':
        plot_part_cesar(run.ensBg[date], run.ensAn[date], run.obsReal[date], cl, varName, run.alpha[date], savefig = savefig)
    elif kindobs == 'Arch':
        plot_part_cesar(run.ensBg[date], run.ensAn[date], run.obsArch[date], cl, varName, run.alpha[date], savefig = savefig)

    else:
        plot_part_cesar(run.ensBg[date], run.ensAn[date], run.obsSynth[date], cl, varName, run.alpha[date], savefig = savefig)


def plot_part_cesar(pb, pa, obsReal, cl, varName, alpha = None, savefig = False):

    # version without color
    x_AN = 3
    x_BG = 1
    plt.figure()
    errors = dictErrors()
    yBG = np.ma.masked_invalid(pb.stack[varName][cl, :])
    yAN = np.ma.masked_invalid(pa.stack[varName][cl, :])
    yOBS = np.ma.masked_invalid(obsReal.data[varName][cl])
    # only length-1 for global:
    if not isinstance(alpha, float):
        alpha = alpha[cl]

    # f varName == 'SWE':
    # plt.fill_between((0,5),yOBS-errObs,yOBS+errObs,color='m',alpha=0.1)
    x = np.linspace(yOBS - 4 * np.sqrt(errors[varName]), yOBS + 4 * np.sqrt(errors[varName]), 100)
    y = np.exp(-0.5 * (x - yOBS)**2 / errors[varName])
    plt.plot(y + x_AN, x, color='m', linewidth=2, zorder=1000)

    if alpha is None or alpha != 1.:
        ya = np.exp(-0.5 * (x - yOBS)**2 * alpha / errors[varName])
        plt.plot(ya + x_AN, x, color='m', linewidth=2, zorder=1000, ls = 'dashed')
    a, binh = np.histogram(yBG.compressed(), bins=100)
    whatmb = 0
    for i, count in enumerate(a):
        for j in range(0, count):
            plt.plot(x_BG + j * 0.1, binh[i], 'o', color=(0.5, 0.5, 0.5), markersize=5, markeredgecolor='b')
            whatmb += 1

    a, binh = np.histogram(yAN.compressed(), bins=100)

    for i, count in enumerate(a):
        for j in range(0, count):
            plt.plot(x_AN + j * 0.1, binh[i], 'o', color=(0.7, 0.7, 0.7), markersize=5, markeredgecolor='g')

    plt.boxplot(yBG, positions=[x_BG - 0.2], showfliers=False, boxprops={'color': 'b'}, medianprops={'color': 'b'}, whiskerprops={'color': 'b'}, capprops={'color': 'b'})
    plt.boxplot(yAN, positions=[x_AN - 0.2], showfliers=False, boxprops={'color': 'g'}, medianprops={'color': 'g'}, whiskerprops={'color': 'g'}, capprops={'color': 'g'})

    plt.grid(True)
    plt.xlim((0, 6))

    plt.ylabel(varName)
    plt.xticks([])
    if savefig:
        plt.savefig('../pie/PartFilter.png')
