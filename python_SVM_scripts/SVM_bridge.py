import matplotlib.pyplot as plt
import sys, os, re
import numpy as np
from operator import itemgetter

binfile = '/home/johannesk/kette_repo/source/SVMmagnetic/SVMxyzA/'

MeVfm = 197.3161329
mn = {
    '137': 938.91852,
    '300': 1053.0,
    '450': 1226.0,
    '510': 1320.0,
    '805': 1634.0
}

mpi = '805'
spin = 1

spinconf = [[[1, 2, 1], [-1, 1, 2]], [[1, 1, 1]]][spin]

h2m = MeVfm**2 / mn[mpi]
Bi = 60 / h2m
Bf = 100 / h2m
Bn = 45

Brangelog = np.logspace(np.log10(Bi + 1e-12), np.log10(Bf), Bn)
Brangelin = np.linspace(Bi, Bf, Bn)

Brange = Brangelin

mm = 20
kk = 20
mn = 70

bmi = 0.1
bma = 6

#plt.plot(Brangelin, np.zeros(Bn) + 0.2, 'bo', alpha=0.5, markersize=2.0)
#plt.plot(Brangelog, np.zeros(Bn), 'ro', alpha=0.5, markersize=2.0)
#plt.show()


def outstr(outf='',
           npar=2,
           xm=1,
           spincfg=[],
           ispincfg=[],
           h2m=23.83,
           irand=-1776,
           ibf=1,
           eB=0.0,
           eBspin=0.0,
           mm0=15,
           kk0=5,
           mnb=45,
           bmin=0.01,
           bmax=7,
           npt=1,
           nop=3,
           vpot1=-396.568,
           vpot3=-8.0610,
           apot=4):
    sout = '\n'
    sout += 'npar=%d\n' % npar

    for p in range(npar):
        sout += 'xm(%1d)=1.0   ' % (p + 1)
    sout += '\n!\n'
    sout += 'nspc=%d\n' % len(spincfg)
    for scfg in range(len(spincfg)):
        sout += 'cspc(%d)=%f  ' % (scfg + 1, spincfg[scfg][0])
        for n in range(npar):
            sout += 'isp(%d,%d)=%d ' % (n + 1, scfg + 1, spincfg[scfg][n + 1])
        sout += '\n!\n'
    sout += 'nisc=%d\n' % len(ispincfg)
    for icfg in range(len(ispincfg)):
        sout += 'cisc(%d)=%f  ' % (icfg + 1, ispincfg[icfg][0])
        for n in range(npar):
            sout += 'iso(%d,%d)=%d ' % (n + 1, icfg + 1, ispincfg[icfg][n + 1])
        sout += '\n!\n'

    sout += 'h2m=%4.4f\n!\n' % h2m
    sout += 'irand=%d  ibf=%d\n!\n' % (irand, ibf)
    sout += 'eB=%8.8f  eBspin=%8.8f\n!\n' % (eB, eBspin)
    sout += 'mm0=%d  kk0=%d  mnb=%d\n' % (mm0, kk0, mnb)
    sout += 'bmin=%5.5f  bmax=%5.5f\n!\n' % (bmin, bmax)
    sout += 'npt=%d  nop=%d\n' % (npt, nop)
    sout += 'vpot(1,1)=%8.8f  apot(1,1)=%4.4f   bpot(1,1)=0  npot(1,1)=0\n' % (
        vpot1, apot)
    sout += 'vpot(1,3)=%8.8f  apot(1,3)=%4.4f   bpot(1,3)=0  npot(1,3)=0\n' % (
        vpot3, apot)

    if outf == '':
        print('no outputfile specified. output -> tmp.inp')

    return sout, outf