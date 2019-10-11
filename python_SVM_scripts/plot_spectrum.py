import sys, os, re
import matplotlib.pyplot as plt
import numpy as np

plt.style.use('bmh')

offs = 3
sysem = 'np_singlet_4_806'
sysem = 'np_3dosci_806'
infile = '/home/kirscher/kette_repo/SVMmagnetic/SVMxyzA/output/%s.txt' % sysem
data = [line for line in open(infile)]
evs = []
e0 = []
for n in range(1, len(data)):
    if re.search('more eigenvalues', data[n]):
        evs.append(data[n].split('=')[1].split())
        e0.append(float(data[n - 1].split('=')[2]))

datay = np.append(e0[-1], np.array([float(ev) for ev in evs[-1]]))

print(np.diff(datay), datay, e0[-1])

fig = plt.figure()
#num=None, figsize=(14, 24), dpi=60, facecolor='w', edgecolor='k')
#fig.suptitle(
#    r'$r_e\leq 2\left(R-\frac{R^2}{a}+\frac{R^3}{3a^2}\right)$', fontsize=16)
ax1 = fig.add_subplot(131)

ax1.plot(
    datay,
    'bo',
    label=r'$E_{(i)}(nn)$',
    linestyle='dashed',
    linewidth=1,
    alpha=0.5)
ax1.set_ylabel(r'$E_n\;\;[MeV]$', fontsize=15)

ax1 = fig.add_subplot(132)

ax1.plot(
    np.diff(datay),
    'bo',
    label=r'$\Delta E_{(i)}(nn)$',
    linestyle='dashed',
    linewidth=1,
    alpha=0.5)
ax1.set_ylabel(r'$E_{n-1}-E_n\;\;[MeV]$', fontsize=15)

ax1 = fig.add_subplot(133)

ax1.plot(
    datay / e0[-1],
    'bo',
    label=r'',
    linestyle='dashed',
    linewidth=1,
    alpha=0.5)
ax1.set_ylabel(r'$E_n/E_0$', fontsize=15)

#ax1.axhline(y=r, xmin=0, xmax=1)
#plt.ylim(-0.5 * r, 2 * r)
#plt.legend(loc='best', fontsize=22)
ax1.set_xlabel(r'$i$', fontsize=15)
plt.show()

fig.savefig("%s.pdf" % sysem, bbox_inches='tight')