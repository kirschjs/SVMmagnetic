from SVM_bridge import *

rundir = binfile + 'output/'

basestr = 'npS%d_%s_B' % (spin, mpi)

data = []
for filen in os.listdir(rundir):
    for bfield in Brange:
        if re.search(basestr + '%4.4f' % bfield, filen):
            lines = [line for line in open(rundir + filen)]
            for n in range(len(lines)):
                if re.search('eB=', lines[n]):
                    xdata = float(lines[n].split('=')[1].split()[0])
                if re.search('h2m=', lines[n]):
                    h2m = float(lines[n].split('=')[1])
            for n in range(len(lines)):
                if re.search('more', lines[-n]):
                    ydata = np.array(
                        [lines[-n].split('=')[1].split()[:12]]).astype(float)
                    break
            data += [[xdata, ydata[0]]]

data = sorted(data, key=itemgetter(0))
xdata = [b[0] * h2m for b in data]

fig = plt.figure()
#num=None, figsize=(14, 24), dpi=60, facecolor='w', edgecolor='k')
fig.suptitle(
    r'neutron-proton system ($m_\pi=%dMeV\;\;\;\xi_0=\vert \sigma_p^z=-\sigma_n^z\rangle$)'
    % (int(mpi)),
    fontsize=16)
ax1 = fig.add_subplot(111)

[
    ax1.plot(
        xdata, [b[1][n] for b in data],
        label=r'$E_{(%d)}$' % n,
        linestyle='dashed',
        linewidth=1,
        alpha=0.5) for n in range(min(15, len(data[0][1])))
]

ax1.set_xlabel(
    r'$\left\vert B_z\right\vert\;\;\left[\frac{m}{\hbar^2}\right]$',
    fontsize=15)
ax1.set_ylabel(r'$E_n\;\;[MeV]$', fontsize=15)

#ax1.axhline(y=r, xmin=0, xmax=1)
#plt.ylim(-0.5 * r, 2 * r)
plt.legend(loc='best', fontsize=18)
plt.show()

fig.savefig("diproton_eB.pdf", bbox_inches='tight')