from SVM_bridge import *

os.chdir(binfile)

spec = []

for bfield in Brange:

    sout, tmp = outstr(
        outf='%s' % ('npS%d_' % spin + mpi + '_B%4.4f' % bfield + '.inp'),
        kk0=kk,
        mm0=mm,
        eB=bfield,
        eBspin=bfield,
        mnb=mn,
        bmin=bmi,
        bmax=bma,
        spincfg=[[1, 1, 2]],
        ispincfg=[[1, 1, 2]])

    outf = binfile + 'input/' + tmp

    with open(outf, 'w') as outfile:
        outfile.write(sout)

    os.system('./svm.x %s' % tmp[:-4])

    ofil = [line for line in open(binfile + 'output/' + tmp[:-4] + '.txt')]
    for n in range(len(ofil)):
        if re.search('more eigenvalues', ofil[-n]):
            spec += [[bfield, ofil[-n].split()[2:5]]]
            break

print(spec)