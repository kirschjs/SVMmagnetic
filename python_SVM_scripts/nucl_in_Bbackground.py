from SVM_bridge import *

os.chdir(binfile)

bfield = 20.0 / h2m

sout, tmp = outstr(
    outf='tmp_pp.inp',
    #vpot1=0.0,
    #vpot3=0.0,
    kk0=25,
    mm0=12,
    eB=bfield,
    eBspin=0.0,
    mnb=70,
    bmin=0.3,
    bmax=8,
    ispincfg=[[1, 1, 2], [-1, 2, 1]],
    spincfg=[[1, 1, 1]])

outf = binfile + 'input/' + tmp

with open(outf, 'w') as outfile:
    outfile.write(sout)

os.system('%s/svm.x %s' % (binfile, tmp[:-4]))

spec = []

ofil = [line for line in open(binfile + 'output/' + tmp[:-4] + '.txt')]
for n in range(len(ofil)):
    if re.search('more eigenvalues', ofil[-n]):
        spec += [[bfield, ofil[-n].split()[2:5]]]
        break

print(spec)