#!/usr/bin/python3
import sys
import string
import os
import matplotlib as plt
from operator import attrgetter
#import matplotlib.pyplot as plt
#from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
#from scipy.optimize    import curve_fit
#import lmfit
#from lmfit import minimize, Parameters

#=============================================================================
class jobclass:
        def __init__(self, dirname,jobname):
                self.dirname = dirname 
                self.name    = jobname
                self.filename = dirname+'/'+jobname
                self.npar    = 0
                self.cutoff  = 0.0
                self.eB      = 0.0
                self.iters   = 0
                self.ee0     = 0.0
                self.ee1     = 0.0
                self.dee0    = 0.0
        def __repr__(self):
                return repr(self.name)
#=============================================================================
def GetObservables(job):
        
    filename=job.filename
    jobfile = ""
    for line in open(filename):
        if (line.find('!')!=-1):
           line = line[:line.index('!')]     
        jobfile += line
    jobfile = jobfile.replace('=',' ')
    jobsplt = jobfile.split()
    job.npar    = int(  jobsplt[jobsplt.index('npar')+1])
    apot        = float(jobsplt[jobsplt.index('apot(1,1)')+1])
    job.cutoff  = 2*np.sqrt(apot)
    job.C1      = float(jobsplt[jobsplt.index('vpot(1,1)')+1])
    job.C3      = float(jobsplt[jobsplt.index('vpot(1,3)')+1])
    job.eB      = float(jobsplt[jobsplt.index('eB')+1])
    job.bmax      = float(jobsplt[jobsplt.index('bmax')+1])
#    job.D       = float(jobsplt[jobsplt.index('vpot3b')+1])
    iters       = [ jobsplt[i+1] for i in range(0,len(jobsplt)) if (jobsplt[i] == 'itr')]
    eners       = [ jobsplt[i+3] for i in range(0,len(jobsplt)) if (jobsplt[i] == 'itr')]
    job.iters   = int(  iters[-1])
    job.ee0     = float(eners[-1])
    print(job.name,job.npar,job.cutoff,job.ee0)

#=============================================================================
def CollectData(mydir,jobname):
    global jobs
    joblist = []
    for filename in os.listdir(mydir):
        fsize=os.path.getsize(mydir+"/"+filename)
        if (fsize == 0): continue
        if (not filename.endswith(".txt")): continue
        if (filename.find(jobname)==-1):   continue
        print(filename)
        job = jobclass(mydir,filename)
        joblist.append(job)
        GetObservables(job)

    jobs = sorted(joblist,key=attrgetter('npar','cutoff', 'eB'))
#=============================================================================
def WriteSikum(jobname):
    sikum = open(jobname+".sikum","w")
    for job in jobs:
        CM=0.5*41.47*job.eB
        fmt = "{:<32s}     Cut={:3.2f}      bmax={:3.2f}       eB={:3.6f}       Iters={:4d}       E_R={:3.6f}       E_CM={:3.6f}       E={:3.6f} "
        line = fmt.format(job.name ,job.cutoff, job.bmax, job.eB,  job.iters,  job.ee0,  CM,  job.ee0+CM)
        print(line)
        sikum.write(line+"\n")
    sikum.close()        
#=============================================================================

mydir='./input/'
jobname='deuteron'
CollectData(mydir,jobname)
WriteSikum(jobname)
jobname='triton'
CollectData(mydir,jobname)
WriteSikum(jobname)
jobname='helium4'
CollectData(mydir,jobname)
WriteSikum(jobname)

