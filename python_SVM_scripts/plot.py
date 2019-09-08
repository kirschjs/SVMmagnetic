#!/usr/bin/python3
import sys
import matplotlib as plt
import matplotlib.pyplot as plt
import string
import os
from operator import attrgetter
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from scipy.interpolate import interp1d
#============================================================================   
list=["deuteron", "triton", "helium4"]
for name in list:  
 jobname=name   
 job = open(jobname+".sikum","r")
 E2=[]
 E4=[]
 E6=[]
 B2=[]
 B4=[]
 B6=[]
 for line in job:
   line = line.replace('=',' ')
   splt = line.split()
   magnetic = float(splt[splt.index('eB')+1])
   cutoff   = float(splt[splt.index('Cut')+1])
   energy   = float(splt[splt.index('E')+1])
   if(cutoff==2.0):
     B2.append(magnetic)
     E2.append(energy)
   if(cutoff==4.0):
     B4.append(magnetic)
     E4.append(energy)
   if(cutoff==6.0):
     B6.append(magnetic)
     E6.append(energy)
   #print("eB= "+str(magnetic)+"   Cut= "+str(cutoff)+"   E= "+str(energy))
 job.close()

 fig=plt.plot(B2,E2, label='$\Lambda=2$')
 fig=plt.plot(B4,E4, label='$\Lambda=4$')
 fig=plt.plot(B6,E6, label='$\Lambda=6$')
 if(jobname=="deuteron"):
   plt.axis([0, 1, -3, 25])
   plt.title("binding energy $E$ as function of magnetic field $B$ for deuteron")
   plt.legend()
   plt.xlabel('$B(fm^{-2})$')
   plt.ylabel('$E(MeV)$')
   plt.show()
   #plt.savefig('/home/moti/Dropbox/SVMmagnetic/deuteron.png')
   plt.close()
 if(jobname=="triton"):
   plt.axis([0, 1, -9, 25])
   plt.title("binding energy $E$ as function of magnetic field $B$ for triton")
   plt.legend()
   plt.xlabel('$B(fm^{-2})$')
   plt.ylabel('$E(MeV)$')
   plt.show()
   #plt.savefig('/home/moti/Dropbox/SVMmagnetic/triton.png')
   plt.close()
 if(jobname=="helium4"):
   plt.axis([0, 1, -25, 25])
   plt.title("binding energy $E$ as function of magnetic field $B$ for helium4")
   plt.legend()
   plt.xlabel('$B(fm^{-2})$')
   plt.ylabel('$E(MeV)$')
   plt.show()
   #plt.savefig('/home/moti/Dropbox/SVMmagnetic/helium4.png')
   plt.close()




