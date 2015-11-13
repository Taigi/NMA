#!/usr/bin/env python
# Author: Asmit Bhowmick (asmit3@gmail.com)
# Initial code by Paul Glenn
# This python script will use nma.py to calculate nma similarity value of a trajectory
# Command line input. The input traj.pdb might depend on whether its generated using cpptraj. Please check before using with other amber versions 
# 1 = Start Model number
# 2 = End Model number
# usage: ./analyze_inms.py [start model] [end model]
# Note that start and end references have been hardcoded. Please change according to protein

highPass = 0.1 # percentage of top structures to admit
import cProfile
import re
import nma
import numpy as np
import sys
import time
import os

#prefix = '/global/scratch2/sd/asmit/cas9/apo'
#prefix2 = '/global/scratch2/sd/asmit/cas9/apo/config_0ns_hopper/'
prefix = os.getcwd()+'/'
#print prefix
#prefix = '/home/asmit/cas9/nma_codes/nma_python/nma_v2.3/'
#prefix = '/global/homes/a/asmit/nma_python/nma_v2.3/'
#prefix = '/Users/paulglen/Dropbox/CompBiophysics/1ptq.pdb'
#prefix = '/Users/paulglen/Dropbox/CompBiophysics/cas9/1aqb.pdb'
#prefix = '/Users/paulglen/github/NMA/Python/sandbox/1vom.pdb'
#prefix = '/Users/paulglen/Dropbox/CompBiophysics/cas9/inma_test/em/em.gro'
#prefix = '/Users/paulglen/Dropbox/CompBiophysics/cas9/inma_test/'

# CHANGE prefix for inms.dat file !!
#ref1 = prefix + '4CMP_50ns.pdb'
#ref2 = prefix + '4CMP_50ns.pdb'
ref1 = prefix + '4CMP_50ns.pdb'
ref2 = prefix + '4UN3_50ns.pdb'
trajFile =prefix + 'traj.pdb'
#trajFile = prefix2 + sys.argv[1] +'traj.pdb'
#em = prefix
#ref = prefix
#trajFile = prefix

startmod = int(sys.argv[1])
endmod = int(sys.argv[2])
tstart = time.clock()
pr = cProfile.Profile()
pr.enable()
X,Y,chains = nma.TrajectoryINMs(ref1,ref2, trajFile,startmod,endmod)
pr.disable()
pr.print_stats()

tdiff = time.clock() - tstart
print 'time taken: ', tdiff/ 60., 'min'

lenX = len(X)

#------------------------------------------------
# print all similarities
ofn = prefix+'inms_'+sys.argv[1]+'_'+sys.argv[2]+'.dat' # outfile name
print ofn
outfile = open(ofn,'w')
for i in range(0,len(chains)):
    outfile.write('%4i %6.3f %6.3f\n'%(chains[i],X[i],Y[i]))
#outfile.write('Checking inms.dat')
outfile.close()

#------------------------------------------------
# print top structures
#itr = 0
#ofn = prefix + 'top_structures_Ref1.dat'
#outfile = open(ofn, 'w')
#numStructures = int( highPass * lenX )
#while itr < numStructures and itr < lenX :
#    M = max(X)
#    frame = X.index(M)
#    outfile.write('%i %.2f\n'%(frame, X[frame]))
#    X.remove(M)
#    itr += 1
#outfile.close()
# ----
#itr = 0
#ofn = prefix + 'top_structures_Ref2.dat'
#outfile = open(ofn, 'w')
#numStructures = int( highPass * lenX )
#while itr < numStructures and itr < lenX :
#    M = max(Y)
#    frame = Y.index(M)
#    outfile.write('%i %.2f\n'%(frame, Y[frame]))
#    Y.remove(M)
#    itr += 1
#outfile.close()

tdiff = time.clock() - tstart
print 'time taken: ', tdiff/ 60., 'min'
