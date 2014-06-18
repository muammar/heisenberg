#!/usr/bin/env python

import itertools as it
import numpy as np

###myarray = [1,1,1,1]
###
###combinations = it.combinations(myarray,2)
###
###pcombinations=[]
###for e in combinations:
###    pcombinations.append(e)
###
###print pcombinations


# This is the basis state for the spin
nelec=14
n=float(nelec)
bc='obc'
preket=list(it.product([1,-1], repeat=nelec))
bra=list(it.product([1,0], repeat=nelec))

# Convert list of tuples to list of lists

"""
Ms partition
"""

brams=[]
for i in bra:
    brams.append(np.sum(i))


print brams

ovm = np.zeros(shape=(len(bra),len(bra)))
print bra
for idx,i in enumerate(bra):
    for idy,j in enumerate(bra):
        if brams[idx] == brams[idy]:
            diff=np.matrix(i)-np.matrix(j)
            fuckit=np.linalg.norm(diff)**2
            if int(fuckit) == 2:
#                print diff
                 diffn= np.squeeze(np.asarray(diff))
                 order=np.flatnonzero(diffn)
#                print order
                 if bc == 'pbc':
                     if order[1]-order[0] == 1 or  order[1]-order[0] == nelec-1: #PBC
                         ovm.itemset((idx,idy),1)
                 else:
                     if order[1]-order[0] == 1:
                         ovm.itemset((idx,idy),1)

"""
Define pairs
"""
def pairs(lst):
    for i in range(1, len(lst)):
        yield lst[i-1], lst[i]
        if bc == 'pbc':
            yield lst[-1], lst[0]

prediagonal=[]
expec=[]
for idx,j in enumerate(bra):
    elemd=[]
    for i1, i2 in pairs(j):
#       print i1, i2
        if i1 == i2:
            elemd.append(0.5)
        else:
            elemd.append(-0.5)
    prediagonal.append(elemd)
for xx in prediagonal:
    expec.append(np.sum(xx))
#print expec


for idd,i in enumerate(prediagonal):
    j=np.sum(i)
    #print i, idd, j
    #print type(i)
    ovm.itemset((idd,idd),j)

#print ovm

from scipy import linalg as LA
e_vals, e_vecs = LA.eigh(ovm)

for idx, i in enumerate(e_vals):
    print e_vals[idx]
    print e_vecs[:,idx]
    print ''

gsev=e_vecs[:,0]

np.set_printoptions(precision=8,suppress=True)
print ('Print ground state vector')
print gsev

idxtomatch=np.flatnonzero(gsev)

presumtps=[]
for i in idxtomatch:
#   print bra[i]
    brarray=np.array(bra[i])
    positions=(np.flatnonzero(brarray)).astype(float)-(n-1.)/2.
#   print positions
    sq=np.sum(positions)**2
    coefsq=gsev[i]**2*sq
#   print coefsq
    presumtps.append(coefsq)

#print presumtps
print ('alfa,alfa contrib')
print ('')
print np.sum(presumtps)



##########
###########print hamiltonian
##########
##########For idx,j in enumerate(bra):
##########    print ('')
##########    print ('')
##########    print ('')
##########    print 'iter '+str(idx+1)
##########    for idy,k in enumerate(bra):
##########        subs=np.matrix(j)-np.matrix(k)
##########        print 'subs'
##########        print subs
##########        print ('producto punto')
##########        print np.dot(j,k)
##########        print ('indices')
##########        print idx,idy
###########       if np.dot(j,k) == 0 and np.sum(j) != nelec and np.sum(j) != -nelec:
##########        #hamiltonian.itemset((idx,idy),1)
##########
##########Print hamiltonian
##########
##########"""
##########Import numpy as np
###########import array as numpy array
##########Npcombinations=np.asarray(pcombinations)
##########Shape=npcombinations.shape
##########
##########Print shape
##########"""
