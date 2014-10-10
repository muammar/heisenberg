#!/usr/bin/env python
# -*- coding: utf-8 -*-

import itertools as it
import numpy as np


# This is the basis state for the spin
nelec=4
j2=1.0   # This is the RATIO j2/j1 (low value high dimerization)
n=float(nelec)
bc='obc'
ms=0
nstate=2    #Number of states to print
preket=list(it.product([0.5,-0.5], repeat=nelec))
prebra=list(it.product([1,0], repeat=nelec))
#### Commented on 10/10/2014 print prebra
# Convert list of tuples to list of lists

welcome = """\
-------------------------------------------------------------------------------


██╗  ██╗███████╗██╗███████╗███████╗███╗   ██╗██████╗ ███████╗██████╗  ██████╗
██║  ██║██╔════╝██║██╔════╝██╔════╝████╗  ██║██╔══██╗██╔════╝██╔══██╗██╔════╝
███████║█████╗  ██║███████╗█████╗  ██╔██╗ ██║██████╔╝█████╗  ██████╔╝██║  ███╗
██╔══██║██╔══╝  ██║╚════██║██╔══╝  ██║╚██╗██║██╔══██╗██╔══╝  ██╔══██╗██║   ██║
██║  ██║███████╗██║███████║███████╗██║ ╚████║██████╔╝███████╗██║  ██║╚██████╔╝
╚═╝  ╚═╝╚══════╝╚═╝╚══════╝╚══════╝╚═╝  ╚═══╝╚═════╝ ╚══════╝╚═╝  ╚═╝ ╚═════╝

A program to calculate the total position spread tensor using a Heisenberg
hamiltonian.

This program has been written by Muammar El Khatib and Edoardo Fertitta.
-------------------------------------------------------------------------------
"""

print (welcome)

"""
Ms partition
"""
bra=[]
for idx,i in enumerate(preket):
    if np.sum(i) == ms:
        bra.append(prebra[idx])

#### Commented on 10/10/2014 print bra

brams=[]
for i in bra:
    brams.append(np.sum(i))


#### Commented on 10/10/2014 print brams

print ('Heisenberg hamiltonian matrix:')
print ('')
ovm = np.zeros(shape=(len(bra),len(bra)))
#### Commented on 10/10/2014 print bra
for idx,i in enumerate(bra):
    for idy,j in enumerate(bra):
        if brams[idx] == brams[idy]:
            diff=np.matrix(i)-np.matrix(j)
            vectornorm=np.linalg.norm(diff)**2
            if int(vectornorm) == 2:
#                print diff
                 diffn= np.squeeze(np.asarray(diff))
                 order=np.flatnonzero(diffn)
                 #### Commented on 10/10/2014 print order
                 if bc == 'pbc':
                     if order[1]-order[0] == 1 or  order[1]-order[0] == nelec-1: #PBC
                         ovm.itemset((idx,idy),1)
                 else:
                     if order[1]-order[0] == 1:
                         if order[0] % 2== 0:
                             ovm.itemset((idx,idy),1)
                         else:
                             ovm.itemset((idx,idy),j2)

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
    suma=0
    for i1, i2 in pairs(j):
#       print suma, i1, i2
        if suma % 2 == 0:
            r = 1.0
        else:
            r = j2
        if i1 == i2:
            elemd.append(0.5*r)
        else:
            elemd.append(-0.5*r)
        suma=suma+1
    prediagonal.append(elemd)
for xx in prediagonal:
    expec.append(np.sum(xx))
#print expec


for idd,i in enumerate(prediagonal):
    j=np.sum(i)
    #print i, idd, j
    #print type(i)
    ovm.itemset((idd,idd),j)

print ovm
print ('')

from scipy import linalg as LA
e_vals, e_vecs = LA.eigh(ovm)

"""
This is let just for debugging reasons.
"""
#for idx, i in enumerate(e_vals):
#    print e_vals[idx]
#    print e_vecs[:,idx]
#    print ''

enumerate(e_vals)

np.set_printoptions(precision=8,suppress=True)

for state in range(nstate):
    print ('Results for STATE '+str(state+1)+': ')
    print ('Eigenvalue= '+ str(e_vals[state]))
    print ('Corresponding state vector:')
    print e_vecs[:,state]
    idxtomatch=np.flatnonzero(e_vecs[:,state])
    presumtps=[]
    for i in idxtomatch:
        #   print bra[i]
        brarray=np.array(bra[i])
        positions=(np.flatnonzero(brarray)).astype(float)-(n-1.)/2.
        #   print positions
        sq=np.sum(positions)**2
        coefsq=e_vecs[:,state][i]**2*sq
        #   print coefsq
        presumtps.append(coefsq)
    #print presumtps
    print ('Λ(αα) contribution = ' + str(np.sum(presumtps)))
    print ('')
    #print ('Λ(αα+ββ) contrib = ' + str(np.sum(presumtps)))
