#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# import sys
# sys.path.append('/path/to/site-packages/')

import itertools as it
import numpy as np
from scipy import linalg as LA      # Removed in favor of scipy.sparse

"""
    __author__ = "Muammar El Khatib; Edoardo Fertitta."
    __copyright__ = "Copyright 2014, Muammar El Khatib; Edoardo Fertitta."
    __credits__ = [""]
    __license__ = "MIT"
    __version__ = ""
    __maintainer__ = "Muammar El Khatib"
    __email__ = "muammarelkhatib@gmail.com"
    __status__ = "Development"


    The MIT License (MIT)

    Copyright (c) 2015 Muammar El Khatib; Edoardo Fertitta.

    Permission is hereby granted, free of charge, to any person obtaining
    a copy of this software and associated documentation files (the
    "Software"), to deal in the Software without restriction, including without
    limitation the rights to use, copy, modify, merge, publish, distribute,
    sublicense, and/or sell copies of the Software, and to permit persons to
    whom the Software is furnished to do so, subject to the following
    conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.
"""

welcome = """\
===============================================================================


██╗  ██╗███████╗██╗███████╗███████╗███╗   ██╗██████╗ ███████╗██████╗  ██████╗
██║  ██║██╔════╝██║██╔════╝██╔════╝████╗  ██║██╔══██╗██╔════╝██╔══██╗██╔════╝
███████║█████╗  ██║███████╗█████╗  ██╔██╗ ██║██████╔╝█████╗  ██████╔╝██║  ███╗
██╔══██║██╔══╝  ██║╚════██║██╔══╝  ██║╚██╗██║██╔══██╗██╔══╝  ██╔══██╗██║   ██║
██║  ██║███████╗██║███████║███████╗██║ ╚████║██████╔╝███████╗██║  ██║╚██████╔╝
╚═╝  ╚═╝╚══════╝╚═╝╚══════╝╚══════╝╚═╝  ╚═══╝╚═════╝ ╚══════╝╚═╝  ╚═╝ ╚═════╝

A program to calculate the total position spread tensor using a Heisenberg
model hamiltonian.

This program has been written by Muammar El Khatib and Edoardo Fertitta.
===============================================================================
"""

print(welcome)

"""
When running python in clusters, sometimes you would need to specify the
site-packages path. If it is your case, just uncomment and configure the two
lines showed below accordingly.
"""


nelec = 4               # This is the basis state for the spin
j2 = 1.0                # This is the RATIO j2/j1 (low value high dimerization)
n = float(nelec)
bc = 'obc'              # obc and pbcfor open and periodic conditions respectively.
ms = 0                  # Set ms = 'all' in order to have all configurations.
nstate = 2              # Set nstate = 'all' in order to print all eigenvalues.


preket = list(it.product([0.5, -0.5], repeat=nelec))
prebra = list(it.product([1, 0], repeat=nelec))

"""
Ms partition
"""
bra = []

if isinstance(ms, str):
    bra = prebra
else:
    if ms != 0:
        divi = nelec / ms
        if ms == divi:
            for idx, i in enumerate(preket):
                if (np.absolute(np.sum(i)) == ms or
                   np.absolute(np.sum(i)) == -ms):
                    bra.append(prebra[idx])
                elif np.sum(i) == ms:
                    bra.append(prebra[idx])
        else:
            for idx, i in enumerate(preket):
                if np.sum(i) == ms:
                    bra.append(prebra[idx])
    else:
        for idx, i in enumerate(preket):
            if np.sum(i) == ms:
                bra.append(prebra[idx])

brams = []

for i in bra:
    brams.append(np.sum(i))

print('           Number of electrons: '+str(nelec))
print('              Number of states: '+str(nstate))
print('                            Ms: '+str(ms))
print('     Configurations for Ms = '+str(ms)+': '+str(len(bra)))
print('Total number of configurations: '+str(len(prebra)))
print('        Degree of dimerization: '+str(j2))

if bc == 'pbc':
    print('Periodic boundary conditions are used.')
else:
    print('Open boundary conditions are used.')
print('')

hmatrix = """\
-------------------------------------------------------------------------------
                        Heisenberg hamiltonian matrix
--------------------------------------------------------------------------------
"""
print(hmatrix)
print('Dimension of Ĥ matrix is: {}' .format(str(len(bra)) + 'x' +
      str(len(bra))))
print('')
shape = (len(bra), len(bra))
ovm = np.zeros(shape)

for idx, i in enumerate(bra):
    for idy, j in enumerate(bra):
        if brams[idx] == brams[idy]:
            diff = np.matrix(i) - np.matrix(j)
            vector_norm = np.square(np.linalg.norm(diff))

            if int(vector_norm) == 2:
                diffn = np.squeeze(np.asarray(diff))
                order = np.flatnonzero(diffn)

                if bc == 'pbc':
                    if (order[1] - order[0] == 1 or
                       order[1] - order[0] == nelec - 1):
                        ovm.itemset((idx, idy), 1)
                else:
                    if order[1] - order[0] == 1:
                        if order[0] % 2 == 0:
                            ovm.itemset((idx, idy), 1)
                        else:
                            ovm.itemset((idx, idy), j2)

"""
Define pairs
"""

if bc == 'obc':
    def pairs(lst):
        for i in range(1, len(lst)):
            yield lst[i - 1], lst[i]
else:
    def pairs(lst):
        n = len(lst)
        for i in range(n):
            yield lst[i], lst[(i + 1) % n]


prediagonal = []
expec = []

for idx, j in enumerate(bra):
    elemd = []
    suma = 0
    for i1, i2 in pairs(j):
        if suma % 2 == 0:
            if j2 < 0:
                r = j2
            else:
                r = 1.0
        else:
            r = j2

        if i1 == i2:
            elemd.append(.5 * r)
        else:
            elemd.append(-.5 * r)
        suma += 1
    prediagonal.append(elemd)


for xx in prediagonal:
    expec.append(np.sum(xx))

for idd, i in enumerate(prediagonal):
    j = np.sum(i)
    ovm.itemset((idd, idd), j)

print(ovm)

# Uncomment the line below to dump ovm matrix to text file.
# np.savetxt('ovm.txt', ovm, delimiter=',',fmt="%.5f")

hmatrixend = """\
-------------------------------------------------------------------------------
"""
print(hmatrixend)

print('')

e_vals, e_vecs = LA.eigh(ovm)


"""
This is let just for debugging reasons.
"""

np.set_printoptions(precision=16, suppress=True)

print('State vector configurations:')
print('')
print(bra)
print('')
if nstate == 'all':         # This is to print all possible eigenvalues
    nstate = len(e_vals)
for state in range(nstate):
    print('Results for STATE {}' .format(str(state + 1) + ': '))
    print('Eigenvalue = {}' .format(str(e_vals[state])))
    print('Corresponding state vector:')
    print(e_vecs[:, state])
    idxtomatch = np.flatnonzero(e_vecs[:, state])
    presumtps = []
    for i in idxtomatch:
        brarray = np.array(bra[i])

        if np.all(brarray == 0):
            brarrayn = brarray - 1
            positions = (np.flatnonzero(brarrayn)).astype(float) - \
                        (n - 1.) / 2.
        else:
            positions = (np.flatnonzero(brarray)).astype(float) - \
                         (n - 1.) / 2.
        sq = np.square(np.sum(positions))
        coefsq = np.square(e_vecs[:, state][i]) * sq
        presumtps.append(coefsq)
    print ('Λ(αα) contribution = ' + str(np.sum(presumtps)))
