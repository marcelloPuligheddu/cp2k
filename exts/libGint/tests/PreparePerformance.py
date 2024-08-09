"""
Copyright (c) 2023 Science and Technology Facilities Council

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import numpy as np
import math
import pyscf
from collections import defaultdict
import time
import sys


def NLco(L):
    return (L+1)*(L+2)//2

def NLsp(L):
    return (2*L+1)


output_file_name = 'performance_test.inp'

#mol_h2o = pyscf.gto.M(atom = 'O 0 0 0; O 0.0 0.0 1.1 ')
#mol_h2o = pyscf.gto.M(atom = 'O 0 2 0;')
#mol_h2o = pyscf.gto.M(atom = 'H 0 2 0; H 0.5 0.5 0; H 2 2 0; H 3.5 0.5 0; H 0 0 0 ; H 1 1 1')
#mol_h2o = pyscf.gto.M(atom = 'O 2 0 0; H 0 2 0; H 0.5 0.5 0; O 2 2 2; H 0 4 2; H 0.5 2.5 2')
mol_h2o = pyscf.gto.M(atom = """
   O 1 0 0.53; H 1 0.1 0.52; H 1 0.9 0.54; O 2 0 0.64; H 2 0.2 0.21; H 2 0.7 0.41;
   O 3 0 0.26; H 3 0.2 0.74; H 3 0.5 0.74; O 4 0 0.81; H 4 0.3 0.47; H 4 0.3 0.13;
   O 4 0 0.26; H 4 0.2 0.74; H 4 0.5 0.74; O 5 0 0.81; H 5 0.3 0.47; H 5 0.3 0.13
""")
#mol_h2o = pyscf.gto.M(atom = 'O 2 0 0; H 0 2 0; H 0.5 0.5 0')
#mol_h2o = pyscf.gto.M(atom = ' Ti 0 0 0; O 0 2 0; O 0.5 0.5 0; Ti 1 0 0; O 1.2 2 0; O 1.25 0.5 0;')
#mol_h2o = pyscf.gto.M(atom = 'O 0 0 0; O 0.1 0.2 0.3; O -0.11 -0.33 -0.22; O 3 2 1; O 3. 3. 3. ')
#mol_h2o = pyscf.gto.M(atom = 'O 0 0 0; O 0.1 0.2 0.3; O -0.11 -0.33 -0.22; O 3 2 1; O 3. 3. 3.; O 1.23 4.21 0.13; O 2. 4. 0. ')
#mol_h2o = pyscf.gto.M(atom = 'O 0 0 0; O 0 0 1; O 0 1 0; O 0 1 1; H 1 0 0; H 1 0 1; H 1 1 0; H 1 1 1;')
#mol_h2o = pyscf.gto.M(atom = 'O 0 0 0; O 0 0 1; O 0 1 0; O 0 1 1; O 1 0 0; O 1 0 1; O 1 1 0; O 1 1 1;')
#mol_h2o = pyscf.gto.M(atom=' O 2 0 0; H 2 0 1; H 2 1 0; O 2 1 1; H 2 0 0; H 2 0 1; O 2 1 0; H 2 1 1; H 3 0 0;')
#mol_h2o = pyscf.gto.M(atom = 'O 0 0 0;')

# basis = {atom_type1:[[angular_momentum
#                       (GTO-exp1, contract-coeff11, contract-coeff12),
#mol_h2o.basis = {'H': [ [0,[1.,1.]]]}
#mol_h2o.basis = {'H': [ [0,[1.,1.],[0.9,0.9]]]}
#mol_h2o.basis = {'O': [ [1,[1.,1.],[0.9,0.9]],  [2,[1.,1.],[0.5,0.5]], [3,[1.,1.],[0.5,0.5]] ]}
#mol_h2o.basis = {'O': [ [1,[1.,1.]] ]}
#mol_h2o.basis = {'O': [ [1,[1.,1.],[2.,2.]] ]}
#mol_h2o.basis = {'O': [ [1,[1.,1.],[2.,2.],[3.,3.]] ]}
#mol_h2o.basis = {'H': [ [0,[1.,1.],[0.9,0.9]], [1,[1.,1.]] ]}
#mol_h2o.basis = {'H': [ [0,[1.,1.],[0.9,0.9]], [1,[1.,1.]], [2,[1.,1.]] ]}
#mol_h2o.basis = {'H': [ [0,[1.,1.],[0.9,0.9]], [1,[1.,1.]], [2,[1.,1.]], [3,[1.,1.]] ]}
#mol_h2o.basis = {'O': [[3,[1.,1.]]]}
#mol_h2o.basis = {'O': [[1,[1.,1.]]]}
#mol_h2o.basis = {'O': [[1,[1.,1.]]]}
#mol_h2o.basis = {'O': [[0,[1.,1.]], [1,[1.,1.]]]}
#mol_h2o.basis = {'O': [[0,[1.,1.]], [1,[1.,1.],[2.,2.]]]}
#mol_h2o.basis = {'O': [[4,[1.,1.]], [2,[1.5,1.],[0.5,0.5]], [3,[1.,1.],[0.1,0.1]]],'H': [[1,[2.5,1.]], [0,[2.,1.],[0.1,1.0],[1.4,0.7]]] }
#mol_h2o.basis = {'O': [[2,[1.,1.],[0.5,0.5],[0.3,0.3]]], 'H': [[0,[1.,1.],[.5,.5]]] }
mol_h2o.basis = '631++g**'
#mol_h2o.basis = '631g'
#mol_h2o.basis = 'ccpvdz'
mol_h2o.build()

#start_time = time.time()
#ijkl = mol_h2o.intor("int2e")
#
# pyscf randomly fails and returns nan
#while np.isnan(ijkl).any():
#   print('pyscf failed. retrying. Ignore timer ! ')
#   ijkl = mol_h2o.intor("int2e")
#
#end_time = time.time()
#execution_time = end_time - start_time
#print( 'Pyscf took ', execution_time*1.e3, "milliseconds" )

with open(output_file_name, "w") as output_file_handle:
   output_file_handle.write("P\n")
   nbas = len(mol_h2o._bas)
   np.savetxt(output_file_handle, np.array([nbas]), fmt="%i")
   np.savetxt(output_file_handle, mol_h2o._bas, fmt="%i")

   np.savetxt(output_file_handle, np.array([len(mol_h2o._atm)]), fmt="%i")
   np.savetxt(output_file_handle, mol_h2o._atm, fmt="%i")

   np.savetxt(output_file_handle, np.array([len(mol_h2o._env)]), fmt="%i")
   np.savetxt(output_file_handle, mol_h2o._env, fmt="%.30f")

   offset_set = defaultdict(int)
   lda_set = defaultdict(int)
   bas = mol_h2o._bas
   istart = 0
   curr_offset = 0

   for set_i in range(len(bas)):
       atom_a = bas[set_i][0]
       li = bas[set_i][1]
       nli = bas[set_i][3]
       di = (2*li+1) * nli
       ifinish = istart + di
       jstart = 0
       for set_j in range(len(bas)):
           atom_b = bas[set_j][0]
           lj = bas[set_j][1]
           nlj = bas[set_j][3]
           dj = (2*lj+1) * nlj
           jfinish = jstart + dj
           if atom_a >= atom_b:
               offset_set[set_i, set_j, atom_a, atom_b] = curr_offset
               lda_set[set_i, set_j, atom_a, atom_b] = dj
               curr_offset += di*dj
           jstart += dj
       istart += di

   ## prints the structure of the density matrix
   np.savetxt(output_file_handle, np.array([len(offset_set)]), fmt="%i")
   for set_i, set_j, atom_a, atom_b in offset_set:
      off = offset_set[set_i, set_j, atom_a, atom_b]
      ld = lda_set[set_i, set_j, atom_a, atom_b]
      output_file_handle.write( " ".join( [ str(x) for x in [set_i, set_j, atom_a, atom_b, off, ld]] ) )
      output_file_handle.write( "\n" )

   np.savetxt(output_file_handle, np.array([curr_offset]), fmt="%i")

   np.savetxt(output_file_handle, np.array([len(bas)**4]), fmt="%i")




