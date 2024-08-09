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


output_file_name = 'reference.dat'

#mol_h2o = pyscf.gto.M(atom = 'O 0 0 0; O 0.0 0.0 1.1 ')
#mol_h2o = pyscf.gto.M(atom = 'O 0 0 0;')
#mol_h2o = pyscf.gto.M(atom = 'H 0 2 0; H 0.5 0.5 0; H 2 2 0; H 3.5 0.5 0; H 0 0 0 ; H 1 1 1')
mol_h2o = pyscf.gto.M(atom = 'O 2 0 0; H 0 2 0; H 0.5 0.5 0; O 2 2 2; H 0 4 2; H 0.5 2.5 2')
#mol_h2o = pyscf.gto.M(atom = 'O 2 0 0; H 0 2 0; H 0.5 0.5 0')
#mol_h2o = pyscf.gto.M(atom = 'O 0 0 0; O 0.1 0.2 0.3; O -0.11 -0.33 -0.22; O 3 2 1 ')
#mol_h2o = pyscf.gto.M(atom = 'O 0 0 0; O 0 0 1; O 0 1 0; O 0 1 1; H 1 0 0; H 1 0 1; H 1 1 0; H 1 1 1;')
#mol_h2o = pyscf.gto.M(atom = 'O 0 0 0; O 0 0 1; O 0 1 0; O 0 1 1; O 1 0 0; O 1 0 1; O 1 1 0; O 1 1 1;')
#mol_h2o = pyscf.gto.M(atom=' O 2 0 0; H 2 0 1; H 2 1 0; O 2 1 1; H 2 0 0; H 2 0 1; O 2 1 0; H 2 1 1; H 3 0 0;')
#mol_h2o = pyscf.gto.M(atom = 'O 0 0 0;')

# basis = {atom_type1:[[angular_momentum
#                       (GTO-exp1, contract-coeff11, contract-coeff12),
#mol_h2o.basis = {'H': [ [0,[1.,1.]]]}
#mol_h2o.basis = {'H': [ [0,[1.,1.],[0.9,0.9]]]}
#mol_h2o.basis = {'O': [ [0,[1.,1.],[0.5,0.5]], [1,[1.,1.],[0.5,0.5]] ]}
#mol_h2o.basis = {'H': [ [1,[1.,1.],[2.,2.]] ]}
#mol_h2o.basis = {'H': [ [0,[1.,1.],[0.9,0.9]], [1,[1.,1.]] ]}
#mol_h2o.basis = {'H': [ [0,[1.,1.],[0.9,0.9]], [1,[1.,1.]], [2,[1.,1.]] ]}
#mol_h2o.basis = {'H': [ [0,[1.,1.],[0.9,0.9]], [1,[1.,1.]], [2,[1.,1.]], [3,[1.,1.]] ]}
#mol_h2o.basis = {'O': [[0,[1.,1.]]]}
#mol_h2o.basis = {'O': [[1,[1.,1.],[0.5,0.5]]]}
#mol_h2o.basis = {'O': [[1,[1.,1.]]]}
#mol_h2o.basis = {'O': [[2,[1.,1.]]]}
#mol_h2o.basis = {'O': [[0,[1.,1.]], [2,[1.,1.]]]}
#mol_h2o.basis = {'O': [[4,[1.,1.]], [2,[1.5,1.],[0.5,0.5]], [3,[1.,1.],[0.1,0.1]]],'H': [[1,[2.5,1.]], [0,[2.,1.],[0.1,1.0],[1.4,0.7]]] }
#mol_h2o.basis = {'O': [[2,[1.,1.],[0.5,0.5],[0.3,0.3]]], 'H': [[0,[1.,1.],[.5,.5]]] }
mol_h2o.basis = '631++g**'
#mol_h2o.basis = '631g'
#mol_h2o.basis = 'ccpvdz'
mol_h2o.build()

start_time = time.time()
ijkl = mol_h2o.intor("int2e")

# pyscf randomly fails and returns nan
while np.isnan(ijkl).any():
   ijkl = mol_h2o.intor("int2e")

end_time = time.time()

execution_time = end_time - start_time


with open(output_file_name, "w") as output_file_handle:
   output_file_handle.write("C\n")
   output_file_handle.write("1\n")
   nbas = len(mol_h2o._bas)

   ## Prints the bas, atm and env vectors
   np.savetxt(output_file_handle, np.array([nbas]), fmt="%i")
   np.savetxt(output_file_handle, mol_h2o._bas, fmt="%i")

   np.savetxt(output_file_handle, np.array([len(mol_h2o._atm)]), fmt="%i")
   np.savetxt(output_file_handle, mol_h2o._atm, fmt="%i")

   np.savetxt(output_file_handle, np.array([len(mol_h2o._env)]), fmt="%i")
   np.savetxt(output_file_handle, mol_h2o._env, fmt="%.30f")

   ## Solve the RHF problem and get the density
   mol_h2o_rhf = pyscf.scf.RHF(mol_h2o)
   mol_h2o_rhf.kernel()
   P_from_pyscf = mol_h2o_rhf.make_rdm1() 
   F_from_pyscf = - np.einsum('ij,ikjl',P_from_pyscf, ijkl)
   ## Computes the sparse transform of the density matrix
   SparseDens = np.empty(0)
   SparseFock = np.empty(0)
   offset_set = defaultdict(int)
   lda_set = defaultdict(int)
   bas = mol_h2o._bas
   istart = 0
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
               offset_set[set_i, set_j, atom_a, atom_b] = len( SparseDens )
               lda_set[set_i, set_j, atom_a, atom_b] = dj
               SparseDens = np.concatenate( (SparseDens, P_from_pyscf[istart:ifinish, jstart:jfinish].flatten()) )
               SparseFock = np.concatenate( (SparseFock, F_from_pyscf[istart:ifinish, jstart:jfinish].flatten()) )
           jstart += dj
       istart += di

   ## prints the structure of the density matrix
   np.savetxt(output_file_handle, np.array([len(offset_set)]), fmt="%i")
   for set_i, set_j, atom_a, atom_b in offset_set:
      off = offset_set[set_i, set_j, atom_a, atom_b]
      ld = lda_set[set_i, set_j, atom_a, atom_b]
      output_file_handle.write( " ".join( [ str(x) for x in [set_i, set_j, atom_a, atom_b, off, ld]] ) )
      output_file_handle.write( "\n" )

   ## Prints the actual sparse density matrix
   np.savetxt(output_file_handle, np.array([len(SparseDens)]), fmt="%i")
   np.savetxt(output_file_handle, SparseDens, fmt="%.30f")

   ## Prints the actual sparse fock matrix
   np.savetxt(output_file_handle, np.array([len(SparseFock)]), fmt="%i")
   np.savetxt(output_file_handle, SparseFock, fmt="%.30f")

#   ## Prints the dense density from pyscf
#   np.savetxt(output_file_handle, np.array([P_from_pyscf.shape[0]]), fmt="%i")
#   np.savetxt(output_file_handle, P_from_pyscf.flatten(), fmt="%.30f")

#   ## Prints the dense Fock from pyscf
#   np.savetxt(output_file_handle, F_from_pyscf.flatten(), fmt="%.30f")

   ## now write the integrals
   np.savetxt(output_file_handle, np.array([len(ijkl.flatten())]), fmt="%i")
   ia0 = 0
   for ibas in range(nbas):
      la = mol_h2o._bas[ibas][1]
      nla = mol_h2o._bas[ibas][3]
      na = (2*la+1) * nla;
      ib0 = 0
      for jbas in range(nbas):
         lb = mol_h2o._bas[jbas][1]
         nlb = mol_h2o._bas[jbas][3]
         nb = (2*lb+1) * nlb;
         ic0 = 0
         for kbas in range(nbas):
            lc = mol_h2o._bas[kbas][1]
            nlc = mol_h2o._bas[kbas][3]
            nc = (2*lc+1) * nlc;
            id0 = 0
            for lbas in range(nbas):
               ld = mol_h2o._bas[lbas][1]
               nld = mol_h2o._bas[lbas][3]
               nd = (2*ld+1) * nld;

               integers = ijkl[ ia0:ia0 + na, ib0:ib0 + nb, ic0:ic0 + nc, id0:id0 + nd ]            
               np.savetxt(output_file_handle, integers.flatten(), fmt="%.30f")

               id0 += nd
            ic0 += nc
         ib0 += nb
      ia0 += na








