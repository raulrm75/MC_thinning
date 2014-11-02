#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  sin título.py
#  
#  Copyright 2013 Raul Reina <raul@RMNet>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

class CellComplex(object):
	cellCount = 0
	_border = {}
	_skeleton = {}
	dim = {}
	vectorField = {}
	criticalCells = []
	maxDim = 0
	
	def __init__(self, cellCount, borderDict, dimDict):
		self.cellCount = cellCount
		self._border = borderDict
		self.dim = dimDict
		self.maxDim = max(self.dim.values())
		self._skeleton = dict([(n,[]) for n in xrange(self.maxDim + 1)])
		for i in xrange(self.cellCount):
			self._skeleton[self.dim[i]].append(i)
		
		for d in self._skeleton:
			self._skeleton[d].sort()
			
		self.vectorField = {}
		self.criticalCells = []
	
	def __call__(self, value):
		return self._skeleton[value]
	
	def __iter__(self):
		for d in self._skeleton:
			for i in self._skeleton[d]:
				yield i
	
	def __getitem__(self, cell):
		if 0 <= cell < self.cellCount:
			return Chain(self, cell)
		else:
			raise IndexError('Wrong index {}'.format(index))
			
	def border(self, i):
		result = Chain(self)
		for j in self._border[i]:
			result += self[j] * self._border[i][j]
		return result

	def facets(self, i):
		if self.dim[i] == self.maxDim:
			return []
		else:
			result = []
			for j in self:
				if i in self._border[j].keys():
					result.append(j)
			return result
	
	def buildVectorField(self):
		self.vectorField = {}
		for i in self:
			for j in self.facets(i):
				if self._border[j][i] == -1:
					if (i not in self.vectorField.keys() and
						i not in self.vectorField.values() and
						j not in self.vectorField.keys() and
						j not in self.vectorField.values()):
							self.vectorField[i] = j
							break
		self.buildCriticalCells()
	
	def buildCriticalCells(self):
		self.criticalCells = []
		for i in self:
			if i not in self.vectorField.keys() and i not in self.vectorField.values():
				self.criticalCells.append(i)
	
	def V(self, chain):
		if isinstance(chain, Chain):
			result = Chain(self)
			for i in chain.coeff:
				result += self.V(i) * chain.coeff[i]
			return result
		elif isinstance(chain, int):
			if chain in self.vectorField:
				return Chain(self, self.vectorField[chain])
			else:
				return Chain(self)
	
	def d(self, chain):
		if isinstance(chain, Chain):
			return chain.border()
		elif isinstance(chain, int):
			return self.d(Chain(self, chain))
	
	def flow(self, chain, repeat=1):
		if repeat == 1:
			if isinstance(chain, Chain):
				return chain + self.V(self.d(chain)) + self.d(self.V(chain))
			elif isinstance(chain, int):
				return self.flow(self[chain])
		else:
			return self.flow(chain, repeat-1)
	
	def Flow(self, chain):
		iterations = 1
		prevFlow = Chain(self)
		currentFlow = self.flow(chain)
		
		while prevFlow != currentFlow:
			prevFlow = currentFlow
			currentFlow = self.flow(prevFlow)
			iterations += 1
		
		return currentFlow

	def criticalComplex(self):
		cellCount = len(self.criticalCells)
		border = {}
		dim = {}
		for i in xrange(cellCount):
			dim[i] = self.dim[self.criticalCells[i]]
			border[i] = {}
			F = self.Flow(self.criticalCells[i]).border()
			for j in F.coeff:
				if F.coeff[j] != 0:
					border[i][self.criticalCells.index(j)] = F.coeff[j]

		return CellComplex(cellCount, border, dim)

	def homology(self, index=None):
		hasZeroBorder = True
		for i in self:
			hasZeroBorder = hasZeroBorder and self.border(i) == 0
		if hasZeroBorder:
			result = dict([(n,len(self(n))) for n in xrange(self.maxDim + 1)])
			if index is None:
				return result
			else:
				if index in result:
					return result[index]
				else:
					return 0
		else:
			raise Exception('Only complexes with null differential can be computed yet.')
		
		
class Chain(object):
	coeff = {}
	cellComplex = None
	dim = -1
	
	def __init__(self, cellComplex, cell=None, coeff=1):
		self.cellComplex = cellComplex
		self.coeff = {}
		if cell != None:
			if isinstance(cell, int):
				self.coeff[cell] = coeff
				self.dim = self.cellComplex.dim[cell]
			
	
	def __eq__(self, other):
		if isinstance(other, Chain):
			return self.coeff == other.coeff
		if isinstance(other, int) and other == 0:
			return self.isZero()
	
	def __ne__(self, other):
		return not self == other
	
	def isZero(self):
		iszero = True
		for i in self.coeff:
			iszero = iszero and self.coeff[i] == 0
			if not iszero:
				break
		return iszero
		
	def __add__(self, other):
		if isinstance(other, Chain):
			if self == 0 and other == 0:
				return Chain(self.cellComplex)
			elif self == 0 and other != 0:
				return other
			elif self != 0 and other == 0:
				return self
			else:
				if self.dim == other.dim:
					result = Chain(self.cellComplex)
					result.dim = self.dim
					S1 = set(self.coeff.keys())
					S2 = set(other.coeff.keys())
					L1 = [(i, self.coeff[i]) for i in S1-S2 if self.coeff[i] != 0]
					L2 = [(i, other.coeff[i]) for i in S2-S1 if other.coeff[i] != 0]
					L3 = [(i, self.coeff[i] + other.coeff[i]) for i in S1.intersection(S2) if self.coeff[i] + other.coeff[i] != 0]
					result.coeff.update(dict(L1))
					result.coeff.update(dict(L2))
					result.coeff.update(dict(L3))
					
					return result
				else:
					raise Exception('Two non zero chains must have the same dimension to be operated.')
		elif isinstance(other, int) and other == 0:
			return self
		else:
			raise Exception('Wrong argument {}'.format(other))
				
	def __mul__(self, other):
		if isinstance(other, int):
			if self.coeff and (other != 0):
				result = Chain(self.cellComplex)
				result.dim = self.dim
				for i in self.coeff:
					if self.coeff[i] != 0:
						result.coeff[i] = self.coeff[i] * other
				
				return result
			else:
				return Chain(self.cellComplex)
	
	def __sub__(self, other):
		return self + other * (-1)
	
	def __str__(self):
		result = ''
		for i in self.coeff:
			c = self.coeff[i]
			result += '{} {} · <{}> '.format(
				'+' if c > 0 else '-',
				abs(c),
				i)
		if not result:
			return '0'
		else:
			return result
		
	def border(self):
		if not self.coeff:
			return Chain(self)
		else:
			result = Chain(self)
			K = self.cellComplex
			L = [K.border(i) * self.coeff[i] for i in self.coeff]
			for c in L:
				result += c
			result.dim = self.dim - 1
			return result
			
# TODO:
# 1. Hay que comprobar que el campo vactorial sea acíclico.
# 2. Hay que ver como obtener los grupos de torsión a partir del opeador borde del complejo crítico.		


if __name__ == '__main__':
	KleinBottle = CellComplex(6, 
		{
			0:{}, 1:{}, 2:{}, 3:{},
			4:{1:1,2:1,3:-1},5:{1:1,2:-1,3:1}
		},
		{0:0,1:1,2:1,3:1,4:2,5:2})
	
	KleinBottle.vectorField[3] = 4
	for d in xrange(3):
		print 'd_{}'.format(d)
		for i in KleinBottle:
			if not i in KleinBottle.vectorField.keys() and not i in KleinBottle.vectorField.values() and KleinBottle.dim[i] == d:
				print '{} --> {}'.format(i, KleinBottle.Flow(i).border())
		
