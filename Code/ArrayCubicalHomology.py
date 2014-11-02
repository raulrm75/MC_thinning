#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  ArrayCubicalHomology.py
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

import numpy as np
from itertools import product
import AbstractHomology as AH
import CubicalHomology as CH

# Class: Cube
class Cube(object):
	'''
	A Cube is made up of intervals.
	'''
	intervals = []
	emb = 0
	dim = -1
	cell = ()
	
	def __init__(self, cellList):
		if isinstance(cellList, list) or isinstance(cellList, tuple):
			self.emb = len(cellList)
			self.dim = sum(i % 2 for i in cellList)
			self.intervals = []
			for c in cellList:
				if c % 2 == 0:
					self.intervals.append((c/2, c/2))
				else:
					self.intervals.append(((c-1)/2, (c+1)/2))
			self.cell = tuple(cellList)
		else:
			raise Exception('Wrong type of argument for Cube contructor {}'.format(cellList))
	
	def drawToSVGCanvas(self, SVGCanvas, **kwargs):
		assert(self.emb <= 2)
		c = self.intervals
		if self.dim == 0:
			y = SVGCanvas.dotSize + 2 + c[0][0] * SVGCanvas.squareSize
			x = SVGCanvas.dotSize + 2 + c[1][0] * SVGCanvas.squareSize

			SVGCanvas.canvas.circle(x, y, SVGCanvas.dotSize, **kwargs)
		elif self.dim == 1:
			y = np.array(c[0]) * SVGCanvas.squareSize + SVGCanvas.dotSize + 2
			x = np.array(c[1]) * SVGCanvas.squareSize + SVGCanvas.dotSize + 2

			SVGCanvas.canvas.line(x[0], y[0], x[1], y[1], **kwargs)
		else:
			y = c[0][0] * SVGCanvas.squareSize + 2 * SVGCanvas.dotSize + 2
			x = c[1][0] * SVGCanvas.squareSize + 2 * SVGCanvas.dotSize + 2
			
			SVGCanvas.canvas.rect(x, y, 
				SVGCanvas.squareSize - 2 * SVGCanvas.dotSize, 
				SVGCanvas.squareSize - 2 * SVGCanvas.dotSize, 
				stroke_width=0, **kwargs)
		
	def border(self):
		n = self.emb
		cell = self.cell
		coeff = [1 if i == 0 else (-1)**sum(cell[j] % 2 for j in xrange(i))
			for i in xrange(n)]
		borderCells = {}
		for i in xrange(n):
			if cell[i] % 2 == 1:
				borderCells[tuple([cell[j] - 1 if j==i else cell[j] for j in xrange(n)])] = -coeff[i]
				borderCells[tuple([cell[j] + 1 if j==i else cell[j] for j in xrange(n)])] = coeff[i]
		
		return Chain(borderCells)
		
	def __lt__(self, other):
		if isinstance(other, Cube):
			if self.emb == other.emb:
				A = np.array(self.cell)
				B = np.array(other.cell)
				return (abs(A - B) <= np.ones(A.shape, A.dtype)).all() and (A != B).all()
			else:
				raise Exception('Only cubes with the same embedding number can be compared.')
		else:
			raise Exception('A cube can only be compared with another cube.')				
	
	def __le__(self, other):
		if isinstance(other, Cube):
			if self.emb == other.emb:
				A = np.array(self.cell)
				B = np.array(other.cell)
				return (abs(A - B) <= np.ones(A.shape, A.dtype)).all()
			else:
				raise Exception('Only cubes with the same embedding number can be compared.')
		else:
			raise Exception('A cube can only be compared with another cube.')				
	
	def __str__(self):
		return ' x '.join([
			'({})'.format(I[0]) if I[0] == I[1] else '({}, {})'.format(*I) 
			for I in self.intervals])
	
	def __eq__(self, other):
		if isinstance(other, Cube):
			if self.emb != other.emb or self.dim != other.dim:
				return False
			else:
				return self.cell == other.cell
	
	def __mul__(self, other):
		if isinstance(other, int):
			chain = Chain()
			if other != 0:
				chain.isZero = False
				chain.cubes.append(self.cell)
				chain.dim = self.dim
				chain.emb = self.emb
				chain.coeff[self.cell] = other
			return chain
		else:
		 raise Exception('Wrong type of operand {}'.format(other))
	
	def __add__(self, other):
		if isinstance(other, Cube):
			if self.emb == other.emb and self.dim == other.dim:
				return Chain(self) + Chain(other)
			else:
				raise Exception('Dimension missmatch.')
		elif isinstance(other, Chain):
			if self.emb == other.emb and self.dim == other.dim:
				return Chain(self) + other
			else:
				raise Exception('Dimension missmatch.')
		else:
			raise Exception('Wrong type of argument {}'.format(other))
	
	def __getitem__(self, index):
		if 0 <= index < self.emb:
			return self.cell[index]
		else:
			raise IndexError('Index error out of bounds {}'.format(index))
			
	def __sub__(self, other):
		return self + other * (-1)

	def incidence(self, cube):
		if isinstance(cube, Cube):
			if cube.dim != self.dim - 1:
				return 0
			else:
				bdr = self.border()
				if cube.cell in bdr.cubes:
					return bdr.coeff[cube.cell]
				else:
					return 0
		elif isinstance(cube, tuple) or isinstance(cube, list):
			return self.incidence(Cube(cube))
		else:
			raise Exception('Wrong type of argument {}'.format(cube))

## Class: Chain
class Chain(object):
	coeff = {}
	cubes = []
	isZero = True
	emb = 0
	dim = -1
	def __init__(self, data=None):
		self.coeff = {}
		self.cubes = []
		self.isZero = True
		self.emb = 0
		self.dim = -1
		if data:
			cube = None
			if isinstance(data, list) or isinstance(data, tuple):
				cube = Cube(data)
			elif isinstance(data, Cube):
				cube = data
			elif isinstance(data, int) and data == 0:
				pass
			elif isinstance(data, dict):
				for c in data:
					if data[c] != 0:
						self.coeff[c] = data[c]
				
				self.emb = Cube(data.keys()[0]).emb
				self.dim = Cube(data.keys()[0]).dim
				self.isZero = not bool(self.coeff)
				self.cubes = self.coeff.keys()
				self.cubes.sort()
			else:
				raise Exception('Wrong type of argument in Chain constructor {}'.format(cube))
			
			if cube:
				self.cubes.append(cube.cell)
				self.coeff[cube.cell] = +1
				self.isZero = False
				self.emb = cube.emb
				self.dim = cube.dim

	def __str__(self):
		if self.isZero:
			return '0'
		else:
			s = ''
			for c in self.cubes:
				coeff = abs(self.coeff[c])
				sgn = '+' if self.coeff[c] > 0 else '-'
				cube = str(Cube(c))
				if coeff == 1:
					s += '{} {}  '.format(sgn, cube)
				else:
					s += '{}{} · {}  '.format(sgn, coeff, cube)
			return s	
			
	def __add__(self, other):
		if isinstance(other, Chain):
			if other.isZero:
				return self
			elif self.isZero:
				return other
			else:
				if self.emb == other.emb and self.dim == other.dim:
					result = Chain()
					result.emb = self.emb
					result.dim = self.dim
					result.cubes = list(set(self.cubes + other.cubes))
					for c in result.cubes:
						coeff = 0
						if c in self.cubes:
							coeff += self.coeff[c]
						if c in other.cubes:
							coeff += other.coeff[c]
						result.coeff[c] = coeff
					cubesToRemove = []
					for c in result.cubes:
						if result.coeff[c] == 0:
							cubesToRemove.append(c)
					for c in cubesToRemove:
						result.cubes.remove(c)
						del result.coeff[c]
						
					if len(result.cubes) == 0:
						result.isZero = True
					else:
						result.isZero = False
					return result
				else:
					raise Exception('Two non zero chains must have the same embeding number and dimension to be operated.')
		elif isinstance(other, Cube):
			if other.emb == self.emb and other.dim == self.dim:
				if other.cell in self.cubes:
					self.coeff[other.cell] += 1
				else:
					self.cubes.append(other.cell)
					self.coeff[other.cell] = +1
			else:
				raise Exception('Two non zero chains must have the same embeding number and dimension to be operated.')
			
	def __sub__(self, other):
		if isinstance(other, Cube):
			return self + (Chain(other) * -1)
		elif isinstance(other, Chain):
			return self + (other * -1)
		else:
			raise Exception('Wrong type of operand {}'.format(other))
		
	def __mul__(self, other):
		if isinstance(other, int):
			if other != 0:
				chain = Chain(0)
				chain.isZero = self.isZero
				chain.emb = self.emb
				chain.dim = self.dim
				chain.cubes = self.cubes
				chain.coeff = self.coeff.copy()
				for c in chain.cubes:
					chain.coeff[c] = self.coeff[c] * other
				return chain
			else:
				return Chain(0)
		
	def __eq__(self, other):
		if isinstance(other, Chain):
			if self.dim != other.dim or self.emb != other.emb:
				return False
			else:
				if set(self.cubes) != set(other.cubes):
					return False
				else:
					eq = True
					for c in self.cubes:
						if self.coeff[c] != other.coeff[c]:
							eq = False
							break
					return eq
		elif isinstance(other, int) and other == 0:
			return self == Chain(0)
		else:
			return False
	
	def __ne__(self, other):
		return not (self == other)

	def __iter__(self):
		for cube in self.cubes:
			yield Cube(cube)
			
	def border(self):
		result = Chain(0)
		for cube in self.cubes:
			result += Cube(cube).border() * self.coeff[cube]
		return result
	
	def drawToSVGCanvas(self, canvas, **kwargs):
		assert(self.emb == 2)
		for cube in self:
			cube.drawToSVGCanvas(canvas, **kwargs)

## Class CubicalComplex
def shift(array, sign, direction):
	n = len(array.shape)
	result = np.zeros_like(array)
	idx = abs(np.sum(np.array(direction) * np.arange(n)))
	if sign > 0:
		S1 = tuple([
			slice(0, array.shape[i]) if i != idx else 
			slice(1, array.shape[i]) for i in xrange(n)])
		S2 = tuple([
			slice(0, array.shape[i]) if i != idx else 
			slice(0, array.shape[i] - 1) for i in xrange(n)])
	else:
		S1 = tuple([
			slice(0, array.shape[i]) if i != idx else 
			slice(0, array.shape[i] - 1) for i in xrange(n)])
		S2 = tuple([
			slice(0, array.shape[i]) if i != idx else 
			slice(1, array.shape[i]) for i in xrange(n)])
	result[S1] = array[S2]	
	
	return result

class CubicalComplex(object):
	emb = 0
	dim = -1
	cubes = None
	vectorField = None
	criticalCells = {}
	cellModels = {}
	dimensionMask = {}
	indexMask = None
	size = ()
	
	def __init__(self, size=None):
		if size != None:
			self.emb = len(size)
			self.dim = -1
			self.size = tuple(2*np.array(size) - 1)
			self.cubes = np.zeros(self.size, dtype=np.int8)
			prod = 1
			for s in self.size:
				prod *= s
			self.indexMask = np.arange(1, prod + 1).reshape(self.size)
			self.vectorField = np.zeros_like(self.indexMask)
			self.criticalCells = np.ones_like(self.cubes)
			
			cellList = list(product([0,1], repeat=self.emb))
			self.cellModels = {}
			for cell in cellList:
				d = sum(c % 2 for c in cell)
				if d in self.cellModels:
					self.cellModels[d].append(cell)
				else:
					self.cellModels[d] = [cell]
			for d in self.cellModels:
				self.cellModels[d].reverse()
				
			n = self.emb
			self.dimensionMask = {}
			for d in xrange(n+1):
				self.dimensionMask[d] = np.zeros_like(self.cubes)
				for cell in self.cellModels[d]:
					S = tuple([slice(0, self.size[i], 2) if cell[i] == 0 else slice(1, self.size[i] - 1, 2) for i in xrange(n)])
					self.dimensionMask[d][S] = 1
				
	def addCube(self, cube):
		newCube = None
		if isinstance(cube, Cube):
			newCube = cube
		elif isinstance(cube, tuple) or isinstance(cube, list):
			newCube = Cube(cube)
		else:
			raise Exception('Wrong type of argument for CubicalComplex contructor {}.'.format(cube))
		
		if newCube:
			if newCube.emb == self.emb:
				self.cubes[newCube.cell] = 1
				self.dim = max(self.dim, newCube.dim)
				for c in newCube.border().cubes:
					self.addCube(c)
			else:
				raise Exception('Wrong embeding number for cube {}'.format(newCube))
		
	def __iter__(self):
		for dim in xrange(self.emb + 1):
			M = self.cubes * self.dimensionMask[dim]
			for idx in product(*[xrange(s) for s in self.size]):
				if M[idx] > 0:
					yield Cube(idx)
				else:
					continue
			
	def __len__(self):
		return self.cubes.sum()
	
	def facets(self, cube):
		if isinstance(cube, list) or isinstance(cube, tuple):
			dim = Cube(cube).dim
			if dim >= self.emb:
				return []
			else:
				return self.facets(Cube(cube))
		elif isinstance(cube, Cube):
			if cube.dim >= self.emb:
				return []
			else:
				cell = cube.cell
				n = cube.emb
				evenIndices = [i for i in xrange(n) if cell[i] % 2 == 0]
				
				facetsCells = []
				for i in evenIndices:
					facetsCells.append(tuple([cell[j] if j != i else cell[i] - 1 for j in xrange(n)]))
					facetsCells.append(tuple([cell[j] if j != i else cell[i] + 1 for j in xrange(n)]))

				result = []
				c0 = np.zeros(n, np.int8)
				c1 = np.array(self.size)
				
				for cell in facetsCells:
					c = np.array(cell)
					validCell = (c0 <= c).all() and (c < c1).all()
					if validCell:
						result.append(cell)
				
				return [c for c in result if self.cubes[c] == 1]
		else:
			raise Exception('Wrong type of argument {}.'.format(cube))
	
	def __call__(self, idx):
		if self.indexMask.min() <= idx <= self.indexMask.max():
			return Cube(map(lambda x: x[0], np.where(self.cubes * self.indexMask == idx)))
		else:
			raise IndexError('Index {} out of bounds')
				
	def buildVectorField(self):
		C = self.cubes.copy()
		vF = np.zeros_like(self.indexMask)
		src = np.zeros_like(self.cubes)
		dst = np.zeros_like(self.cubes)
		for dim in xrange(self.emb):
			for direction in self.cellModels[1]:		
				K1 = shift(C * self.dimensionMask[dim], +1, direction)*self.dimensionMask[dim+1]*self.indexMask*C
				dst += (K1 > 0).astype(np.int8)
				K0 = shift(K1, -1, direction)*self.dimensionMask[dim]
				src += (K0 > 0).astype(np.int8)
				vF += K0
				C = self.cubes - src - dst
		
		self.vectorField = vF
		self.criticalCells = C
		
	def isCritical(self, cube):
		if isinstance(cube, list) or isinstance(cube, tuple):
			return self.isCritical(Cube(cube))
		elif isinstance(cube, Cube):
			return self.criticalCells[cube.cell] > 0
		else:
			raise Exception('Wrong type of argument {}'.format(cube))
	
	def V(self, chain):
		if isinstance(chain, Cube):
			if self.vectorField[chain.cell] == 0:
				return Chain(0)
			else:
				facet = self(self.vectorField[chain.cell])
				return Chain(facet) * facet.incidence(chain) * (-1)
				
		if isinstance(chain, Chain):
			result = Chain(0)
			for c in chain.cubes:
				if self.vectorField[c] > 0:
					result += Chain(self(self.vectorField[c])) * chain.coeff[c]
			return result
		else:
			raise Exception('Wrong type of argument {}'.format(chain))
	
	def d(self, chain):
		if isinstance(chain, Chain) or isinstance(chain, Cube):
			return chain.border()
		else:
			raise Exception('Wrong type of argument {}'.format(chain))
		
	def flow(self, cube):
		if isinstance(cube, list) or isinstance(cube, tuple):
			return self.flow(Cube(cube))
		elif isinstance(cube, Cube):
			return self.flow(Chain(cube))
		elif isinstance(cube, Chain):
			return cube + self.V(self.d(cube)) + self.d(self.V(cube))
		else:
			raise Exception('Wrong type of argument {}'.format(cube))
	
	def Flow(self, cube):
		prevFlow = Chain(0)
		currentFlow = self.flow(cube)
		while prevFlow != currentFlow:
			prevFlow = currentFlow
			currentFlow = self.flow(prevFlow)
		return currentFlow
		
	def criticalBorder(self, cube):
		if isinstance(cube, list) or isinstance(cube, tuple):
			return self.criticalBorder(Cube(cube))
		elif isinstance(cube, Cube):
			return self.Flow(cube)[0].border()
		else:
			raise Exception('Wrong type of argument {}'.format(cube))
						
	def translateChain(self, chain, cplxName='K'):
		if isinstance(chain, str):
			if chain == '0':
				return '0'
			chainInfo = map(lambda x: x.split(' · '), chain.split('  '))
			if chainInfo == [[chain]]:
				return '{}[{}]'.format(cplxName, self.index(Cube(chain)))
			elif len(chainInfo) == 1:
				return '0'
			else:
				chainInfo = chainInfo[:-1]
				if len(chainInfo[0]) == 2:
					for i in xrange(len(chainInfo)):
						chainInfo[i][1] = '{}[{}]'.format(cplxName, self.index(Cube(chainInfo[i][1])))
					result = map(lambda x: ' · '.join(x), chainInfo)
					return '  '.join(result)
				elif len(chainInfo[0]) == 1:
					for i in xrange(len(chainInfo)):
						sgn, cube = chainInfo[i][0][0], Cube(chainInfo[i][0][2:])
						chainInfo[i][0] = '{} {}[{}]'.format(sgn, cplxName, self.index(cube))
					return '  '.join(map(lambda x: x[0], chainInfo))
		elif isinstance(chain, Chain) or isinstance(chain, Cube):
			return self.translateChain(str(chain))
		else:
			raise Exception('Translate error with {}'.format(chain))
	
	def fromArray(self, array):
		self.__init__(array.shape)
		self.cubes = np.zeros(self.size, np.int8)
		K = np.zeros_like(self.cubes)
		n = self.emb
		K[tuple([slice(0,K.shape[i],2) for i in xrange(n)])] = array
		for d in xrange(1, n + 1):
			for c in self.cellModels[d]:
				K += (shift(K,+1,c) + shift(K,-1,c))*self.dimensionMask[d]/2
						
		self.cubes = K
	
	def index(self, cube):
		if isinstance(cube, list) or isinstance(cube, tuple):
			return self.index(Cube(cube))
		elif isinstance(cube, Cube):
			cell = cube.cell
			if (self.cubes * self.indexMask)[cell] > 0:
				return self.indexMask[cell]
			else:
				raise ValueError('{} is not in the cubical complex.'.format(cube))
		else:
			raise Exception('Wrong argument {}'.format(cube))
	
	def __getitem__(self, idx):
		if isinstance(idx, int):
			length = len(self)
			if idx < 0:
				return self[length + idx]
			elif idx >= length:
				raise IndexError('Wrong index {}'.format(idx))
			else:
				B = self.cubes * self.indexMask
				availableIndexes = list(B.flatten())
				nonzeroIndexes = [l for l in availableIndexes if l != 0]
				return Cube(map(lambda x: x[0], np.where(B == nonzeroIndexes[idx])))
		else:
			raise IndexError('Wrong type of index {}.'.format(idx))	
	
	def toCellComplex(self):
		cellCount = len(self)
		borderDict = {}
		dimDict = {}
		
		i = 0
		for idx in product(*[xrange(s) for s in self.size]):
			if self.cubes[idx] > 0:
				chain = Cube(idx).border()
				borderDict[i] = dict([(K.index(cube), chain.coeff[cube]) for cube in chain.cubes])
			dimDict[i] = K[i].dim
		
		return CellComplex(cellCount, borderDict, dimDict) 
	
	def toCubicComplexArray(self):
		L = []
		for c in self:
			idx = [2 * c[i][0] if c[i][0] == c[i][1] else 2*c[i][0] + 1 for i in xrange(self.emb)]
			L.append(idx)
		size = [1 + max([i[n] for i in L]) for n in xrange(self.emb)]
		A = np.zeros(size, np.int8)
		for idx in L:
			A[tuple(idx)] = 1
		return A
	
	def createSVGCanvas(self, fileName='temp.svg', squareSize=20, dotSize=3):
		return SVGCubical2DComplexCanvas(self, fileName, squareSize, dotSize)
		
	def drawToSVGCanvas(self, SVGCanvas, **kwargs):
		colorScheme = ['red', 'green', 'blue']
		for c in self:
			c.drawToSVGCanvas(SVGCanvas, fill=colorScheme[c.dim], **kwargs)
	
	def drawVectorToSVGCanvas(self, cube, SVGCanvas, **kwargs):
		if isinstance(cube, Cube):
			if cube.intervals in self.vectorField:
				p = cube.center + Cube(self.vectorField[cube.intervals]).center
				p = map(
					lambda x: 
						int(x * SVGCanvas.squareSize + SVGCanvas.dotSize + 2), 
					p)
					
				SVGCanvas.canvas.arrow(
					p[1], p[0], p[3], p[2], 
					SVGCanvas.arrowWidth, stroke_width=2, **kwargs) 
		elif isinstance(cube, list) or isinstance(cube, tuple):
			self.drawVectorToSVGCanvas(Cube(cube), SVGCanvas, **kwargs)
			
	def drawVectorFieldToSVGCanvas(self, SVGCanvas, **kwargs):
		for c in self.vectorField:
			self.drawVectorToSVGCanvas(c, SVGCanvas, fill='red', **kwargs)
	
	def criticalCubicalComplex(self):
		CCC = CubicalComplex(self.emb)
		criticalDims = self.criticalCells.keys()
		criticalDims.reverse()
		for dim in criticalDims:
			for cube in self.criticalCells[dim]:
				CCC.addCube(cube)
		return CCC
		
	def criticalComplex(self):
		cellCount = self.criticalCells.sum()
		criticalCells = [c for c in np.unique(self.criticalCells * self.indexMask) if c > 0]
		
		borderDict = {}
		dimDict = {}
		for i in criticalCells:
			criticalBorder = K.Flow(K(i)).border()
			borderDict[criticalCells.index(i)] = dict(
				[(criticalCells.index(K.index(Cube(c))),criticalBorder.coeff[c]) 
					for c in criticalBorder.cubes])
			dimDict[criticalCells.index(i)] = K(i).dim
			
		return AH.CellComplex(cellCount, borderDict, dimDict)
	
	def homology(self, index=None):
		C = self.criticalComplex()
		C.buildVectorField()
		return C.criticalComplex().homology(index)
	
	def toCubicalComplex(self):
		K = CH.CubicalComplex(self.emb)
		for cube in self:
			K.addCube(CH.Cube(cube.intervals))
		return K
	
if __name__ == '__main__':
	K = CubicalComplex((3,3,2))
	K.addCube((1,1,1))
	K.addCube((3,1,1))
	K.addCube((1,3,1))	
	K.buildVectorField()
	print K.homology()
