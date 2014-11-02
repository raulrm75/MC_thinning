# -*- coding: utf-8 -*-
#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Cubicalomology.py
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
from itertools import product
import numpy as np
import AbstractHomology as AH
import CubicComplex as CC
from Language import Language
from Membranes import Membrane
from SVGCubicalComplex2DCanvas import SVGCubical2DComplexCanvas

class Interval(object):
	'''
	An Interval is made up of two integers a, b.
	If b = a the interval is degenerated.
	If b = a+1 the interval is non degenerated.
	'''
	
	def __init__(self, a, b=None):
		if b == None and isinstance(a, int):
			self.interval = (a, a)
			self.isDegenerated = True
		elif isinstance(b, int) and isinstance(a, int):
			if a == b:
				self.interval = (a, a)
				self.isDegenerated = True
			elif a == b - 1:
				self.interval = (a, b)
				self.isDegenerated = False
			else:
				raise Exception('Wrong arguments for Interval constructor: {}.'.format((a,b)))
		else:
			raise Exception('Wrong arguments for Interval constructor: {}.'.format((a,b)))
		
		self.dim = self.interval[1] - self.interval[0]
	
	def __str__(self):
		if self.isDegenerated:
			return '({})'.format(self.interval[0])
		else:
			return '({}, {})'.format(*self.interval)
	
	def __mul__(self, other):
		if isinstance(other, Interval):
			return Cube([self.interval, other.interval])
		elif isinstance(other, int):
			chain = Chain()
			if other != 0:
				chain.cubes.append((self.interval,))
				chain.coeff[(self.interval,)] = other
				chain.isZero = False
				chain.emb = self.emb
				chain.dim = self.dim
			return chain
		elif isinstance(other, Cube):
			return Cube([self.interval] + other.intervals)
		elif isinstance(other, Chain):
			chain = Chain()
			if not other.isZero:
				for c in other.cubes:
					cube = self * Cube(c)
					chain.cubes.append(cube.intervals)
					chain.coeff[cube.intervals] = other.coeff[c]
				chain.isZero = False
				chain.emb = other.emb + 1
				chain.dim = other.dim + self.dim
			return chain
		else:
			raise Exception('Wrong type of operand {}.'.format(other))
	
class Cube(object):
	'''
	A Cube is made up of intervals.
	'''
	intervals = []
	emb = 0
	dim = -1
	center = ()
	
	def __init__(self, intervals):
		if isinstance(intervals, list) or isinstance(intervals, tuple):
			self.intervals = tuple([Interval(*I).interval for I in intervals])
			self.emb = len(self.intervals)
			self.dim = sum(I[1] - I[0] for I in self.intervals)
			self.center = tuple([(I[0] + I[1])/2. for I in self.intervals])
		elif isinstance(intervals, str):
			intervalTuple = map(
				lambda x: x if isinstance(x, tuple) else (x,x), 
				map(eval, intervals.split(' x '))
				)
			
			self.__init__(intervalTuple)
		else:
			raise Exception('Wrong type of argument for Cube contructor {}'.format(intervals))
	
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
		chain = Chain(0)
		if self.dim == 0:
			return Chain(0)
		else:
			# We are sure that dim > 0
			if self.emb == 1:
				I = self.intervals[0]
				chain.isZero = False
				chain.emb = 1
				chain.dim = 0
				chain.cubes = [((I[0],I[0]),), ((I[1],I[1]),)]
				chain.coeff[((I[0],I[0]),)] = -1
				chain.coeff[((I[1],I[1]),)] = +1
				return chain
			else:
				I = Cube([self.intervals[0]])
				Q = Cube(self.intervals[1:])
				A = I.border() * Q
				B = I * Q.border()
				if I.dim % 2 == 0:
					return A + B
				else:
					return A - B
				
	def __lt__(self, other):
		if isinstance(other, Cube):
			if self.emb == other.emb:
				result = True
				for i in xrange(self.emb):
					sI = self.intervals[i]
					oI = other.intervals[i]
					result = result and (sI[0] in oI) and (sI[1] in oI) and (sI != oI)
					if not result:
						break
				return result
			else:
				raise Exception('Only cubes with the same embedding number can be compared.')
		else:
			raise Exception('A cube can only be compared with another cube.')				
	
	def __le__(self, other):
		if isinstance(other, Cube):
			if self.emb == other.emb:
				result = True
				for i in xrange(self.emb):
					sI = self.intervals[i]
					oI = other.intervals[i]
					result = result and (sI[0] in oI) and (sI[1] in oI)
					if not result:
						break
				return result
			else:
				raise Exception('Only cubes with the same embedding number can be compared.')
		else:
			raise Exception('A cube can only be compared with another cube.')				
	
	def __str__(self):
		return ' x '.join([str(Interval(*I)) for I in self.intervals])
	
	def __eq__(self, other):
		if isinstance(other, Cube):
			if self.emb != other.emb or self.dim != other.dim:
				return False
			else:
				return self.intervals == other.intervals
	
	def __mul__(self, other):
		if isinstance(other, Interval):
			return Cube(self.intervals + (other.interval,))
		elif isinstance(other, Cube):
			return Cube(self.intervals + other.intervals)
		elif isinstance(other, int):
			if other != 0:
				chain = Chain()
				chain.isZero = False
				chain.cubes.append(self.intervals)
				chain.dim = self.dim
				chain.emb = self.emb
				chain.coeff[self.intervals] = other
			return chain
		elif isinstance(other, Chain):
			chain = Chain()
			if not other.isZero:
				chain.isZero = False
				for c in other.cubes:
					chain.cubes.append(self.intervals + c)
					chain.coeff[self.intervals + c] = other.coeff[c]
				chain.emb = self.emb + other.emb
				chain.dim = self.dim + other.dim
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
			return self.intervals[index]
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
				if cube.intervals in bdr.cubes:
					return bdr.coeff[cube.intervals]
				else:
					return 0
		elif isinstance(cube, tuple) or isinstance(cube, list):
			return self.incidence(Cube(cube))
		else:
			raise Exception('Wrong type of argument {}'.format(cube))
			
class Chain(object):
	'''
	A Chain is a linear integer combination of cubes with the same 
	embeding number and dimension.
	However, for technical reasons some operations can be performed
	with chains of diferent dimensions (to enable the border calculation).
	
	The only way of making a chain is by operating cubes, intervals and
	integers in the right way.
	'''
	coeff = {}
	cubes = []
	isZero = True
	emb = 0
	dim = -1
	def __init__(self, cube=None):
		self.coeff = {}
		self.cubes = []
		self.isZero = True
		self.emb = 0
		self.dim = -1
		if cube:
			c = None
			if isinstance(cube, list) or isinstance(cube, tuple):
				c = Cube(cube)
			elif isinstance(cube, Cube):
				c = cube
			elif isinstance(cube, int) and cube == 0:
				pass
			else:
				raise Exception('Wrong type of argument in Chain constructor {}'.format(cube))
			
			if c:
				self.cubes.append(c.intervals)
				self.coeff[c.intervals] = +1
				self.isZero = False
				self.emb = c.emb
				self.dim = c.dim

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
				if other.intervals in self.cubes:
					self.coeff[other.intervals] += 1
				else:
					self.cubes.append(other.intervals)
					self.coeff[other.intervals] = +1
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
		elif isinstance(other, Cube):
			if self.isZero:
				return Chain(0)
			else:
				result = Chain(0)
				result.isZero = False
				result.emb, result.dim = self.emb + other.emb, self.dim + other.dim
				for c in self.cubes:
					cube = (Cube(c)*other).intervals
					result.cubes.append(cube)
					result.coeff[cube] = self.coeff[c]
				return result
		elif isinstance(other, Chain):
			if self.isZero or other.isZero:
				return Chain(0)
			else:
				products = []
				for c1 in self.cubes:
					for c2 in other.cubes:
						factor = Chain(0)
						factor.dim = self.dim + other.dim
						factor.emb = self.emb + other.emb
						factor.isZero = False
						cubeFactor = (Cube(c1)*Cube(c2)).intervals
						factor.cubes.append(cubeFactor)
						factor.coeff[cubeFactor] = self.coeff[c1] * other.coeff[c2]
						products.append(factor)
				result = Chain(0)
				for p in products:
					result += result + p
				return result
	
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
		
class CubicalComplex(object):
	emb = 0
	dim = -1
	cubes = {}
	vectorField = {}
	criticalCells = {}
	
	def __init__(self, emb):
		self.emb = emb
		self.dim = -1
		self.cubes = dict([(dim, []) for dim in xrange(emb + 1)])
		self.vectorField = {}
		self._acumDim = {}
		self.criticalCells = {}
		self.size = [0] * emb
		
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
				if newCube.dim in self.cubes:
					if newCube.intervals not in self.cubes[newCube.dim]:
						self.cubes[newCube.dim].append(newCube.intervals)
						
						self.dim = max(self.dim, newCube.dim)
						self.size = [max(newCube[i][1], self.size[i]) for i in xrange(self.emb)]
						for c in newCube.border().cubes:
							self.addCube(c)
				else:
					self.cubes[newCube.dim] = [newCube.intervals]
					self.dim = max(self.dim, newCube.dim)
					for c in newCube.border().cubes:
							self.addCube(c)
				
				for d in self.cubes:
					self.cubes[d].sort()
					self._acumDim[d] = sum(len(self.cubes[i]) for i in xrange(d+1))
				self._acumDim[-1] = 0
			else:
				raise Exception('Wrong embeding number for cube {}'.format(newCube))
		
	def __iter__(self):
		for d in self.cubes:
			for cube in self.cubes[d]:
				yield Cube(cube)
	
	def __getitem__(self, idx):
		if isinstance(idx, int):
			length = len(self)
			if idx < 0:
				return self[length + idx]
			elif idx >= length:
				raise IndexError('Wrong index {}'.format(idx))
			else:
				for d in xrange(self.dim + 1):
					if self._acumDim[d - 1] <= idx < self._acumDim[d]:
						return Cube(self.cubes[d][idx - self._acumDim[d - 1]])
		else:
			raise IndexError('Wrong type of index {}.'.format(idx))
	
	def __len__(self):
		return sum(len(self.cubes[d]) for d in self.cubes)
	
	def index(self, cube):
		if isinstance(cube, Cube):
			if cube.dim in self.cubes:
				if cube.intervals in self.cubes[cube.dim]:
					return self._acumDim[cube.dim - 1] + self.cubes[cube.dim].index(cube.intervals)
				else:
					raise ValueError('{} is not in the cubical complex.'.format(cube))
			else:
				raise ValueError('{} is not in the cubical complex.'.format(cube))
		else:
			if isinstance(cube, tuple) or isinstance(cube, list):
				return self.index(Cube(cube))

	def membraneLanguage(self):
		alphabet = {}
		alphabet['s'] = (len(self),)
		alphabet['d+'] = (len(self), len(self))
		alphabet['d-'] = (len(self), len(self))
		alphabet['U+'] = (len(self), len(self), self.emb)
		alphabet['U-'] = (len(self), len(self), self.emb)
		alphabet['C'] = (len(self),)
		alphabet['V+'] = (len(self), len(self))
		alphabet['V-'] = (len(self), len(self))
		alphabet['V'] = (len(self),)
		alphabet['F'] = (len(self), 2**self.emb+1)
		alphabet['D'] = (len(self), self.dim + 1)
		
		language = Language(alphabet)
		for name in language.objects:
			language.environment.append(name)
			
		return language
	
	def populateMembrane(self, membranes, index):
		for cube in self:
			cubeIdx = self.index(cube)
			membranes[index]['s'][cubeIdx] = 1
			for bdr in cube.border():
				bdrIdx = self.index(bdr)
				inc = cube.incidence(bdr)
				sgn = '+' if inc > 0 else '-'
				name = 'd{}'.format(sgn)
				membranes[index][name][bdrIdx, cubeIdx] = 1
				
				u = 2*(np.array(cube.center) - np.array(bdr.center))
				k = int((np.array(xrange(self.emb)) * u).sum())
				sgn = '+' if u[k] > 0 else '-'
				membranes[index]['U' + sgn][bdrIdx, cubeIdx, abs(k)] = 1
				
			membranes[index]['F'][cubeIdx, len(self.facets(cube))] = 1
			membranes[index]['D'][cubeIdx, cube.dim] = 1
			membranes[index]['C'][cubeIdx] = 1

	def facets(self, cube):
		if isinstance(cube, list) or isinstance(cube, tuple):
			dim = Cube(cube).dim
			if dim >= self.emb:
				return []
			else:
				return [Cube(c) for c in self.cubes[Cube(cube).dim + 1] if Cube(cube) <= Cube(c)]
		elif isinstance(cube, Cube):
			if cube.dim >= self.emb:
				return []
			else:
				return [Cube(c) for c in self.cubes[cube.dim + 1] if cube <= Cube(c)]
		else:
			raise Exception('Wrong type of argument {}.'.format(cube))

	def buildVectorField(self):
		self.vectorField = {}
		for cube in self:
			vectors = {}
			for facet in self.facets(cube):
				direction = map(
					lambda (x,y): int(2*x - 2*y),
					zip(facet.center, cube.center))
				dirIndex = sum(map(
					lambda (x,y): x*y,
					zip(xrange(self.emb), direction)))
				
				if dirIndex in vectors:
					vectors[dirIndex].append(facet)
				else:
					vectors[dirIndex] = [facet]

			for k in xrange(self.emb):
				if k in vectors:
					for facet in vectors[k]:
						if (cube.intervals not in self.vectorField.keys() and 
							cube.intervals not in  self.vectorField.values() and
							facet.intervals not in self.vectorField.keys() and
							facet.intervals not in self.vectorField.values()):
							self.vectorField[cube.intervals] = facet.intervals
							break
		
		self.buildCriticalCells()
		
	def buildCriticalCells(self):
		self.criticalCells = {}
		for dim in xrange(self.emb + 1):
			for cube in self.skeleton(dim):
				if self.isCritical(cube):
					if dim in self.criticalCells:
						if cube not in self.criticalCells:
							self.criticalCells[dim].append(cube.intervals)
					else:
						self.criticalCells[dim] = [cube.intervals]
						
	def skeleton(self, dim):
		if dim in self.cubes:
			for cube in self.cubes[dim]:
				yield Cube(cube)
		else:
			for cube in []:
				yield cube
			
	def isCritical(self, cube):
		if isinstance(cube, list) or isinstance(cube, tuple):
			return self.isCritical(Cube(cube))
		elif isinstance(cube, Cube):
			notCriticalCubes = self.vectorField.keys()
			notCriticalCubes.extend(self.vectorField.values())
			return cube.intervals not in notCriticalCubes
		else:
			raise Exception('Wrong type of argument {}'.format(cube))
	
	def V(self, chain):
		if isinstance(chain, Cube):
			if chain.intervals not in self.vectorField:
				return Chain(0)
			else:
				facet = Cube(self.vectorField[chain.intervals])
				return Chain(facet) * facet.incidence(chain) * (-1)
				
		if isinstance(chain, Chain):
			result = Chain(0)
			for c in chain.cubes:
				if c in self.vectorField:
					result += Chain(self.vectorField[c]) * chain.coeff[c]
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
		B = CC.fromBinaryArray(array)
		cells = CC.cells(B, graded=False)
		for cell in cells:
			self.addCube([(c/2,c/2) if c % 2 == 0 else ((c-1)/2, (c+1)/2) for c in cell])
	
	def __call__(self, value):
		if isinstance(value, int):
			return self.skeleton(value)
		else:
			raise Exception('Wrong type of argument {}'.format(value))
	
	def toCellComplex(self):
		cellCount = len(self)
		borderDict = {}
		dimDict = {}
		
		for i in xrange(cellCount):
			chain = K[i].border()
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
		
	def criticalComplex(self):
		cellCount = sum(len(self.criticalCells[d]) for d in self.criticalCells)
		criticalCells = []
		for d in self.criticalCells:
			for cube in self.criticalCells[d]:
				criticalCells.append(self.index(cube))
		borderDict = {}
		dimDict = {}
		for i in criticalCells:
			criticalBorder = K.Flow(K[i])[0].border()
			borderDict[criticalCells.index(i)] = dict(
				[(criticalCells.index(K.index(Cube(c))),criticalBorder.coeff[c]) 
					for c in criticalBorder.cubes])
			dimDict[criticalCells.index(i)] = K[i].dim
			
		return AH.CellComplex(cellCount, borderDict, dimDict)
	
	def homology(self, index=None):
		C = self.criticalComplex()
		C.buildVectorField()
		return C.criticalComplex().homology(index)	
		
if __name__ == '__main__':
	K = CubicalComplex(2)
	K.fromArray(np.array(
		[
			[1,1,1,1,1,0],
			[1,1,1,1,1,1],
			[1,1,0,1,1,1],
			[1,1,0,1,1,0],
			[1,1,1,1,1,0],
			[1,1,1,1,1,0],
			[1,1,1,0,0,0]
		]))
		
	K.buildVectorField()
	canvas = K.createSVGCanvas(squareSize=40)
	K.drawVectorFieldToSVGCanvas(canvas, opacity=0.2)
	K.drawToSVGCanvas(canvas, opacity=0.05)
	for d in K.criticalCells:
		for cube in K.criticalCells[d]:
			F = K.Flow(cube)
			F.drawToSVGCanvas(canvas, fill='green', opacity=1)
	canvas.display()
