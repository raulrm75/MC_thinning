#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Rules.py
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
from Language import Language, Index, MultiSet
from itertools import product, combinations
import pdb

class Rule(object):
	language = None
	promoters = None
	inhibitors = None
	membraneIndex1 = None
	membraneIndex2 = None
	commData1 = None
	commData2 = None
	indexes = []
	availableIndexes = []
	membranes = {}
	index = None
	instance = {}
	constraints = []
	
	def __init__(self, language,
		promoters=None, inhibitors=None, 
		membraneIndex1=-1, commData1=None, commData2=None, membraneIndex2=-1,
		index=None, instance={}, constraints=[]):
			
		self.language = language
		self.index = index
		self.instance = instance
		self.constraints = constraints
		
		if promoters == None:
			self.promoters = MultiSet(language)
		else:
			self.promoters = promoters
		
		if inhibitors == None:
			self.inhibitors = MultiSet(language)
		else:
			self.inhibitors = inhibitors
		
		if membraneIndex1 == -1:
			raise SyntaxError('Argument missing')
		else:
			self.membraneIndex1 = membraneIndex1
		
		if membraneIndex2 == -1:
			raise SyntaxError('Argument missing')
		else:
			self.membraneIndex2 = membraneIndex2
		
		if commData1 == None:
			self.commData1 = MultiSet(language)
		else:
			self.commData1 = commData1
		
		if commData2 == None:
			self.commData2 = MultiSet(language)
		else:
			self.commData2 = commData2
		
		self.indexes = []
		for attrib in self.__dict__:
			if isinstance(self.__dict__[attrib], MultiSet):
				mset = self.__dict__[attrib]
				for index in mset.indexes:
					if index not in self.indexes:
						self.indexes.append(index)
		
		self.availableIndexes = list(product(*[xrange(*idx.interval) for idx in self.indexes]))
		
		
	def addConstraint(self, constraint):
		self.constraints.append(constraint)
	
	def evalConstraints(self, values):
			result = True
			valueDict = dict([(idx.name, values[idx]) for idx in values])
			for constraint in self.constraints:
				result = result and eval(constraint, valueDict)
				if not result:
					break
			return result
				
	def __mul__(self, other):
		if isinstance(other, Rule):
			s1 = set([self.membraneIndex1, self.membraneIndex2])
			s2 = set([other.membraneIndex1, other.membraneIndex2])
			if s1.intersection(s2) == set([]):
				return True
			else:
				for i in [1, 2]:
					for j in [1, 2]:
						selfMI = self.__dict__['membraneIndex{}'.format(i)]
						otherMI = other.__dict__['membraneIndex{}'.format(j)]
						if  selfMI == otherMI:
							selfCD = self.__dict__['commData{}'.format(i)]
							otherCD = other.__dict__['commData{}'.format(j)]
							collide = False
							for name, index in selfCD:
								if (name, index) in otherCD:
									collide = True
									break
						if collide:
								break
					if collide:
						break
				return collide
						
	def __iter__(self):
		for index in self.availableIndexes:
			D = dict([(self.indexes[i], index[i]) for i in xrange(len(index))])
			if self.evalConstraints(D):
				yield self.eval(D)
			else:
				continue
	
	def eval(self, values):
		result = Rule(
			self.language,
			promoters=self.promoters.eval(values),
			inhibitors=self.inhibitors.eval(values),
			membraneIndex1=self.membraneIndex1,
			commData1=self.commData1.eval(values),
			commData2=self.commData2.eval(values),
			membraneIndex2=self.membraneIndex2,
			index=self.index, 
			instance=values,
			constraints=self.constraints)
		
		#~ indexEv = []
		#~ for idx in self.index:
			#~ if isinstance(idx, int):
				#~ indexEv.append(idx)
			#~ elif isinstance(idx, Index):
				#~ indexEv.append(values[idx])
		#~ result.index = tuple(indexEv)
		
		return result
	
	def isConstant(self):
		return self.indexes == []
	
	def __getitem__(self, index):
		return self.eval(index)
	
	def shortStr(self):
		return 'R[{}][{}]'.format(
			self.index,
			map(lambda x: '{} = {}'.format(x[0].name, x[1]), self.instance.items()))
		
	def __str__(self):
		if self.promoters != 0:
			proStr = '{{{}}}'.format(self.promoters)
		else:
			proStr = ''
		if self.inhibitors != 0:
			inhStr = '¬{{{}}}'.format(self.inhibitors)
		else:
			inhStr = ''
		
		if proStr != '' and inhStr != '':
			result = 'R[{}] = ({} {} | {}, {} / {}, {})'.format(
				self.index, 
				proStr,
				inhStr,
				self.membraneIndex1,
				self.commData1,
				self.commData2,
				self.membraneIndex2)
		elif proStr != '' and inhStr == '':
			result = 'R[{}] = ({}| {}, {} / {}, {})'.format(
				self.index,
				proStr,
				self.membraneIndex1,
				self.commData1,
				self.commData2,
				self.membraneIndex2)
		elif proStr == '' and inhStr != '':
			result = 'R[{}] = ({}| {}, {} / {}, {})'.format(
				self.index,
				inhStr,
				self.membraneIndex1,
				self.commData1,
				self.commData2,
				self.membraneIndex2)
		else:
			result = 'R[{}] = ({}, {} / {}, {})'.format(
				self.index,
				self.membraneIndex1,
				self.commData1,
				self.commData2,
				self.membraneIndex2)
		
		if self.indexes == []:
			return result
		else:
			idxStr = ', '.join(['{} <= {} < {}'.format(
				idx.interval[0], 
				idx.name, 
				idx.interval[1]) for idx in self.indexes])
				
			result = result + ' for ' + idxStr
			if self.constraints:
				result = result + ' and ' + ', '.join(self.constraints)
			return result
	
	
		
if __name__ == '__main__':
	import CubicalHomology as CH
	K = CH.CubicalComplex(2)
	K.addCube(CH.Cube([(0,1),(0,1)]))
	K.addCube(CH.Cube([(0,1),(1,2)]))
	K.addCube(CH.Cube([(1,2),(0,1)]))
	
	L = K.membraneLanguage()
	i1 = Index('i1', (0, len(K)))
	j1 = Index('j1', (0, len(K)))
	i2 = Index('i2', (0, len(K)))
	j2 = Index('j2', (0, len(K)))
	R = Rule(L,
	promoters=L['C'][j1],
	inhibitors=None,
	membraneIndex1=1,
	commData1=L['s'][i1] + L['s'][j1] + L['s'][i2] + L['s'][j2] + L['d-'][i1,j1] + L['d-'][i1,j2] + L['d+'][i2,j2] + L['V+'][i1,j2],
	commData2=L['s'][j1] + L['s'][i2] + L['d-'][i2,j1],
	membraneIndex2=0,
	index=6,
	constraints=map(
		lambda x: ' != '.join(x),
		combinations(['{}{}'.format(n, i) for n in ['i', 'j'] for i in [1,2]], 2)
		))
	
	for idx in  R.indexes:
		print idx
