#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Language.py
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
from itertools import product

class LanguageObject(object):
	name = ''
	shape = ()
	
	def __init__(self, name, shape):
		self.name = name
		self.shape = shape
	
	def __str__(self):
		if len(self.shape) == 0:
			indexStr = ''
		elif len(self.shape) == 1:
			indexStr = 'i'
		else:
			indexStr = ', '.join(['i{}'.format(k+1) for k in xrange(len(self.shape))])
		
		if indexStr == '':
			return name
		else:
			if len(self.shape) == 1:
				return '{0}[{1}] for 0 <= {1} < {2}'.format(
					self.name, indexStr, self.shape[0])
			else:
				return '{0}[{1}] for {2} <= [{1}] < {3}'.format(
					self.name, indexStr, [0]*len(self.shape), list(self.shape))
	
	def __getitem__(self, index):
		if (isinstance(index, int) or isinstance(index, Index)) and len(self.shape) == 1:
			return MultiSet(self, {(self.name, index):1})
		elif isinstance(index, tuple):
			if len(index) != len(self.shape):
				raise SystaxError('Index len mismatch.')
			else:
				return MultiSet(self, {(self.name, index):1})
		else:
			raise Exception('Wrong type of argument {}.'.format(index))
			
class Language(object):
	objects = {}
	environment = []
	
	def __init__(self, objects, environment=[]):
		if isinstance(objects, dict):
			for name in objects:
				shape = objects[name]
				self.objects[name] = LanguageObject(name, shape)
				if name in environment:
					self.environment.append(name)
		else:
			raise Exception('Wrong type of argument {}'.format(objects))
	
	def __getitem__(self, name):
		if name in self.objects:
			return self.objects[name]
		else:
			raise KeyError('{} is not a Language Object'.format(name))
	
	def __iter__(self):
		for name in self.objects:
			yield name

	
class Index(object):
	name = ''
	interval = ()
	
	def __init__(self, name, interval):
		self.name = name
		self.interval = interval
	
	def __str__(self):
		return str(self.name)
	
	def __iter__(self):
		for i in range(*self.interval):
			yield i

class MultiSet(object):
	language = None
	data = {}
	indexes = []
	
	def __init__(self, language, data={}):
		self.language = language
		self.data = data
		self.indexes = []
		for name, index in self.data:
			if isinstance(index, Index) and index not in self.indexes:
				self.indexes.append(index)
			elif isinstance(index, tuple):
				for idx in index:
					if isinstance(idx, Index) and idx not in self.indexes:
						self.indexes.append(idx)
	
	def __mul__(self, other):
		if isinstance(other, int):
			result = {}
			if other > 0:
				for name, index in self.data:
					result[name, index] = other * self.data[name, index]
			return MultiSet(self.language, result) 
		else:
			raise SyntaxError('Wrong type of operand.')
	
	def eval(self, values):
		result = {}
		for name, index in self.data:
			if isinstance(index, int):
				if (name, index) in result:
					result[name, index] += self.data[name, index]
				else:
					result[name, index] = self.data[name, index]
			elif isinstance(index, Index):
				if index in values:
					if (name, values[index]) in result:
						result[name, values[index]] += self.data[name, index]
					else:
						result[name, values[index]] = self.data[name, index]
				else:
					result[name, index] = self.data[name, index]
			elif isinstance(index, tuple):
				indexEv = []
				for idx in index:
					if isinstance(idx, int):
						indexEv.append(idx)
					elif isinstance(idx, Index):
						if idx in values:
							indexEv.append(values[idx])
						else:
							indexEv.append(idx)
				result[name, tuple(indexEv)] = self.data[name, index]
				
		return MultiSet(self.language, result)
		
	def __add__(self, other):
		if isinstance(other, MultiSet):
			result = {}
			
			for name, index in self.data:
				result[name, index] = self.data[name, index]
				
			
			for name, index in other.data:
				if (name, index) in result:
					result[name, index] += other.data[name, index]
				else:
					result[name, index] = other.data[name, index]
			
			return MultiSet(self.language, result)
		else:
			raise SyntaxError('Wrong type of operand.')
	
	def __iter__(self):
		indexes = product(*[xrange(*idx.interval) for idx in self.indexes])
		for index in indexes:
			yield self.eval(dict([(self.indexes[i], index[i]) for i in xrange(len(index))]))
	
	def isConstant(self):
		return len(self.indexes) == 0
	
	def __nonzero__(self):
		return bool(self.data)
		
	def __str__(self):
		if self.data == []:
			return '0'
		else:
			L = []
			orderedData = self.data.keys()
			orderedData.sort()
			for name, idx in orderedData:
				mult = self.data[name, idx]
				if isinstance(idx, int):
					idxStr = str(idx)
				elif isinstance(idx, Index):
					idxStr = idx.name
				elif isinstance(idx, tuple):
					idxStr = ', '.join([str(i) if isinstance(i, int) else i.name for i in idx])
				L.append('{}[{}]{}'.format(
					name, idxStr,
					'' if mult == 1 else ' * {}'.format(mult)))
			if len(L) == 1:
				return L[0]
			else:
				return ' '.join(L)
	def __eq__(self, other):
		if isinstance(other, int) and other == 0:
			return self.data == {}
		elif isinstance(other, MultiSet):
			return self.data == other.data
		else:
			return False
	
	def __ne__(self, other):
		return not (self == other)
	
	def __iter__(self):
		for (name, index) in self.data:
			yield (name, index)
	
	def __getitem__(self, index):
		return self.data[index]
	
	def fromMembrane(self, membrane):
		for name in membrane:
			for index in product(*[xrange(s) for s in membrane[name].shape]):
				if membrane[name][index] > 0:
					self.data[name, index] = membrane[name][index]
					
if __name__ == '__main__':
	L = Language({'a':(10,),'b':(10,10)})
	i = Index('i', (0,10))
	j = Index('j', (0,10))
	print map(str, L['b'][j,i].eval({'j':1}).indexes)
