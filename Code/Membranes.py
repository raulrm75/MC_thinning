#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Membranes.py
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
import numpy as np
from itertools import product

def indexes(ndarray, onlyNonZero=True):
	for idx in product(*[xrange(s) for s in ndarray.shape]):
		if onlyNonZero:
			if ndarray[idx] != 0:
				yield idx
			else:
				continue
		else:
			yield idx

	
class Membrane(object):
	language = None
	objects = {}
	index = -1
	isEnvironment = False
	
	def __init__(self, language, index, isEnvironment=False):
		self.language = language
		self.index = index
		self.isEnvironment = isEnvironment
		# Initialize multiplicities
		for name in self.language.objects:
			if self.isEnvironment and (name in self.language.environment):
				shape = self.language.objects[name].shape
				self.objects[name] = np.ones(shape, np.int16)
			else:
				shape = self.language.objects[name].shape
				self.objects[name] = np.zeros(shape, np.int16)
				
	def  __getitem__(self, index):
		if isinstance(index, str):
			name = index
			return self.objects[name]
		elif isinstance(index, tuple):
			name, idx = index
			return self.objects[name][idx]
		else:
			raise IndexError('Index {} not found.'.format(index))
	
	def  __setitem__(self, index, value):
		if isinstance(index, str):
			name = index
			self.objects[name] = value
		elif isinstance(index, tuple):
			name, idx = index
			self.objects[name][idx] = value
		else:
			raise IndexError('Index {} not found.'.format(index))
	
	def __contains__(self, multiset):
		if isinstance(multiset, MultiSet):
			if bool(multiset):
				if multiset.isConstant():
					result = True
					for name, index in multiset.data:
						mult = multiset.data[name, index]
						result = (result and (
							(self.isEnvironment and name in self.language.environment)) 
							or self.objects[name][index] >= mult)
						if not result:
							break
					return result
				else:
					result = True
					for m in multiset:
						result = result and m in self
					return result
			else:
				return True
				
	def __iter__(self):
		for name in self.objects:
			for idx in indexes(self.objects[name]):
				yield MultiSet(self.language, {(name, idx): self.objects[name][idx]})
			
					
