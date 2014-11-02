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
import numpy as np
from itertools import product, combinations
from random import choice
import CubicalHomology as CH
from Language import Language, Index, MultiSet
from Rules import Index, Rule
from TopSort import topSort
import pdb

class TissuePSystem(object):
	language = None
	membranes = {}
	rules = {}
	priorities = []
	
	def __init__(self, language, membraneCount):
		self.language = language
		self.membranes = {}
		for i in xrange(membraneCount + 1):
			self.membranes[i] = {}
			for name in self.language:
				shape = self.language.objects[name].shape
				if i == 0 and name in self.language.environment:
					self.membranes[i][name] = np.ones(shape, np.int16)
				else:
					self.membranes[i][name] = np.zeros(shape, np.int16)
		self.rules = {}
		self.priorities = []
	
	def addRule(self, rule):
		self.rules[rule.index] = rule
			
	
	def ruleMultiplicity(self, rule):
		if rule.isConstant():
			test = {
				'promoters': True,
				'inhibitors': True,
				'commData1': True,
				'commData2': True
			}
			pro = rule.promoters
			inh = rule.inhibitors
			m1 = rule.membraneIndex1
			m2 = rule.membraneIndex2
			c1 = rule.commData1
			c2 = rule.commData2
			
			# Test promoters
			for name, index in pro:
				mult = pro[name, index]
				if (name in self.language.environment) and m1 == 0:
					continue
				else:
					test['promoters'] = test['promoters'] and self.membranes[m1][name][index] >= mult
				if not test['promoters']:
					return 0
			
			# Test inhibitors
			for name, index in inh:
				mult = inh[name, index]
				if (name in self.language.environment) and m1 == 0:
					return 0
				else:
					test['inhibitors'] = test['inhibitors'] and self.membranes[m1][name][index] < mult
				if not test['inhibitors']:
					return 0
			
			# Test comm. data 1
			test['commData1'] = None
			for name, index in c1:
				mult = c1[name, index]
				if (name in self.language.environment) and m1 == 0:
					continue
				else:
					if test['commData1'] == None:
						test['commData1'] = self.membranes[m1][name][index] / mult
					else:
						test['commData1'] = min(test['commData1'], self.membranes[m1][name][index] / mult)
				if test['commData1'] != None and test['commData1'] == 0:
					return 0
			
			# Test comm. data 2
			test['commData2'] = None
			for name, index in c2:
				mult = c2[name, index]
				if (name in self.language.environment) and m2 == 0:
					continue
				else:
					if test['commData2'] == None:
						test['commData2'] = self.membranes[m2][name][index] / mult
					else:
						test['commData2'] = min(test['commData2'], self.membranes[m2][name][index] / mult)
				if test['commData2'] != None and test['commData2'] == 0:
					return 0
			
			if test['commData1'] != None and test['commData2'] != None:
				return min(test['commData1'], test['commData2'])
			elif test['commData1'] == None:
				return test['commData2']
			elif test['commData2'] == None:
				return test['commData1']
			else:
				return 0
		else:
			raise Exception('Only can test selection of constant rules.')
			
	def ruleSelection(self):
		
		#~ def allPos(membranes):
			#~ return np.array([
				#~ (membranes[i][name] >= 0).all() for i in membranes 
				#~ for name in self.language]).all()
		
		availableRules = {}
		print '*'*80
		print 'Calculating available rules...'
		for idx in self.rules:
			print '\t... for rule index {}'.format(idx)
			for rule in self.rules[idx]:
				mult = self.ruleMultiplicity(rule)
				if mult > 0:
					if idx in availableRules:
						availableRules[idx][rule] = mult
					else:
						availableRules[idx] = {rule:mult}
		print '{} available rules found.'.format(len(availableRules.items()))
		
		selectedRules = {}
		priorities = topSort(self.priorities)
		print '*'*80
		print 'Selecting rules not afected by priority...'
		for idx in set(availableRules.keys()) - set(priorities):
			print '\t... for rule index {}'.format(idx)
			for rule in availableRules[idx]:
				mult = availableRules[idx][rule]
				for appMult in xrange(1, mult + 1):
						applied = self.apply(rule)
						if not applied:
							appMult -= 1
							break
				
				if appMult > 0:
					print '{} * {}'.format(rule.shortStr(), appMult)
					raw_input('>>')	
					if idx in selectedRules:
						selectedRules[idx][rule] = appMult
					else:
						selectedRules[idx] = {rule:appMult}
		print '*'*80
		print 'Selecting rules afected by priority...'
		for idx in priorities:
			print '\t... for rule index {}'.format(idx)
			if idx in availableRules:
				for rule in availableRules[idx]:
					mult = availableRules[idx][rule]
					for appMult in xrange(1, mult + 1):
						applied = self.apply(rule)
						if not applied:
							appMult -= 1
							break
					
					if appMult > 0:
						print '{} * {}'.format(rule.shortStr(), appMult)
						raw_input('>>')	
						if idx in selectedRules:
							selectedRules[idx][rule] = appMult
						else:
							selectedRules[idx] = {rule:appMult}
		
		return selectedRules
	
	def apply(self, rule, ruleMult=1):
		m1 = rule.membraneIndex1
		m2 = rule.membraneIndex2
		c1 = rule.commData1
		c2 = rule.commData2

		if self.ruleMultiplicity(rule) > 0:
			for name, index in c1:
				mult = c1[name, index]
				if not (m1 == 0 and name in self.language.environment):
					#~ print 'Removing {}[{}] * {} from {}'.format(
						#~ name, index, mult, m1)
					self.membranes[m1][name][index] -= mult
				if not (m2 == 0 and name in self.language.environment):
					#~ print 'Adding {}[{}] * {} to {}'.format(
						#~ name, index, mult, m2)
					self.membranes[m2][name][index] += mult
			
			for name, index in c2:
				mult = c2[name, index]
				if not (m2 == 0 and name in self.language.environment):
					#~ print 'Removing {}[{}] * {} from {}'.format(
						#~ name, index, mult, m2)
					self.membranes[m2][name][index] -= mult
				if not (m1 == 0 and name in self.language.environment):
					#~ print 'Adding {}[{}] * {} to {}'.format(
						#~ name, index, mult, m1)
					self.membranes[m1][name][index] += mult
			return True
		else:
			return False
		
	def unapply(self, rule, ruleMult=1):
		m1 = rule.membraneIndex1
		m2 = rule.membraneIndex2
		c1 = rule.commData1
		c2 = rule.commData2
		
		for name, index in c1:
			mult = c1[name, index]
			if not (m1 == 0 and name in self.language.environment):
				#~ print 'Adding {}[{}] * {} to {}'.format(
					#~ name, index, mult, m1)
				self.membranes[m1][name][index] += mult
			if not (m2 == 0 and name in self.language.environment):
				#~ print 'Removing {}[{}] * {} from {}'.format(
					#~ name, index, mult, m2)
				self.membranes[m2][name][index] -= mult
		
		for name, index in c2:
			mult = c2[name, index]
			if not (m2 == 0 and name in self.language.environment):
				#~ print 'Adding {}[{}] * {} to {}'.format(
					#~ name, index, mult, m2)
				self.membranes[m2][name][index] += mult
			if not (m1 == 0 and name in self.language.environment):
				#~ print 'Removing {}[{}] * {} from {}'.format(
					#~ name, index, mult, m1)
				self.membranes[m1][name][index] -= mult
		
	def run(self):
		selections = []
		Halt = False
		step = 1
		while not Halt:
			sel = P.ruleSelection()
			if sel:
				selections.append(sel)
			else:
				Halt = True
		
		if Halt:
			return selections
		
					
if __name__ == '__main__':
	K = CH.CubicalComplex(2)
	K.addCube(CH.Cube([(0,0), (0,1)]))
	K.addCube(CH.Cube([(0,1), (1,2)]))
	K.addCube(CH.Cube([(0,1), (0,0)]))
	K.addCube(CH.Cube([(1,1), (0,1)]))
	K.addCube(CH.Cube([(1,2), (0,0)]))
	K.addCube(CH.Cube([(2,2), (0,1)]))
	K.addCube(CH.Cube([(1,2), (1,1)]))
	
	L = K.membraneLanguage()
	P = TissuePSystem(L, 1)
	K.populateMembrane(P.membranes, 1)
	
	i = Index('i', (0, len(K)))
	j = Index('j', (0, len(K)))
	i1 = Index('i1', (0, len(K)))
	j1 = Index('j1', (0, len(K)))
	i2 = Index('i2', (0, len(K)))
	j2 = Index('j2', (0, len(K)))
	P.addRule(Rule(L,
		promoters=None,
		inhibitors=L['V'][i] + L['V'][j],
		membraneIndex1=1,
		commData1=L['s'][i] + L['s'][j] + L['d-'][i,j] + L['U+'][i,j,0] + L['D'][i,0],
		commData2=L['s'][i] + L['s'][j] + L['d-'][i,j] + L['U+'][i,j,0] + L['D'][i,0] + L['V+'][i,j] + L['V'][i] + L['V'][j],
		membraneIndex2=0,
		index=1,
		constraints=['i < j']))
	
	P.addRule(Rule(L,
		promoters=None,
		inhibitors=L['V'][i] + L['V'][j],
		membraneIndex1=1,
		commData1=L['s'][i] + L['s'][j] + L['d-'][i,j] + L['U+'][i,j,1] + L['D'][i,0],
		commData2=L['s'][i] + L['s'][j] + L['d-'][i,j] + L['U+'][i,j,1] + L['D'][i,0] + L['V+'][i,j] + L['V'][i] + L['V'][j],
		membraneIndex2=0,
		index=2,
		constraints=['i < j']))
	
	P.addRule(Rule(L,
		promoters=None,
		inhibitors=L['V'][i] + L['V'][j],
		membraneIndex1=1,
		commData1=L['s'][i] + L['s'][j] + L['d-'][i,j] + L['U+'][i,j,0] + L['D'][i,1],
		commData2=L['s'][i] + L['s'][j] + L['d-'][i,j] + L['U+'][i,j,0] + L['D'][i,1] + L['V+'][i,j] + L['V'][i] + L['V'][j],
		membraneIndex2=0,
		index=3,
		constraints=['i < j']))
	
	P.addRule(Rule(L,
		promoters=None,
		inhibitors=L['V'][i] + L['V'][j],
		membraneIndex1=1,
		commData1=L['s'][i] + L['s'][j] + L['d+'][i,j] + L['U+'][i,j,1] + L['D'][i,1],
		commData2=L['s'][i] + L['s'][j] + L['d+'][i,j] + L['U+'][i,j,1] + L['D'][i,1] + L['V-'][i,j] + L['V'][i] + L['V'][j],
		membraneIndex2=0,
		index=4,
		constraints=['i < j']))
	
	P.addRule(Rule(L,
	promoters=L['V'][i],
	inhibitors=None,
	membraneIndex1=1,
	commData1=L['C'][i],
	commData2=None,
	membraneIndex2=0,
	index=5))

	P.addRule(Rule(L,
	promoters=None,
	inhibitors=L['C'][i1] + L['C'][j2],
	membraneIndex1=1,
	commData1=L['s'][i1] + L['s'][j1] + L['s'][i2] + L['s'][j2] + L['d-'][i1,j1] + L['d-'][i1,j2] + L['d+'][i2,j2] + L['V+'][i1,j2],
	commData2=L['s'][j1] + L['s'][i2] + L['d-'][i2,j1],
	membraneIndex2=0,
	index=6,
	constraints=['i1 < j1', 'i2 < j2', 'i1 != i2', 'j1 != j2']))
	
	P.addRule(Rule(L,
	promoters=None,
	inhibitors=L['C'][i1] + L['C'][j2],
	membraneIndex1=1,
	commData1=L['s'][i1] + L['s'][j1] + L['s'][i2] + L['s'][j2] + L['d+'][i1,j1] + L['d-'][i1,j2] + L['d+'][i2,j2] + L['V+'][i1,j2],
	commData2=L['s'][j1] + L['s'][i2] + L['d+'][i2,j1],
	membraneIndex2=0,
	index=7,
	constraints=['i1 < j1', 'i2 < j2', 'i1 != i2', 'j1 != j2']))
	
	P.priorities = [(1,2), (1,4), (2,3), (3,4)]
	
	selections = P.run()
	for step in xrange(len(selections)):
		print 'STEP {}'.format(step+1)
		for i in selections[step]:
			for rule in selections[step][i]:
				mult = selections[step][i][rule]
				print '{} ...\n {} times'.format(rule.shortStr(), mult)
				raw_input('>>')
