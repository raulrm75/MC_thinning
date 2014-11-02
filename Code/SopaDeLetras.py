#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np

sopa = [
	'RTRFHJKLJKVBM',
	'SVENUSVCMMNHE',
	'DRTGSFCUSJJFR',
	'NFIBOTVVAADGC',
	'EBPHNFICTRXNU',
	'PHUMARTEURJBT',
	'TNJYRFVCRGFAI',
	'UKHGURVHNRDSO',
	'NOTULPVJOBADT',
	'OKBGSDVMMCHFR']

S = np.array([c for c in ''.join(sopa)]).reshape((len(sopa), len(sopa[0])))

planetas = ['MERCURIO', 'VENUS', 'TIERRA', 'MARTE', 'JUPITER', 
	'SATURNO', 'URANO', 'NEPTUNO', 'PLUTON']

for planeta in planetas:
	for row in xrange(10):
		for col in xrange(13 - len(planeta) + 1):
			found = ''.join(S[row, col:col+len(planeta)])
			
			if found == planeta:
				print '{} encontrado en la fila {}, columna {} en horizontal.'.format(planeta, row + 1, col + 1)
	
	for col in xrange(13):
		for row in xrange(10 - len(planeta) + 1):
			found = ''.join(S[row: row + len(planeta), col])
			
			if found == planeta:
				print '{} encontrado en la fila {}, columna {} en vertical.'.format(planeta, row + 1, col + 1)
