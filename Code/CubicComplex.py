#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from itertools import product
import simpleSVG

# Una celda es un np.array, con shape == (n,) 

def emb(cell):
	if isinstance(cell, tuple):
		return len(cell)
	else:
		return cell.shape[0]

def dim(cell):
	if isinstance(cell, tuple):
		c = np.array(cell, dtype=np.int16)
	else:
		c = cell
		
	return (c % 2).sum()

# Un complejo celular es un np.array con shape % 2 == (0, ..., 0) y dtype == np.bool

def cells(cellComplex, graded=False):
	indices = product(*[xrange(n) for n in cellComplex.shape])
	if not graded:
		CELLS = []
		for i in indices:
			if cellComplex[i]:
				CELLS.append(i)
		return tuple(CELLS)
	else:
		CELLS = {}
		for i in indices:
			c = np.array(i, dtype=np.int16)
			d = dim(c)
			if cellComplex[i]:
				if d in CELLS:
					CELLS[d].append(i)
				else:
					CELLS[d] = [i]
		
		return CELLS

def border(cell, cellComplex=None):
	oddIndices = []
	n = emb(cell)
	for i in xrange(n):
		if cell[i] % 2 == 1:
			oddIndices.append(i)
	
	borderCells = []
	for i in oddIndices:
		borderCells.append(tuple([cell[j] if j != i else cell[i] - 1 for j in xrange(n)]))
		borderCells.append(tuple([cell[j] if j != i else cell[i] + 1 for j in xrange(n)]))
	
	if cellComplex == None:
		return borderCells
	else:
		return [c for c in borderCells if cellComplex[c] == 1]

def facets(cell, shape, cellComplex=None):
	evenIndices = []
	n = emb(cell)
	evenIndices = [i for i in xrange(n) if cell[i] % 2 == 0]
	
	facetsCells = []
	for i in evenIndices:
		facetsCells.append(tuple([cell[j] if j != i else cell[i] - 1 for j in xrange(n)]))
		facetsCells.append(tuple([cell[j] if j != i else cell[i] + 1 for j in xrange(n)]))

	result = []
	c0 = np.zeros(n, np.int16)
	c1 = np.array(shape)
	
	for cell in facetsCells:
		c = np.array(cell)
		validCell = np.alltrue(c0 <= c) and np.alltrue(c < c1)
		if validCell:
			result.append(cell)
	
	if cellComplex == None:
		return result
	else:
		return [c for c in result if cellComplex[c] == 1]

def generators(cell):
	if isinstance(cell, tuple):
		_cell = np.array(cell)
	else:
		_cell = cell
	if dim(_cell) == 0:
		return [tuple(_cell)]
	else:
		n,d = emb(_cell), dim(_cell)
		base = [[0]*n for i in xrange(d)]
		oddIndices = [i for i in xrange(n) if _cell[i] % 2 == 1]
		
		for i in xrange(d):
			base[i][oddIndices[i]] = 1
		
		P = product(*[product([-1,1], [b]) for b in base])
		G = [sum(p[i][0] * np.array(p[i][1]) for i in xrange(d)) for p in P]
		
		return [tuple(_cell + g) for g in G] 
		
def fromBinaryArray(I):
	K = np.zeros(tuple(2*np.array(I.shape)-1), dtype=np.bool)
	n = len(K.shape)
	K[tuple([slice(0,K.shape[i],2) for i in xrange(n)])] = I
	cellModels = calcCellModels(n)
	dimensionMask = calcDimensionMask(K)
	M = K.astype(np.uint8)
	
	for p in xrange(1, n + 1):
		M += sum(shift(M,+1,c) + shift(M,-1,c) for c in cellModels[p])*dimensionMask[p]/2
					
	return M.astype(np.bool)

def fromImage(imageFileName, threshold):
	img = mahotas.imread(imageFileName)
	return fromBinaryArray(img < threshold)

def drawCubicComplex(cubicComplex, fileName='temp.svg', squareSize=10):
	R = 2
	n = len(cubicComplex.shape)
	assert(n == 2)
	
	CELLS = cells(cubicComplex, graded=True)
	canvas = simpleSVG.svg_class(fname=fileName, bbx=squareSize + squareSize*cubicComplex.shape[1]/2, bby=squareSize + squareSize*cubicComplex.shape[0]/2)
	canvas.scale()
	
	for d in CELLS:
		if d == 0:
			for cell in CELLS[0]:
				y = R + 2 + cell[0] / 2 * squareSize 
				x = R + 2 + cell[1] / 2 * squareSize
				canvas.circle(x, y, R,fill='red')
		if d == 1:
			for cell in CELLS[1]:
				generator = generators(cell)
				y = [R + 2 + g[0] / 2 * squareSize for g in generator]
				x = [R + 2 + g[1] / 2 * squareSize for g in generator]
				if x[0] == x[1]:
					y.sort()
					y[0] += R
					y[1] -= R
				if y[0] == y[1]:
					x.sort()
					x[0] += R
					x[1] -= R
					
				canvas.line(x[0], y[0], x[1], y[1], fill='green')

		if d == 2:
			for cell in CELLS[2]:
				generator = generators(cell)
				y = [R + 2 + g[0] / 2 * squareSize for g in generator]
				x = [R + 2 + g[1] / 2 * squareSize for g in generator]
				bx = sum(x)/4
				by = sum(y)/4
				px = bx - squareSize/2 + R
				py = by - squareSize/2 + R
				
				canvas.rect(px, py, squareSize - 2*R, squareSize - 2*R, fill='blue')
	
	canvas.close()
	
def drawPixels(cubicComplex, fileName):
	I = np.zeros(np.array(cubicComplex.shape)/2, dtype=np.int8)
	I = 255 - 255 * cubicComplex[::2,::2]
	mahotas.imsave(fileName, I.astype(np.uint8))

def shift(cubicComplex, sign, direction):
	n = len(cubicComplex.shape)
	result = np.zeros(cubicComplex.shape, dtype=cubicComplex.dtype)
	idx = abs(np.sum(np.array(direction) * np.arange(n)))
	if sign > 0:
		S1 = tuple([
			slice(0,cubicComplex.shape[i]) if i != idx else 
			slice(1, cubicComplex.shape[i]) for i in xrange(n)])
		S2 = tuple([
			slice(0,cubicComplex.shape[i]) if i != idx else 
			slice(0, cubicComplex.shape[i] - 1) for i in xrange(n)])
	else:
		S1 = tuple([
			slice(0,cubicComplex.shape[i]) if i != idx else 
			slice(0, cubicComplex.shape[i] - 1) for i in xrange(n)])
		S2 = tuple([
			slice(0, cubicComplex.shape[i]) if i != idx else 
			slice(1, cubicComplex.shape[i]) for i in xrange(n)])
	result[S1] = cubicComplex[S2]	
	
	return result

def calcCellModels(dimension):
	cellList = list(product([0,1], repeat=dimension))
	cellModels = {}
	for cell in cellList:
		d = dim(cell)
		if d in cellModels:
			cellModels[d].append(cell)
		else:
			cellModels[d] = [cell]
	
	return cellModels

def calcDimensionMask(cubicComplex):
	n = len(cubicComplex.shape)
	cellModels = calcCellModels(n)
	dimensionMask = {}
	for p in xrange(n+1):
		dimensionMask[p] = np.zeros(cubicComplex.shape, dtype=np.bool)
		for cell in cellModels[p]:
			S = tuple([slice(0, cubicComplex.shape[i], 2) if cell[i] == 0 else slice(1, cubicComplex.shape[i] - 1, 2) for i in xrange(n)])
			dimensionMask[p][S] = 1
			
	
	return dimensionMask

if __name__ == "__main__":
	I = np.array([
		[0,1,1,1,1,1,1],
		[0,1,0,0,0,1,1],
		[0,1,0,0,1,1,1],
		[0,1,0,0,1,1,1],
		[0,1,0,0,1,1,0],
		[0,1,0,0,1,1,0],
		[0,1,1,1,1,0,0]])
	
	K = np.zeros(tuple(2*np.array(I.shape)-1), dtype=np.bool)
	n = len(K.shape)
	K[tuple([slice(0,K.shape[i],2) for i in xrange(n)])] = I
	cellModels = calcCellModels(n)
	dimensionMask = calcDimensionMask(K)
	M = K.astype(np.uint8)
	for p in xrange(1, n + 1):
		M += sum(shift(M,+1,c) + shift(M,-1,c) for c in cellModels[p])*dimensionMask[p]/2
	print M		
	drawCubicComplex(M)
