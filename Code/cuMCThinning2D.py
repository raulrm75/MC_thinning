#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import mahotas
import CubicComplex as CC
import time
from PIL import Image
from math import ceil

import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule

import pdb, sys

class mcTSystem(object):
	shape = None
	cells = []
	cellArray = None
	cR = None
	cI = None
	cF = None
	mS = None
	mI = None
	mR = None
	sR = True
	dimension = 0
	module = None
	
	def __init__(self, shape=None):
		if shape != None:
			self.shape = shape
			self.cells = []
			self.cellArray = np.zeros(shape, np.int32)
			self.cR = np.ones(shape, np.int32)
			self.cI = np.zeros(shape, np.int32)
			self.cF = np.zeros(shape, np.int32)
			self.mS = {
				'src': np.zeros(shape, np.int32), 
				'dst': np.zeros(shape, np.int32)
				}
			self.mI = np.zeros(shape, np.int32)
			self.mR = {
				'src': np.zeros(shape, np.int32), 
				'dst': np.zeros(shape, np.int32)
				}
			self.sR = True
	
	def fromBinaryArray(self, binArray):
		#~ pdb.set_trace()
		self.cellArray = CC.fromBinaryArray(binArray).astype(np.int32)
		self.shape = self.cellArray.shape
		self.cR = np.ones(self.shape, np.int32)
		self.cI = np.zeros(self.shape, np.int32)
		self.cF = np.zeros(self.shape, np.int32)
		self.mS = {
			'src': np.zeros(self.shape, np.int32), 
			'dst': np.zeros(self.shape, np.int32)
			}
		self.mI = np.zeros(self.shape, np.int32)
		self.mR = {
			'src': np.zeros(self.shape, np.int32), 
			'dst': np.zeros(self.shape, np.int32)
			}
		self.sR = True
		
		print "Loading CUDA kernel sources..."
		with open('kernels.cu', 'r') as f:
			src = f.read()
			moduleSrc = src % {"rows": self.cellArray.shape[0], 
				"cols": self.cellArray.shape[1]}
			st = time.time()
			self.module = SourceModule(moduleSrc)
			print "CUDA kernel sources loaded in %.3f s" % (time.time() - st)
			
		self.calcFacetCounter()
		
	def fromImage(self, fileName, threshold):
		img = mahotas.imread(fileName)
		binArray = (img < threshold).astype(np.int32)
		self.fromBinaryArray(binArray)
			
	def calcFacetCounter(self):
		sumFacetsKernel = self.module.get_function("sumFacets")
		self.cF = np.zeros_like(self.cellArray)
		blockSize = (32, 32, 1)
		dimGrid = (1 + self.shape[0] / blockSize[0], 1 + self.shape[1] / blockSize[1])
		sumFacetsKernel(cuda.In(self.cellArray), cuda.InOut(self.cF), block=blockSize, grid=dimGrid)

			
	def markIsolatedCells(self):
		markIsolatedCellsKernel = self.module.get_function("markIsolatedCells")
		blockSize = (32, 32, 1)
		dimGrid = (1 + self.shape[0] / blockSize[0], 1 + self.shape[1] / blockSize[1])
		newIsolatedCells = np.zeros_like(self.cellArray)
		markIsolatedCellsKernel(cuda.In(self.cellArray), cuda.In(self.cF), cuda.InOut(self.mI), cuda.InOut(newIsolatedCells), block=blockSize, grid=dimGrid)
		return newIsolatedCells
	
	def getSimpleFacets(self):
		getSimpleFacetsKernel = self.module.get_function("getSimpleFacets")
		sFacets = np.zeros_like(self.cellArray)
		blockSize = (32, 32, 1)
		dimGrid = (1 + self.shape[0] / blockSize[0], 1 + self.shape[1] / blockSize[1])
		getSimpleFacetsKernel(cuda.In(self.cellArray), cuda.In(self.cF), cuda.InOut(sFacets), block=blockSize, grid=dimGrid)
		
		self.mS['dst'] = sFacets

	def getSimpleFaces(self):
		getSimpleFacesKernel = self.module.get_function("getSimpleFaces")
		sFaces = np.zeros_like(self.cellArray)
		blockSize = (32, 32, 1)
		dimGrid = (1 + self.shape[0] / blockSize[0], 1 + self.shape[1] / blockSize[1])
		getSimpleFacesKernel(cuda.In(self.cF), cuda.In(self.mS['dst']), 
			cuda.InOut(sFaces), block=blockSize, grid=dimGrid)
		self.mS['src'] = sFaces

	def markSimplePairs(self):
		self.getSimpleFacets()
		self.getSimpleFaces()
	
	def updateR(self):
		updateRKernel = self.module.get_function("updateR")
		blockSize = (32, 32, 1)
		dimGrid = (1 + self.shape[0] / blockSize[0], 1 + self.shape[1] / blockSize[1])
		updateRKernel(cuda.In(self.cellArray), cuda.InOut(self.cR), 
			block=blockSize, grid=dimGrid)
	
	def updateI(self, newIsolatedCells):
		updateIKernel = self.module.get_function("updateI")
		blockSize = (32, 32, 1)
		dimGrid = (1 + self.shape[0] / blockSize[0], 1 + self.shape[1] / blockSize[1])
		updateIKernel(cuda.In(newIsolatedCells), cuda.In(self.cR), 
			cuda.InOut(self.cI), block=blockSize, grid=dimGrid)
	
	def calcMedialPersistence(self):
		calcMedialPersistenceKernel = self.module.get_function("calcMedialPersistence")
		blockSize = (32, 32, 1)
		dimGrid = (1 + self.shape[0] / blockSize[0], 1 + self.shape[1] / blockSize[1])
		MPabs = np.zeros(self.shape, np.int32)
		MPrel = np.zeros(self.shape, np.float32)
		calcMedialPersistenceKernel(cuda.In(self.cR), cuda.In(self.cI), 
			cuda.InOut(MPabs), cuda.InOut(MPrel), 
			block=blockSize, grid=dimGrid)
		
		return MPabs, MPrel
		
	def markRemove(self, MPabs=None, MPrel=None, absT=None, relT=None):
		markRemoveKernel = self.module.get_function("markRemove")
		if MPabs != None and MPrel != None and absT != None and relT != None:
			self.mR['src'] = np.zeros(self.shape, dtype=np.int32)
			self.mR['dst'] = np.zeros(self.shape, dtype=np.int32)
			sR = np.zeros((1,), dtype=np.int32)
			blockSize = (32, 32, 1)
			dimGrid = (1 + self.shape[0] / blockSize[0], 1 + self.shape[1] / blockSize[1])
			markRemoveKernel(
				cuda.In(self.mS['dst']), cuda.In(self.mS['src']),
				cuda.In(self.cI), 
				cuda.In(MPabs), cuda.In(MPrel), np.int32(absT),
				np.float32(relT), np.int32(True),
				cuda.InOut(self.mR['dst']), cuda.InOut(self.mR['src']),
				cuda.InOut(sR), block=blockSize, grid=dimGrid)
		else:
			self.mR['src'] = np.zeros(self.shape, dtype=np.int32)
			self.mR['dst'] = np.zeros(self.shape, dtype=np.int32)
			MPabs = np.zeros(self.shape, np.int32)
			MPrel = np.ones(self.shape, np.float32)
			sR = np.zeros((1,), dtype=np.int32)
			blockSize = (32, 32, 1)
			dimGrid = (1 + self.shape[0] / blockSize[0], 1 + self.shape[1] / blockSize[1])
			markRemoveKernel(
				cuda.In(self.mS['dst']), cuda.In(self.mS['src']),
				cuda.In(self.cI),
				cuda.In(MPabs), cuda.In(MPrel), np.int32(0),
				np.float32(1), np.int32(False),
				cuda.InOut(self.mR['dst']), cuda.InOut(self.mR['src']),
				cuda.InOut(sR), block=blockSize, grid=dimGrid)

		self.sR = (sR[0] == 1)
	
	def doRemove(self):
		doRemoveKernel = self.module.get_function("doRemove")
		blockSize = (32, 32, 1)
		dimGrid = (1 + self.shape[0] / blockSize[0], 1 + self.shape[1] / blockSize[1])
		doRemoveKernel(cuda.In(self.mR['dst']), cuda.In(self.mR['src']), 
			cuda.InOut(self.cellArray), block=blockSize, grid=dimGrid)
		self.calcFacetCounter()
		
	def calcRI(self):
		it = 0
		newIsolatedArray = self.markIsolatedCells()
		while self.sR:
			it += 1
			print "Iteration: %d" % it
			self.markSimplePairs()
			self.markRemove()
			self.doRemove()
			newIsolatedArray = self.markIsolatedCells()
			self.updateI(newIsolatedArray)
			self.updateR()
	
	def doThinning(self, absT, relT, maxIter=None):
		#~ pdb.set_trace()
		it = 0
		st = time.time()
		newIsolatedArray = self.markIsolatedCells()
		while self.sR and (maxIter == None or it < maxIter):
			it += 1
			print "Iteration: %d" % it
			self.markSimplePairs()
			MPabs, MPrel = self.calcMedialPersistence()
			self.markRemove(MPabs, MPrel, absT, relT)
			self.doRemove()
			newIsolatedArray = self.markIsolatedCells()
			self.updateI(newIsolatedArray)
			self.updateR()
		print "Skeleton calculated in %.3f s" % (time.time() - st)
		print "Absolute threshold: %d, Relative threshold: %.3f" % (absT, relT)

if __name__ == '__main__':
	if len(sys.argv) < 2:
		print """
Usage: cuMCThinning2D fileName absT relT.
	fileName: Image file name.
	absT: Integer number used as absolute threshold.
	relT: Float number used as relative threshold.
		"""
	else:
		fileName = sys.argv[1]
		absT = None
		relT = None
		if len(sys.argv) == 4:
			absT = int(sys.argv[2])
			relT = float(sys.argv[3])
		elif len(sys.argv) == 3:
			absT = int(sys.argv[2])
			relT = 1.0
			
		img = Image.open(fileName).convert('L')
		imgArray = np.asarray(img)
		if absT is None:
			absT = 5 * max(imgArray.shape) / 100
			relT = 0.5
		T = (np.max(imgArray) + np.min(imgArray))/2
		imgArrayT = (imgArray < T).astype(np.uint8)
		P = mcTSystem()
		P.fromBinaryArray(imgArrayT)
		P.doThinning(absT, relT)
		result = P.cellArray[::2,::2].astype(np.uint8)*255 + 255 - 255*imgArrayT
		Image.fromarray(result).show()
	
	
	
			
	
