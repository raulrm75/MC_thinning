#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import ArrayCubicalHomology as ACH
import time
import os
from subprocess import call
import binvox_rw as bv
import visvis as vv

import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule

import pdb, sys

def imToGS(im):
		img = np.zeros(im.shape[:-1], dtype=im.dtype)
		img = 0.21 * im[...,0] + 0.71 * im[...,1] + 0.07 * im[...,2]
		return img
		
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
	voxelModel = None
	
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
		K = ACH.CubicalComplex()
		K.fromArray(binArray.astype(np.int8))
		self.cellArray = K.cubes.astype(np.int32)
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
		N1 = self.shape[0]
		N2 = self.shape[1]
		N3 = self.shape[2] if len(self.shape) == 3 else 1
		
		with open('kernels-nD.cu', 'r') as f:
			src = f.read()
			moduleSrc = src.format(N1=N1, N2=N2, N3=N3)
			st = time.time()
			self.module = SourceModule(moduleSrc)
			print "CUDA kernel sources loaded in %.3f s" % (time.time() - st)
			
		self.calcFacetCounter()
		
	def fromImage(self, fileName, threshold):
		im = vv.imread(fileName)
		img = imToGS(im)
		binArray = (img < threshold).astype(np.int32)
		self.fromArray(binArray)
	
	def fromMesh(self, fileName, resolution=128):
		# Voxelize the model
		call('optirun binvox', '-t {}'.format(resolution), fileName)
		voxelFile = ''.join(os.path.splitext(fileName)[:-1] + ('.binvox',))
		with open(voxelFile, 'rb') as f:
			self.voxelModel = bv.read_as_3d_array(f)
		
		self.fromBinaryArray(self.voxelModel.data)
			
	def calcFacetCounter(self):
		sumFacetsKernel = self.module.get_function("sumFacets")
		self.cF = np.zeros_like(self.cellArray)
		blockSize = (512, 1, 1)
		N = reduce(lambda x,y: x*y, self.shape)
		dimGrid = (1 + N / blockSize[0], 1)
		sumFacetsKernel(cuda.In(self.cellArray), cuda.InOut(self.cF), block=blockSize, grid=dimGrid)

			
	def markIsolatedCells(self):
		markIsolatedCellsKernel = self.module.get_function("markIsolatedCells")
		blockSize = (512, 1, 1)
		N = reduce(lambda x,y: x*y, self.shape)
		dimGrid = (1 + N / blockSize[0], 1)
		newIsolatedCells = np.zeros_like(self.cellArray)
		markIsolatedCellsKernel(cuda.In(self.cellArray), cuda.In(self.cF), cuda.InOut(self.mI), cuda.InOut(newIsolatedCells), block=blockSize, grid=dimGrid)
		return newIsolatedCells
	
	def getSimpleFacets(self):
		getSimpleFacetsKernel = self.module.get_function("getSimpleFacets")
		sFacets = np.zeros_like(self.cellArray)
		blockSize = (512, 1, 1)
		N = reduce(lambda x,y: x*y, self.shape)
		dimGrid = (1 + N / blockSize[0], 1)
		getSimpleFacetsKernel(cuda.In(self.cellArray), cuda.In(self.cF), cuda.InOut(sFacets), block=blockSize, grid=dimGrid)
		
		self.mS['dst'] = sFacets

	def getSimpleFaces(self):
		getSimpleFacesKernel = self.module.get_function("getSimpleFaces")
		sFaces = np.zeros_like(self.cellArray)
		blockSize = (512, 1, 1)
		N = reduce(lambda x,y: x*y, self.shape)
		dimGrid = (1 + N / blockSize[0], 1)
		getSimpleFacesKernel(cuda.In(self.cF), cuda.In(self.mS['dst']), 
			cuda.InOut(sFaces), block=blockSize, grid=dimGrid)
		self.mS['src'] = sFaces

	def markSimplePairs(self):
		self.getSimpleFacets()
		self.getSimpleFaces()
	
	def updateR(self):
		updateRKernel = self.module.get_function("updateR")
		blockSize = (512, 1, 1)
		N = reduce(lambda x,y: x*y, self.shape)
		dimGrid = (1 + N / blockSize[0], 1)
		updateRKernel(cuda.In(self.cellArray), cuda.InOut(self.cR), 
			block=blockSize, grid=dimGrid)
	
	def updateI(self, newIsolatedCells):
		updateIKernel = self.module.get_function("updateI")
		blockSize = (512, 1, 1)
		N = reduce(lambda x,y: x*y, self.shape)
		dimGrid = (1 + N / blockSize[0], 1)
		updateIKernel(cuda.In(newIsolatedCells), cuda.In(self.cR), 
			cuda.InOut(self.cI), block=blockSize, grid=dimGrid)
	
	def calcMedialPersistence(self):
		calcMedialPersistenceKernel = self.module.get_function("calcMedialPersistence")
		blockSize = (512, 1, 1)
		N = reduce(lambda x,y: x*y, self.shape)
		dimGrid = (1 + N / blockSize[0], 1)
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
			blockSize = (512, 1, 1)
			N = reduce(lambda x,y: x*y, self.shape)
			dimGrid = (1 + N / blockSize[0], 1)
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
			blockSize = (512, 1, 1)
			N = reduce(lambda x,y: x*y, self.shape)
			dimGrid = (1 + N / blockSize[0], 1)
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
		blockSize = (512, 1, 1)
		N = reduce(lambda x,y: x*y, self.shape)
		dimGrid = (1 + N / blockSize[0], 1)
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
Usage: cuMCThinning2D fileName -a absT -r relT -t resolution.
	fileName: Image file name.
	absT: Integer number used as absolute threshold.
	relT: Float number used as relative threshold.
	resolution: resolution of the voxelization (up to 1024)
		"""
	else:
		fileName = sys.argv[1]
		absT = None
		relT = None
		resolution = None
		for arg in sys.argv[2:]:
			cmd, val = arg.lower().split(' ')
			if cmd == '-a':
				absT = int(val)
			if cmd == '-r':
				relT = float(val)
			if cmd == '-t':
				resolution = int(val)
		
		if resolution is None:
			resolution = 128
		
		supportedImageFormats = ['.jpg', '.png']
		supportedMeshFormats = ['.stl', '.ply', '.off', '.obj']
		ext = os.path.splitext(fileName)[-1]
		if ext.lower() in supportedImageFormats:
			im = vv.imread(fileName)
			imgArray = imToGS(im)
			if absT is None:
				absT = 5 * max(imgArray.shape) / 100
				relT = 0.5
			T = (np.max(imgArray) + np.min(imgArray))/2
			imgArrayT = (imgArray < T).astype(np.uint8)
			vv.figure()
			ax1 = vv.subplot(121)
			ax2 = vv.subplot(122)
			t1 = vv.imshow(im, axes=ax1)
			t2 = vv.imshow(imgArray, axes=ax2)
			app = vv.use()
			app.Run()
			
			P = mcTSystem()
			P.fromBinaryArray(imgArrayT)
			P.doThinning(absT, relT)
			result = P.cellArray[::2,::2].astype(np.uint8)*255 + 255 - 255*imgArrayT
			Image.fromarray(result).show()
		elif ext.lower() in supportedMeshFormats:
			P = mcTSystem()
			bm = vv.meshRead(fileName)
			m = vv.mesh(bm)
			app = vv.use()
			app.Run()
			P.fromMesh(fileName, resolution)
			if absT is None:
				absT = 5 * max(P.shape) / 100
				relT = 0.5
			P.doThinning(absT, relT)
			tFileName = ''.join(os.path.splitext(fileName)[:-1] + ('-thinned.binvox'))
			tVoxelModel = binvox_rw.Voxels(
				data=P.cellArray[::2, ::2, ::2].astype(np.bool),
				dims=self.voxelModel.dims,
				translate=self.voxelModel.translate,
				scale=self.voxelModel.scale,
				axis_order=self.voxelModel.axis_order
				)
			with open(tFileName, 'wb') as f:
				tVoxelModel.write(f)
		else:
			print 'Unsupported file format.'
