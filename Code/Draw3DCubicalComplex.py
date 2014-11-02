# -*- coding: utf-8 -*-
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
import CubicalHomology as CH
import ArrayCubicalHomology as ACH
from itertools import product
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.colors import colorConverter
from random import choice

def line3D(P, Q, axes, color='g'):
	axes.plot(*[[P[i], Q[i]] for i in xrange(3)], color='g')

def polygon3D(verts, axes, color='b', alpha=0.6):
	x = [verts[i][0] for i in xrange(len(verts))]
	y = [verts[i][1] for i in xrange(len(verts))]
	z = [verts[i][2] for i in xrange(len(verts))]
	poly = Poly3DCollection(
		[zip(x, y, z)], 
		facecolor=[colorConverter.to_rgba(color, alpha)])
	axes.add_collection3d(poly)

def homothetic(center, ratio, point):
	return tuple(np.array(point) * ratio + (1-ratio)*np.array(center))

def draw3DCube(cube, axes, color):
	r1 = 0.85
	r3 = 0.8
	a1 = 0.9
	a2 = 0.8

	if isinstance(cube, CH.Cube):
		if cube.dim == 0:
			points = map(lambda x:x[0], cube.intervals)
			print cube
			raw_input('>>')
			x, y, z = map(np.array, zip(*points))
			axes.scatter(y, z, -x, c=color, marker='o')
		elif cube.dim == 1:
			points = list(product(*map(set, cube.intervals)))
			V = [homothetic(cube.center, r1, P) for P in points]
			x, y, z = map(np.array, zip(*V))
			axes.plot(y, z, -x, color, linewidth=2, )
		elif cube.dim == 2:
			points = list(product(*map(set, cube.intervals)))
			V = [homothetic(cube.center, r1, P) for P in points]
			x, y, z = map(np.array, zip(*V))
			X, Y, Z = map(lambda x: x.copy(), [x, y, z])
			X[2], X[3] = x[3], x[2]
			Y[2], Y[3] = y[3], y[2]
			Z[2], Z[3] = z[3], z[2]
			poly = Poly3DCollection(
				[zip(Y, Z, -X)], 
				facecolor=colorConverter.to_rgba(color, a2))
			axes.add_collection3d(poly)
		elif cube.dim == 3:
			for face in cube.border():
				points = list(product(*map(set, face.intervals)))
				V = [homothetic(cube.center, r3, P) for P in points]
				x, y, z = map(np.array, zip(*V))
				X, Y, Z = map(lambda x: x.copy(), [x, y, z])
				X[2], X[3] = x[3], x[2]
				Y[2], Y[3] = y[3], y[2]
				Z[2], Z[3] = z[3], z[2]
				poly = Poly3DCollection(
					[zip(Y, Z, -X)], 
					facecolor=colorConverter.to_rgba(color, 1))
				axes.add_collection3d(poly)

def draw3DComplex(K, axes, dimensions = [0,1,2,3]):
	if isinstance(K, CH.CubicalComplex):
		if 0 in dimensions:
			points = [map(lambda x:x[0], cube.intervals) for cube in K(0) if not K.isCritical(cube)]
			if points:
				x, y, z = map(np.array, zip(*points))
				axes.scatter(y, z, -x, c='b', marker='o')
			criticalPoints = [map(lambda x:x[0], cube.intervals) for cube in K(0) if K.isCritical(cube)]
			if criticalPoints:
				x, y, z = map(np.array, zip(*criticalPoints))
				axes.scatter(y, z, -x, c='r', marker='o')
			
		r1 = 0.85
		r3 = 0.8
		a1 = 0.9
		a2 = 0.8

		if 1 in dimensions:
			for edge in K(1):
				c = 'g-' if K.isCritical(edge) else 'b--'
				points = list(product(*map(set, edge.intervals)))
				V = [homothetic(edge.center, r1, P) for P in points]
				x, y, z = map(np.array, zip(*V))
				axes.plot(y, z, -x, c, linewidth=2, )

		if 2 in dimensions:
			for face in K(2):
				points = list(product(*map(set, face.intervals)))
				V = [homothetic(face.center, r1, P) for P in points]
				x, y, z = map(np.array, zip(*V))
				X, Y, Z = map(lambda x: x.copy(), [x, y, z])
				X[2], X[3] = x[3], x[2]
				Y[2], Y[3] = y[3], y[2]
				Z[2], Z[3] = z[3], z[2]
				poly = Poly3DCollection(
					[zip(Y, Z, -X)], 
					facecolor=colorConverter.to_rgba('g', a2))
				axes.add_collection3d(poly)

		if 3 in dimensions:
			for cube in K(3):
				for face in cube.border():
					points = list(product(*map(set, face.intervals)))
					V = [homothetic(cube.center, r3, P) for P in points]
					x, y, z = map(np.array, zip(*V))
					X, Y, Z = map(lambda x: x.copy(), [x, y, z])
					X[2], X[3] = x[3], x[2]
					Y[2], Y[3] = y[3], y[2]
					Z[2], Z[3] = z[3], z[2]
					poly = Poly3DCollection(
						[zip(Y, Z, -X)], 
						facecolor=colorConverter.to_rgba('r', 1))
					axes.add_collection3d(poly)
	elif isinstance(K, ACH.CubicalComplex):
		draw3DComplex(K.toCubicalComplex(), axes, dimension)		

def draw3DVectorField(vectorField, axes, dimensions = [0,1,2]):
	if isinstance(vectorField, dict):
		for srcC in vectorField:
			if CH.Cube(srcC).dim in dimensions:
				dstC = vectorField[srcC]
				src = CH.Cube(srcC).center
				dst = CH.Cube(dstC).center
				axes.scatter(src[1], src[2], -src[0], c='r', marker='o')
				axes.scatter(dst[1], dst[2], -dst[0], c='r', marker='x')
				x, y, z = map(np.array, zip(src, dst))
				axes.plot(y, z, -x, color='r', linewidth=3)
	elif isinstance(vectorField, np.array):
		raise NotImplemented('Only clasical vector fields can be drawn at the moment.')	

if __name__ == '__main__':
	fig = plt.figure()
	ax1 = fig.add_subplot(121, projection='3d')
	plt.axis('off')
	ax2 = fig.add_subplot(122, projection='3d')
	plt.axis('off')
	
	#~ V = list(product([0,1], repeat=8))
	#~ A = np.array([1]*8, dtype=np.int8).reshape((2,2,2))
	K = CH.CubicalComplex(3)
	#~ K.fromArray(A)
	K.addCube(CH.Cube([(0,1), (0,1), (0,1)]))
	K.addCube(CH.Cube([(0,1), (1,2), (0,1)]))
	K.addCube(CH.Cube([(1,2), (0,1), (0,1)]))
	K.buildVectorField()
	draw3DComplex(K, ax1, [0,1])
	draw3DVectorField(K.vectorField, ax1)
	color = {0:'r', 1:'g-', 2:'y', 3:'r'}
	for d in K.criticalCells:
		for cube in K.criticalCells[d]:
			draw3DCube(CH.Cube(cube), ax2, color[d])
			
	plt.show()
