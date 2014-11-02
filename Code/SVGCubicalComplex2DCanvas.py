# -*- coding: utf-8 -*-
#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  SVGCubicalComplex2dCanvas.py
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
import simpleSVG

class SVGCubical2DComplexCanvas(object):
	canvas = None
	squareSize = 20
	dotSize = 2
	arrowWidth = 6
	cubicalComplex = None
	fileName = 'temp.svg'
	_isClosed = False
	
	def __init__(self, cubicalComplex, fileName='temp.svg', squareSize=20, dotSize=2):
		self.cubicalComplex = cubicalComplex
		self.squareSize = squareSize
		self.dotSize = dotSize
		self.arrowwidth = 2 * dotSize
		self._isClosed = False
		assert(self.cubicalComplex.emb == 2)
		self.canvas = simpleSVG.svg_class(
			fname=fileName,
			bbx=self.squareSize * (1 + self.cubicalComplex.size[1]), 
			bby=self.squareSize *( 1 + self.cubicalComplex.size[0]))
		self.canvas.scale()
	
	def close(self):
		if not self._isClosed:
			self.canvas.close()
	
	def display(self):
		if not self._isClosed:
			self.close()
			self.canvas.display()
