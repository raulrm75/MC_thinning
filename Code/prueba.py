#!/usr/bin/env python
# -*- coding: utf-8 -*-

from itertools import product

def F(I, N=[4]*3):
	return (I[2]-1)*(2*N[0]-1)*(2*N[1]-1) + (I[0]-1)*(2*N[1]-1) + I[1]

def f0(i, N=[4]*3):
	I = [2*j-1 for j in i]
	return F(I, N)

def f1(d, i, N=[4]*3):
	I = [2*i[j]-1 if j == d-1 else 2*i[j] for j in xrange(len(i))]
	return F(I, N)

def f2(d, i, N=[4]*3):
	I = [2*i[j]-1 if j != d-1 else 2*i[j] for j in xrange(len(i))]
	return F(I, N)

def f3s(d, i, N=[4]*3):
	I = [2*j for j in i]
	return F(I, N)

for i in product(xrange(1, 5), repeat=3):
	print i, '-->', f0(i)
