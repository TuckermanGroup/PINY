#!/usr/bin/env python

import piny


sim = piny.PINYEmpirical('calc-source', 'calc', 'water.input', 'piny.out')

n_atoms = sim.get_n_atoms()
m = sim.get_m()
x = sim.get_x()
sim.update()
F = -sim.get_dV_dx()
V = sim.get_V()
names = sim.get_names()

print 'number of atoms =', n_atoms
print 'atom names'
print names
print
print 'potential energy =', V
print 'masses'
print m
print
print 'old positions'
print x
print
print 'old forces'
print F
print

# change positions of first three atoms of first replica
x[:3, 0, 0] += 0.1
sim.set_x(x)

# update interactions
sim.update()

print 'new positions'
print sim.get_x()
print
print 'new forces'
print - sim.get_dV_dx()
print
