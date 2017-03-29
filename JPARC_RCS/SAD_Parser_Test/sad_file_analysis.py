#!/usr/bin/env python

#---------------------------------------------------------
#This script will read SAD file and analyze it
#---------------------------------------------------------

import sys
import math
from orbit.parsers.sad_parser import SAD_Parser, SAD_LattElement, SAD_LattLine

def getLength(elem):
	L = 0.
	if(elem.hasParameter("L")):
		L = elem.getParameter("L")
	return L

if( len(sys.argv) != 2 ):
	print "Usage: >python sad_file_analysis.py <name of SAD file>"
	print "Example: >python sad_file_analysis.py ./data/rcs_lat.sad"
	sys.exit(1)

sad_file = sys.argv[1]

parser = SAD_Parser()
parser.parse(sad_file)


lines = parser.getSAD_Lines()
elems = parser.getSAD_Elements()
variables = parser.getSAD_Variables()

print "================================================"
print "The whole SAD file includes:"
print "Number of lattice accelerator lines     =",len(lines)
print "Number of lattice accelerator elements  =",len(elems)
print "Number of lattice accelerator variables =",len(variables)
print "================================================"

#get SAD lines dictionary
linesDict = parser.getSAD_LinesDict()
lines = parser.getSAD_Lines()

ring_lines = []

for ln in lines:
	elems = ln.getElements()
	length = 0.
	for elem in elems:
		length = length + getLength(elem)
	if(math.fabs(length - 348.333) < 50.):
		print "line:",ln.getName()," L=",length," nElem=",len(elems)
		ring_lines.append(ln)


print "================================================="
line = linesDict["RING"]
print "Lattice Line:",line.getName()
elems = line.getElements()
angle = 0.
for elem in elems:
	if(elem.getType() == "BEND" and elem.hasParameter("ANGLE")):
		angle = angle + elem.getParameter("ANGLE")

print "Total angle in PI = ",angle/math.pi
print "Quads statistics:"

k1_values = {}
for elem in elems:
	if(elem.getType() == "QUAD" and elem.hasParameter("K1")):
		if(elem.getParameter("K1") != 0.):
			key = math.floor(math.fabs(elem.getParameter("K1"))*1000000.)
			if(k1_values.has_key(key)):
				k1_values[key].append(elem)
			else:
				k1_values[key] = []
				k1_values[key].append(elem)
i = 0
for key in k1_values.keys():
	i = i + 1
	print "i=",i," K1=",math.fabs(k1_values[key][0].getParameter("K1")),
	print " nQuads=",len(k1_values[key])," nm=",k1_values[key][0].getName()

print "================================================="	
print "RF Cavities Positions"
s = 0.
for elem in elems:
	if(elem.hasParameter("L")):
		s = s + elem.getParameter("L")
	if(elem.getType() == "CAVI"):
		print "cav=",elem.getName()," pos=",s
	

print "Done."
sys.exit(0)
