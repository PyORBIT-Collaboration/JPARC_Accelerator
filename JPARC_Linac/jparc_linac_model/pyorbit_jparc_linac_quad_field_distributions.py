#! /usr/bin/env python

"""
This script will print out the quadrupole field gradient
distribution in accelerator sequences in the JPARC Linac.
The files will include results for the case of hard edge quads 
and the distributed quadrupole fields described by Enge's and
PMQ Trace3D functions.
"""

import sys
import math
import random
import time

from orbit.py_linac.linac_parsers import JPARC_LinacLatticeFactory

from orbit.py_linac.lattice import GetGlobalQuadGradient

from orbit.py_linac.lattice_modifications import Replace_BaseRF_Gap_to_AxisField_Nodes
from orbit.py_linac.lattice_modifications import Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes
from orbit.py_linac.lattice_modifications import Replace_Quads_to_OverlappingQuads_Nodes

from orbit.py_linac.overlapping_fields import JPARC_EngeFunctionFactory

directory_name = "./quad_fields/"

#---- accSeq_names_arr = [[names,file_out_name]]
accSeq_names_arr = []
#---- add "LI_MEBT1","LI_DTL1","LI_DTL2","LI_DTL3"
accSeq_names_arr.append([["LI_MEBT1",],directory_name+"LI_MEBT1"+"_g_fields.dat"])
accSeq_names_arr.append([["LI_DTL1",],directory_name+"LI_DTL1"+"_g_fields.dat"])
accSeq_names_arr.append([["LI_DTL2",],directory_name+"LI_DTL2"+"_g_fields.dat"])
accSeq_names_arr.append([["LI_DTL3",],directory_name+"LI_DTL3"+"_g_fields.dat"])

#---- for SDTL and ACS cavities A and B combined in one
#---- add SDTL cavities
names_tmp = []
for ind in range(16):
	seq_name = "LI_S"+"%02d"%(ind+1)
	names_tmp.append(seq_name)
accSeq_names_arr.append([names_tmp,directory_name+"SDTL"+"_g_fields.dat"])

#---- add MEBT2
accSeq_names_arr.append([["LI_MEBT2",],directory_name+"LI_MEBT2"+"_g_fields.dat"])

#---- add ACS cavities
names_tmp = []
for ind in range(21):
	seq_name = "LI_ACS"+"%02d"%(ind+1)
	names_tmp.append(seq_name)
accSeq_names_arr.append([names_tmp,directory_name+"ACS"+"_g_fields.dat"])

#---- create the factory instance
jparc_linac_factory = JPARC_LinacLatticeFactory()
jparc_linac_factory.setMaxDriftLength(0.1)

#---- the XML file name with the structure
xml_file_name = "../jparc_linac_lattice_xml/jparc_linac.xml"

for [names,fl_name_out] in accSeq_names_arr:
	#---- make lattice from XML file 
	accLattice = jparc_linac_factory.getLinacAccLattice(names,xml_file_name)
	print "Linac lattice is ready. L=",accLattice.getLength()," fl_out=",fl_name_out

	#---- magn_field_arr[[z,g0,g1],...]
	magn_field_arr = []
	step = 0.003
	n_points = int(accLattice.getLength()/step)
	step = accLattice.getLength()/(n_points - 1)
	for ip in range(n_points):
		z = step*ip
		g0 = GetGlobalQuadGradient(accLattice,z)
		g1 = 0.
		magn_field_arr.append([z,g0,g1])
	
	#------------------------------------------------------------------
	#---- Hard edge quads will be replaced by distributed fields for specified sequences 
	#------------------------------------------------------------------
	#---- longitudinal step along the distributed fields lattice
	z_step = 0.001

	Replace_Quads_to_OverlappingQuads_Nodes(accLattice,z_step,names,[],JPARC_EngeFunctionFactory)

	print "Linac new lattice is ready. L=",accLattice.getLength()
	
	for ip in range(n_points):
		[z,g0,g1] = magn_field_arr[ip]
		g1 = GetGlobalQuadGradient(accLattice,z)
		magn_field_arr[ip] = [z,g0,g1]
	
	fl_out = open(fl_name_out,"w")
	fl_out.write("z[m]        G0[T/m]       G1[T/m] \n")
	
	for ip in range(n_points):
		[z,g0,g1] = magn_field_arr[ip]
		fl_out.write(" %14.8f    %12.6f  %12.6f "%(z,g0,g1)+"\n")
		
	fl_out.close()	
	print "debug ======================================================="
	
print "Stop."
	
