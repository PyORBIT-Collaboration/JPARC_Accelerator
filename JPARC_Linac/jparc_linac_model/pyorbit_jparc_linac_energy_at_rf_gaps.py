#! /usr/bin/env python

"""
This script will track the bunch through the JPARC Linac
to generate output file with the energy in and out for each 
RF gap in the lattice.
"""

import sys
import math
import random
import time

from orbit.py_linac.linac_parsers import JPARC_LinacLatticeFactory

# from linac import the C++ RF gap classes
from linac import BaseRfGap, MatrixRfGap, RfGapTTF

from bunch import Bunch
from orbit.lattice import AccLattice, AccNode, AccActionsContainer


#---- Define the list of subsecuences for the lattice
names = ["LI_MEBT1","LI_DTL1","LI_DTL2","LI_DTL3"]

#---- This part will define the sequence names 
#---- for SDTL and ACS cavities A and B combined in one

#---- add SDTL cavities
for ind in range(16):
	seq_name = "LI_S"+"%02d"%(ind+1)
	names.append(seq_name)
	
#---- add MEBT2
names.append("LI_MEBT2")

#---- add ACS cavities
for ind in range(21):
	seq_name = "LI_ACS"+"%02d"%(ind+1)
	names.append(seq_name)

"""
#---- add L3BT
names.append("LI_L3BT")
"""

#---- create the factory instance
jparc_linac_factory = JPARC_LinacLatticeFactory()

#---- the XML file name with the structure
xml_file_name = "../jparc_linac_lattice_xml/jparc_linac.xml"

#---- make lattice from XML file 
accLattice = jparc_linac_factory.getLinacAccLattice(names,xml_file_name)

seq_names = names

print "Linac lattice is ready. L=",accLattice.getLength()

#----set up RF Gap Model -------------
#---- There are three available models at this moment
#---- BaseRfGap  uses only E0TL*cos(phi)*J0(kr) with E0TL = const
#---- MatrixRfGap uses a matrix approach like envelope codes
#---- RfGapTTF uses Transit Time Factors (TTF) like PARMILA
cppGapModel = BaseRfGap
#cppGapModel = MatrixRfGap
#cppGapModel = RfGapTTF
rf_gaps = accLattice.getRF_Gaps()
for rf_gap in rf_gaps:
	rf_gap.setCppGapModel(cppGapModel())

#----- the LI_MEBT1:BNCH01 and LI_MEBT1:BNCH02 have different E0TL in the Trace3D
#----- file (see ./optics/jparc_matched2_40mA100_trace3d_rfgaps.dat file)
rf_gaps[0].setParam("E0TL",0.135408/1000.)
rf_gaps[2].setParam("E0TL",0.147204/1000.)

#set H- mass
#self.bunch.mass(0.9382723 + 2*0.000511)
bunch = Bunch()
bunch.mass(0.939294)
bunch.charge(-1.0)
syncPart = bunch.getSyncParticle()
syncPart.kinEnergy(0.003)


#prepare to track design bunch through the lattice 
rf_gaps_arr = []
paramsDict = {"rf_gaps_arr":rf_gaps_arr}
actionContainer = AccActionsContainer("Design Bunch Tracking")

def action_entrance(paramsDict):
	node = paramsDict["node"]
	bunch = paramsDict["bunch"]
	rf_gaps_arr = paramsDict["rf_gaps_arr"]
	if(node.isRFGap()):
		eKin_in = bunch.getSyncParticle().kinEnergy()*1000.
		rf_gaps_arr.append([node,eKin_in,0.])

def action_exit(paramsDict):
	node = paramsDict["node"]
	bunch = paramsDict["bunch"]
	rf_gaps_arr = paramsDict["rf_gaps_arr"]
	if(node.isRFGap()):
		n_gaps = len(rf_gaps_arr)
		rf_gaps_arr[n_gaps-1][2] = bunch.getSyncParticle().kinEnergy()*1000.

actionContainer.addAction(action_entrance, AccActionsContainer.ENTRANCE)
actionContainer.addAction(action_exit, AccActionsContainer.EXIT)

#set up design - arrival times at each fist gaps of all RF cavities
accLattice.trackDesignBunch(bunch,paramsDict = paramsDict, actionContainer = actionContainer)

print "Design tracking completed."

for [node,eKin_in,eKin_out] in rf_gaps_arr:
	print "name=",node.getName()," eKin_in,eKin_out[MeV] = %10.6f  %10.6f"%(eKin_in,eKin_out)
	
mass = bunch.mass()

file_out = open("pyorbit_jparc_rf_gaps_ekin_in_out.dat","w")
file_out.write(" rf_gap   eKinIn[MeV] eKinOut[MeV]  beta_avg \n")
for [node,eKin_in,eKin_out] in rf_gaps_arr:
	eKin_avg = ((eKin_in+eKin_out)/2.0)/1000.
	momentum = syncPart.energyToMomentum(eKin_avg)
	beta_avg = momentum/(mass+eKin_avg)
	file_out.write(node.getName()+"  "+"  %10.6f  %10.6f     %10.8f"%(eKin_in,eKin_out,beta_avg)+"\n")

file_out.close()	



