#! /usr/bin/env python

"""
This script will calculate longitudinal acceptance for the JPARC Linac.
This example is for bunch at the start of the DTL1. If you want to calculate
the acceptance of  different parts of the linac you have to modify the
initial energy of the synchronous particle in this script.

Please note that you have to use at least BaseRfGap or RfGapTTF c++ gap models
to get reasonable results. The MatrixRfGap which is using the linear matrix for
the tracking through the RF gap will not work for this type of calculations.
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

from orbit.py_linac.lattice_modifications import Add_quad_apertures_to_lattice
from orbit.py_linac.lattice_modifications import Add_rfgap_apertures_to_lattice

from orbit.py_linac.lattice import LinacPhaseApertureNode

from orbit_utils import bunch_utils_functions
from bunch_utils_functions import copyCoordsToInitCoordsAttr
from bunch_utils_functions import swapInitCoordsAttrAndCoords

random.seed(100)

#---- Define the list of subsecuences for the lattice
names = ["LI_MEBT1","LI_DTL1","LI_DTL2","LI_DTL3"]
names = ["LI_DTL1","LI_DTL2","LI_DTL3"]

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

#======= Bunch generation as uniform rectangular region in the Longitudinal Space 
e_kin_ini = 0.003 # in [GeV]
mass = 0.939294    # in [GeV]
gamma = (mass + e_kin_ini)/mass
beta = math.sqrt(gamma*gamma - 1.0)/gamma
print "relativistic gamma=",gamma
print "relativistic  beta=",beta
frequency = 324.0e+6
v_light = 2.99792458e+8  # in [m/sec]

bunch = Bunch()
syncPart = bunch.getSyncParticle()
#set H- mass
#self.bunch.mass(0.9382723 + 2*0.000511)
bunch.mass(0.939294)
bunch.charge(-1.0)
syncPart.kinEnergy(e_kin_ini)

lambda_rf = v_light/frequency
phase_to_z_coeff = - (beta*lambda_rf)/360.

energy_spread = 0.002
nPhasePoints = 361*2
nEnergyPoints = 1000

phase_step = 360./(nPhasePoints-1)
energy_step = energy_spread/(nEnergyPoints-1)

x = 0.
xp = 0.
y = 0.
yp = 0.

for phase_ind in range(nPhasePoints):
	phase = phase_step*phase_ind - 180.
	z = phase_to_z_coeff*phase
	for energy_ind in range(nEnergyPoints):
		dE = energy_step*energy_ind - energy_spread/2
		bunch.addParticle(x,xp,y,yp,z,dE)

#==== memeorize initail coordinates
bunch.addPartAttr("ParticleInitialCoordinates")
copyCoordsToInitCoordsAttr(bunch)

print "Bunch Generation completed. nParts=",bunch.getSizeGlobal()


node_pos_dict = accLattice.getNodePositionsDict()

aprtNodes = []
for rf_gap in rf_gaps:
	phaseAperture = LinacPhaseApertureNode(frequency,rf_gap.getName()+":phaseAprt")
	pos = node_pos_dict[rf_gap][0]
	phaseAperture.setPosition(pos)
	phaseAperture.setMinMaxPhase(-180.,+180.)
	rf_gap.addChildNode(phaseAperture,AccNode.EXIT)
	aprtNodes.append(phaseAperture)
"""
for node in aprtNodes:
	print "aprt=",node.getName()," pos =",node.getPosition()
"""
print "===== Aperture Nodes Added ======= n Aprts=",len(aprtNodes)

#set up design - arrival times at each fist gaps of all RF cavities
accLattice.trackDesignBunch(bunch)

print "Design tracking completed."

#prepare to track through the lattice 
paramsDict = {"old_pos":-1.,"count":0,"pos_step":0.03}
actionContainer = AccActionsContainer("Test Design Bunch Tracking")

pos_start = 0.
print " N node position eKin Nparts "

#---- Here we define action that will analyze and print out beam parameters
#---- at the entrance and exit of each accelerator element (parents and children)

def action_entrance(paramsDict):
	node = paramsDict["node"]
	bunch = paramsDict["bunch"]
	pos = paramsDict["path_length"]
	if(paramsDict["old_pos"] == pos): return
	if(paramsDict["old_pos"] + paramsDict["pos_step"] > pos): return
	paramsDict["old_pos"] = pos
	paramsDict["count"] += 1
	nParts = bunch.getSizeGlobal()
	eKin = bunch.getSyncParticle().kinEnergy()*1.0e+3
	s_prt = " %5d  %35s  %4.5f "%(paramsDict["count"],node.getName(),pos+pos_start)
	s_prt += "  %10.6f   %8d "%(eKin,nParts)
	print s_prt
	
def action_exit(paramsDict):
	action_entrance(paramsDict)
	
#actionContainer.addAction(action_entrance, AccActionsContainer.ENTRANCE)
actionContainer.addAction(action_exit, AccActionsContainer.EXIT)

#---- This is actual tracking of the bunch
time_start = time.clock()

accLattice.trackBunch(bunch, paramsDict = paramsDict, actionContainer = actionContainer)

time_exec = time.clock() - time_start
print "time[sec]=",time_exec

#---- print out the initial phase and energy (dE) of particles
#---- that came through the lattice
file_out = open("pyorbit_init_phase_energy.dat","w")

s = "# phase[deg] energy[keV] "
file_out.write(s+"\n")

nParts = bunch.getSize()
for ind in range(nParts):
	z = bunch.partAttrValue("ParticleInitialCoordinates", ind, 4)
	dE = bunch.partAttrValue("ParticleInitialCoordinates", ind, 5)
	phase = z/phase_to_z_coeff
	dE = dE*1.0e+6
	s = "  %10.4f   %10.3f "%(phase,dE)
	file_out.write(s+"\n")

file_out.close()

#---- print real coordinates at the end of the linac
beta_new = bunch.getSyncParticle().beta()
phase_to_z_coeff_new = - (beta_new*lambda_rf)/360.
file_out = open("pyorbit_final_phase_energy.dat","w")

s = "# phase[deg] energy[keV] "
file_out.write(s+"\n")

nParts = bunch.getSize()
for ind in range(nParts):
	phase = bunch.z(ind)/phase_to_z_coeff_new
	dE = bunch.dE(ind)*1.0e+6
	s = "  %10.4f   %10.3f "%(phase,dE)
	file_out.write(s+"\n")

file_out.close()

