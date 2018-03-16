#! /usr/bin/env python

"""
This script will scan MEBT BNCH02 amplitude. It will track the bunch
to the end of MEBT and caculate the bunch length as a function of 
BNCH02 amplitude.
Then it will use this function and transport matrices to find the longitudinal
Twiss parameters at the entrance of the MEBT01.
"""

import sys
import math
import random
import time

from orbit.py_linac.linac_parsers import JPARC_LinacLatticeFactory

# from linac import the C++ RF gap classes
from linac import BaseRfGap, MatrixRfGap, RfGapTTF

from orbit.bunch_generators import TwissContainer
from orbit.bunch_generators import WaterBagDist3D, GaussDist3D, KVDist3D

from bunch import Bunch
from bunch import BunchTwissAnalysis

from orbit.lattice import AccLattice, AccNode, AccActionsContainer

from orbit.py_linac.lattice_modifications import Add_quad_apertures_to_lattice
from orbit.py_linac.lattice_modifications import Add_rfgap_apertures_to_lattice

from orbit.py_linac.lattice_modifications import Replace_BaseRF_Gap_to_AxisField_Nodes
from orbit.py_linac.lattice_modifications import Replace_BaseRF_Gap_and_Quads_to_Overlapping_Nodes
from orbit.py_linac.lattice_modifications import Replace_Quads_to_OverlappingQuads_Nodes

import orbit_utils
from orbit_utils import bunch_utils_functions
from bunch_utils_functions import copyCoordsToInitCoordsAttr
from bunch_utils_functions import transportMtrxFromInitCoords

from orbit_utils import Matrix, PhaseVector

from orbit.py_linac.overlapping_fields import JPARC_EngeFunctionFactory

from jparc_linac_bunch_generator import JPARC_Linac_BunchGenerator

random.seed(100)


#---- Define the list of subsecuences for the lattice
names = ["LI_MEBT1",]

#---- create the factory instance
jparc_linac_factory = JPARC_LinacLatticeFactory()
jparc_linac_factory.setMaxDriftLength(0.01)

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

node_pos_dict = accLattice.getNodePositionsDict()

#----- Read quads from the external file and set the lattice quads' fields accordingly
quad_fl_in = open("./optics/jparc_matched2_40mA100_trace3d_quads.dat","r")
lns = quad_fl_in.readlines()[1:]
quad_fl_in.close()
quad_names_field_dict = {}
for ln in lns:
	res_arr = ln.split()
	if(len(res_arr) > 2):
		name = res_arr[0].strip()
		field = float(res_arr[1])
		quad_names_field_dict[name] = field
		#print "debug trace3d quad =",name," field=",field

quads = accLattice.getQuads()
for quad in quads:
	if(quad_names_field_dict.has_key(quad.getName())):
		field = quad_names_field_dict[quad.getName()]
		quad.setParam("dB/dr",field)
		#print "debug quad =",quad.getName()," G=",quad.getParam("dB/dr")

#----- Print out quads fields and positions ==================  START
quads = accLattice.getQuads()
quad_fl_out = open("quad_params.dat","w")
quad_fl_out.write(" Quad  G[T/m]   length[m]   pos[m]    \n")
for quad in quads:
	[pos_start,pos_end] = node_pos_dict[quad]
	pos = (pos_start+pos_end)/2
	length = quad.getLength()
	field = quad.getParam("dB/dr")
	quad_fl_out.write(quad.getName()+ "  %+10.6f"%field+"  %+10.6f"%length  + "  %+12.6f"%pos  +"\n" )
quad_fl_out.close()
#----- Print out quads fields and positions ==================  END

#------------------------------------------------------------------
#---- Hard edge quads will be replaced by distributed fields for specified sequences 
#------------------------------------------------------------------
#---- longitudinal step along the distributed fields lattice
z_step = 0.001

Replace_Quads_to_OverlappingQuads_Nodes(accLattice,z_step,["LI_MEBT1",],[],JPARC_EngeFunctionFactory)
#Replace_Quads_to_OverlappingQuads_Nodes(accLattice,z_step,["LI_MEBT1","LI_DTL1"],[],JPARC_EngeFunctionFactory)
#Replace_Quads_to_OverlappingQuads_Nodes(accLattice,z_step,seq_names,[],JPARC_EngeFunctionFactory)
print "Linac new lattice is ready. L=",accLattice.getLength()
#-----------------------------------------------------
# Set up Space Charge Acc Nodes
#-----------------------------------------------------
from orbit.space_charge.sc3d import setSC3DAccNodes, setUniformEllipsesSCAccNodes
from spacecharge import SpaceChargeCalcUnifEllipse, SpaceChargeCalc3D
sc_path_length_min = 0.003

print "Set up Space Charge nodes. "

# set of uniformly charged ellipses Space Charge
nEllipses = 1
calcUnifEllips = SpaceChargeCalcUnifEllipse(nEllipses)
space_charge_nodes = setUniformEllipsesSCAccNodes(accLattice,sc_path_length_min,calcUnifEllips)

"""
# set FFT 3D Space Charge
sizeX = 64
sizeY = 64
sizeZ = 64
calc3d = SpaceChargeCalc3D(sizeX,sizeY,sizeZ)
space_charge_nodes =  setSC3DAccNodes(accLattice,sc_path_length_min,calc3d)
"""

max_sc_length = 0.
min_sc_length = accLattice.getLength()
for sc_node in space_charge_nodes:
	scL = sc_node.getLengthOfSC()
	if(scL > max_sc_length): max_sc_length = scL
	if(scL < min_sc_length): min_sc_length = scL
print "maximal SC length =",max_sc_length,"  min=",min_sc_length

print "===== Aperture Nodes START  ======="
aprtNodes = Add_quad_apertures_to_lattice(accLattice)
aprtNodes = Add_rfgap_apertures_to_lattice(accLattice,aprtNodes)

"""
for node in aprtNodes:
	print "aprt=",node.getName()," pos =",node.getPosition()
"""

print "===== Aperture Nodes Added ======="

#-----TWISS Parameters at the entrance of MEBT1 ---------------
# transverse emittances are unnormalized and in pi*mm*mrad
# longitudinal emittance is in pi*eV*sec
e_kin_ini = 0.003 # in [GeV]
mass = 0.939294    # in [GeV]
gamma = (mass + e_kin_ini)/mass
beta = math.sqrt(gamma*gamma - 1.0)/gamma
print "relat. gamma=",gamma
print "relat.  beta=",beta
frequency = 324.0e+6
v_light = 2.99792458e+8  # in [m/sec]

#-----------------------------------------------------------------
#------ Twiss parameters from Trace3D 
#------ x,x'  y,y'  mm,mrad 
#------ Long Twiss: emittance deg*keV , beta deg/keV 
#------ All emittances are 5 times bigger the RMS values
#------ The Trace3D file is 
#------ ../jparc_linac_xml_generator/jparc_trace3d_lattice/jparc_matched2_40mA100.t3d

(alphaX,betaX,emittX) =  (-2.153,   0.18313, (13.4606/5)*1.0E-6)
(alphaY,betaY,emittY) =  ( 1.8253,  0.15578, (13.3520/5)*1.0E-6)
(alphaZ,betaZ,emittZ) =  (-0.08328, 0.91909, (538.69/5))

"""
 From Trace3D input file ()
 === for 40 mA =====
 ER= 939.3014, Q= -1.0, W =  3.0, XI=  50.0,XI=15,XI=40
 EMITI =     13.460568358239064 13.352006807621624 538.68785
 BEAMI =    -2.153	.18313	1.8253	.15578	-.08328	.9190900000000001

"""
#---- Translate long. Twiss from Trace3D to PyORBIT (emittance m*GeV , beta m/GeV)
betaZ  = betaZ *(v_light*beta/(360.*frequency))*1.0E+6
emittZ = emittZ*(v_light*beta/(360.*frequency))/1.0E+6

print " ========= PyORBIT Twiss ==========="
print " aplha beta emitt[mm*mrad] X= %6.4f %6.4f %6.4f "%(alphaX,betaX,emittX*1.0e+6)
print " aplha beta emitt[mm*mrad] Y= %6.4f %6.4f %6.4f "%(alphaY,betaY,emittY*1.0e+6)
print " aplha beta emitt[mm*MeV] Z= %6.4f %6.4f %6.4f "%(alphaZ,betaZ,emittZ*1.0e+6)

twissX = TwissContainer(alphaX,betaX,emittX)
twissY = TwissContainer(alphaY,betaY,emittY)
twissZ = TwissContainer(alphaZ,betaZ,emittZ)

print "Start Bunch Generation."
bunch_gen = JPARC_Linac_BunchGenerator(twissX,twissY,twissZ)

#set the initial kinetic energy in GeV
bunch_gen.setKinEnergy(e_kin_ini)

#set the beam peak current in mA
bunch_gen.setBeamCurrent(40.0)

#bunch_in = bunch_gen.getBunch(nParticles = 10000, distributorClass = WaterBagDist3D)
bunch_in = bunch_gen.getBunch(nParticles = 10000, distributorClass = GaussDist3D)
#bunch_in = bunch_gen.getBunch(nParticles = 10000, distributorClass = KVDist3D)

print "Bunch Generation completed."

bnch02 = accLattice.getNodeForName("LI_MEBT1:BNCH02:Rg01")
bnch02_ind = accLattice.getNodeIndex(bnch02)

#set up design - arrival times at each fist gaps of all RF cavities
accLattice.trackDesignBunch(bunch_in)
print "Design tracking completed."

#real tracking of the bunch. Now we have the bunch at the BNCH02 entrance
accLattice.trackBunch(bunch_in,None,None,0,bnch02_ind-1)

#---- memorize initial coordinates of particles in PartAttributes
copyCoordsToInitCoordsAttr(bunch_in)

twiss_analysis = BunchTwissAnalysis()

bunch = bunch_in
twiss_analysis.analyzeBunch(bunch)
x_rms = math.sqrt(twiss_analysis.getTwiss(0)[1]*twiss_analysis.getTwiss(0)[3])*1000.
y_rms = math.sqrt(twiss_analysis.getTwiss(1)[1]*twiss_analysis.getTwiss(1)[3])*1000.
z_rms = math.sqrt(twiss_analysis.getTwiss(2)[1]*twiss_analysis.getTwiss(2)[3])*1000.
(alphaZ,betaZ,emittZ) = (twiss_analysis.getTwiss(2)[0],twiss_analysis.getTwiss(2)[1],twiss_analysis.getTwiss(2)[3])
gammaZ = (1.0+alphaZ**2)/betaZ
dE_rms = math.sqrt(gammaZ*emittZ*1.0e+6)*1000.
z_to_phase_coeff = bunch_gen.getZtoPhaseCoeff(bunch)
z_rms_deg = z_to_phase_coeff*z_rms/1000.0
print "========== Bunch RMS sizes at BNCH02 entrance:"
print "(Sx[mm],Sy[mm],Sz[deg]) = %5.3f  %5.3f   %5.3f "%(x_rms,y_rms,z_rms_deg)
print "(alphaZ, betaZ[deg/MeV], emittZ[deg*keV] = (%5.4f  %5.4f  %12.5g )"%(alphaZ,betaZ*z_to_phase_coeff/1000.,emittZ*z_to_phase_coeff*1000.*1000.)
print " z_rms[deg] = %5.2f "%(z_rms_deg)
print " corr_z_dE  = %5.2f "%(-alphaZ*emittZ*z_to_phase_coeff*1000.*1000.)
print " dE_rms[keV]= %5.2f "%(dE_rms)

#==== Now we track the bunch from BNCH02 etrance to the end of the MEBT
bunch = Bunch()
bunch_in.copyBunchTo(bunch)
accLattice.trackBunch(bunch,None,None,bnch02_ind)

twiss_analysis.analyzeBunch(bunch)
x_rms = math.sqrt(twiss_analysis.getTwiss(0)[1]*twiss_analysis.getTwiss(0)[3])*1000.
y_rms = math.sqrt(twiss_analysis.getTwiss(1)[1]*twiss_analysis.getTwiss(1)[3])*1000.
z_rms = math.sqrt(twiss_analysis.getTwiss(2)[1]*twiss_analysis.getTwiss(2)[3])*1000.
z_to_phase_coeff = bunch_gen.getZtoPhaseCoeff(bunch)
z_rms_deg = z_to_phase_coeff*z_rms/1000.0
print "Bunch RMS sizes at the end of MEBT:"
print "(Sx,Sy,Sz) = %5.3f  %5.3f   %5.3f "%(x_rms,y_rms,z_rms_deg)

#==== Now we repeat tracking bunch from BNCH02 to the end of the MEBT
#==== for different values of the E0TL 

#coeff_E0TL_to_E0L = bnch02.getParam("E0L")/bnch02.getParam("E0TL")


E0L_E0TL_arr = []
n_points = 10
amp_coeff = 1.5
E0TL_max = bnch02.getParam("E0TL")*amp_coeff 
E0L_max = bnch02.getParam("E0L")*amp_coeff

E0TL_min = bnch02.getParam("E0TL")*0.3


for ind in range(n_points):
	E0TL = E0TL_max*(1.0*(ind+1))/n_points
	E0L = E0L_max*(1.0*(ind+1))/n_points
	if(E0TL > E0TL_min):
		E0L_E0TL_arr.append([E0L,E0TL])
"""
E0L_E0TL_arr = []
E0L_E0TL_arr.append([1.,0.138/1000.])
E0L_E0TL_arr.append([1.,0.154/1000.])
E0L_E0TL_arr.append([1.,0.170/1000.])
E0L_E0TL_arr.append([1.,0.178/1000.])
E0L_E0TL_arr.append([1.,0.186/1000.])
E0L_E0TL_arr.append([1.,0.196/1000.])

for arr in E0L_E0TL_arr:
	arr[0] = arr[1]*coeff_E0TL_to_E0L
	
z_rms_expr_arr = [ 13.95859799,12.80266955,11.30070564,10.92338666,10.59683631,10.79909468]

#E0TL in MV
#[ 0.13816452,  0.154     ,  0.16983548,  0.17775321,  0.18567095,0.19358869]
#
#dphi_rms in deg. (324MHz)
#[ 13.95859799,12.80266955,11.30070564,10.92338666,10.59683631,10.79909468]
"""
n_points = len(E0L_E0TL_arr)	

use_twiss_weight_x = True
use_twiss_weight_y = True
use_twiss_weight_z = True

res_sizes_trMtrx_arr = []

count = 0
for [E0L,E0TL] in E0L_E0TL_arr:
	trMtrx = Matrix(7,7)
	bnch02.setParam("E0TL",E0TL)
	bnch02.setParam("E0L",E0L)
	bunch = Bunch()
	bunch_in.copyBunchTo(bunch)
	accLattice.trackBunch(bunch,None,None,bnch02_ind)
	transportMtrxFromInitCoords(bunch,trMtrx,use_twiss_weight_x, use_twiss_weight_y, use_twiss_weight_z)	
	twiss_analysis.analyzeBunch(bunch)
	x_rms = math.sqrt(twiss_analysis.getTwiss(0)[1]*twiss_analysis.getTwiss(0)[3])*1000.
	y_rms = math.sqrt(twiss_analysis.getTwiss(1)[1]*twiss_analysis.getTwiss(1)[3])*1000.
	z_rms = math.sqrt(twiss_analysis.getTwiss(2)[1]*twiss_analysis.getTwiss(2)[3])*1000.
	z_to_phase_coeff = bunch_gen.getZtoPhaseCoeff(bunch)
	z_rms_deg = z_to_phase_coeff*z_rms/1000.0	
	#----transformation from [m] to [deg] and from [GeV] to [MeV] for z, phi and dE
	trMtrx.set(4,4+1,trMtrx.get(4,4+1)*z_to_phase_coeff/1000.)
	trMtrx.set(4+1,4,trMtrx.get(4+1,4)*1000./z_to_phase_coeff)
	mz11 = trMtrx.get(4,4)
	mz12 = trMtrx.get(4,4+1)
	#----------set experimental length----------------
	#z_rms_deg = z_rms_expr_arr[count]
	#----------------------------------------------------------
	det_z = trMtrx.get(0+4,0+4)*trMtrx.get(1+4,1+4) - trMtrx.get(1+4,0+4)*trMtrx.get(0+4,1+4)	
	res_sizes_trMtrx_arr.append([E0L,E0TL,(x_rms,y_rms,z_rms_deg),mz11,mz12])
	print "E0TL[keV] =  %5.1f (Sx,Sy,Sz) =   %5.3f  %5.3f   %5.1f  det(TrMtrxZ) = %5.4f "%(E0TL*1000.*1000.,x_rms,y_rms,z_rms_deg,det_z)
	#-------------------
	count += 1


#--------------------------------------------------------------
#----------- LSQ method for Longitudinal Twiss parameters
#--------------------------------------------------------------
n_cases = len(res_sizes_trMtrx_arr)
matrx_M = Matrix(n_cases,3)
z_rms_2_vctr = PhaseVector(n_cases)

line_ind = 0
for [E0L,E0TL,(x_rms,y_rms,z_rms_deg),mz11,mz12] in res_sizes_trMtrx_arr:
	#print "E0TL[keV] =  %5.1f (Sx,Sy,Sz) =   %5.3f  %5.3f   %5.1f "%(E0TL*1000.*1000.,x_rms,y_rms,z_rms_deg)
	matrx_M.set(line_ind,0, mz11**2)
	matrx_M.set(line_ind,1, 2*mz11*mz12)
	matrx_M.set(line_ind,2, mz12**2)
	#----------------------------
	z_rms_2_vctr.set(line_ind,z_rms_deg**2)
	#----------------------------
	line_ind += 1
	
#---- error matrix   sigma_z in deg and sigma_z_rel - relative
sigma_z = 1.0
sigma_z_rel = 0.05
print "======== Error of longitudinal rms size in deg from DTL scans ===="
print "z sigma[deg] = ",sigma_z
print "z sigma relative[%] = ",sigma_z_rel*100.
print "==================================================================="
matrx_W = Matrix(n_cases,n_cases)
for ind0 in range(n_cases):
	for ind1 in range(n_cases):
		matrx_W.set(ind0,ind1,0.)
		if(ind0 == ind1):
			matrx_W.set(ind0,ind1,1./(2*z_rms_2_vctr.get(ind0)*sigma_z_rel)**2)
			#matrx_W.set(ind0,ind1,1./(2*math.sqrt(z_rms_2_vctr.get(ind0))*sigma_z)**2)
			
#---- (MT*W*M)^(-1)
matrx_MT = Matrix(matrx_M)
matrx_MT.transpose()
MTWM_m1 = ((matrx_MT.mult(matrx_W)).mult(matrx_M)).invert()


long_corr_results_vector = MTWM_m1.mult(matrx_MT.mult(matrx_W.mult(z_rms_2_vctr)))
phase2_rms_res = long_corr_results_vector.get(0)
phase_ekin_corr_res = long_corr_results_vector.get(1)
ekin2_rms_res = long_corr_results_vector.get(2)
phase_rms_res = math.sqrt(phase2_rms_res)
ekin_rms_res = math.sqrt(ekin2_rms_res)

#---- Error estimation
phase_rms_err = math.sqrt(MTWM_m1.get(0,0))/(2*phase_rms_res)
phase_ekin_corr_err = math.sqrt(MTWM_m1.get(1,1))
ekin_rms_err = math.sqrt(MTWM_m1.get(2,2))/(2*ekin_rms_res)

print "================BNCH02 entrance================"
print "z_rms[deg]            = %6.2f +- %6.2f "%(phase_rms_res,phase_rms_err)
print "corr(phi,dE)[deg*keV] = %6.2f +- %6.2f"%(phase_ekin_corr_res*1000.,phase_ekin_corr_err*1000.)
print "dE_rms[keV]           = %6.2f +- %6.2f"%(ekin_rms_res*1000.,ekin_rms_err*1000.)

#========================================================
emittance = math.sqrt(phase2_rms_res*ekin2_rms_res - phase_ekin_corr_res**2)

emitt_err2  =  ((phase_rms_res*ekin2_rms_res/emittance)*phase_rms_err)**2
emitt_err2 +=  ((ekin_rms_res*phase2_rms_res/emittance)*ekin_rms_err)**2
emitt_err2 +=  ((phase_ekin_corr_res/emittance)*phase_ekin_corr_err)**2
emitt_err = math.sqrt(emitt_err2)

print "emitt[deg*keV]        = %6.2f +- %6.2f"%(emittance*1000.,emitt_err*1000.)



