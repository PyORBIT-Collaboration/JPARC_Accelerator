#! /usr/bin/env python

"""
This script will scan MEBT BNCH02 amplitude. It will track the bunch
to the end of MEBT and caculate the bunch length as a function of 
BNCH02 amplitude.
It will try to fit longitudinal parameters at the BNCH02 entrance 
to reproduce the measured bunch length for different E0TL.
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

#Replace_Quads_to_OverlappingQuads_Nodes(accLattice,z_step,["LI_MEBT1",],[],JPARC_EngeFunctionFactory)
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
(alphaX,betaX,emittX) = (twiss_analysis.getTwiss(0)[0],twiss_analysis.getTwiss(0)[1],twiss_analysis.getTwiss(0)[3])
(alphaY,betaY,emittY) = (twiss_analysis.getTwiss(1)[0],twiss_analysis.getTwiss(1)[1],twiss_analysis.getTwiss(1)[3])
(alphaZ,betaZ,emittZ) = (twiss_analysis.getTwiss(2)[0],twiss_analysis.getTwiss(2)[1],twiss_analysis.getTwiss(2)[3])
gammaZ = (1.0+alphaZ**2)/betaZ
dE_rms = math.sqrt(gammaZ*emittZ*1.0e+6)*1000.
z_to_phase_coeff = bunch_gen.getZtoPhaseCoeff(bunch)
z_rms_deg = z_to_phase_coeff*z_rms/1000.0
print "========== Bunch RMS sizes at BNCH02 entrance:"
print "(Sx[mm],Sy[mm],Sz[deg]) = %5.3f  %5.3f   %5.3f "%(x_rms,y_rms,z_rms_deg)
print "(alphaZ, betaZ[mm/MeV], emittZ[mm*MeV] = (%5.4f  %12.5g  %12.5g )"%(alphaZ,betaZ,emittZ*1.0e+6)
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
#==== for different values of the E0TL that were used in the measurements.
#==== z_rms[deg] measured bunch lengths from DTL acceptace scans

#==== E0TL_arr[[E0TL[GeV] , z_rms[deg]]
E0TL_Z_RMS_arr = []
E0TL_Z_RMS_arr.append([0.138/1000.,13.96])
E0TL_Z_RMS_arr.append([0.154/1000.,12.80])
E0TL_Z_RMS_arr.append([0.170/1000.,11.30])
E0TL_Z_RMS_arr.append([0.178/1000.,10.92])
E0TL_Z_RMS_arr.append([0.186/1000.,10.60])
E0TL_Z_RMS_arr.append([0.196/1000.,10.80])

n_points = len(E0TL_Z_RMS_arr)	

use_twiss_weight_x = True
use_twiss_weight_y = True
use_twiss_weight_z = True

res_sizes_trMtrx_arr = []

for [E0TL,z_rms_exp] in E0TL_Z_RMS_arr:
	trMtrx = Matrix(7,7)
	bnch02.setParam("E0TL",E0TL)
	bnch02.setParam("E0L",E0TL)
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
	#----------------------------------------------------------
	det_z = trMtrx.get(0+4,0+4)*trMtrx.get(1+4,1+4) - trMtrx.get(1+4,0+4)*trMtrx.get(0+4,1+4)	
	res_sizes_trMtrx_arr.append([E0TL,z_rms_deg,z_rms_exp,mz11,mz12])
	print "E0TL[keV] =  %5.1f (z_rms,z_rms_exp) =  %5.3f  %5.1f    det(TrMtrxZ) = %5.4f "%(E0TL*1.0e+6,z_rms_deg,z_rms_exp,det_z)
	#-------------------



#--------------------------------------------------------------
#----------- LSQ method for Longitudinal Twiss parameters
#--------------------------------------------------------------
n_cases = len(res_sizes_trMtrx_arr)
matrx_M = Matrix(n_cases,3)
z_rms_2_vctr = PhaseVector(n_cases)

line_ind = 0
for [E0TL,z_rms_deg,z_rms_exp,mz11,mz12] in res_sizes_trMtrx_arr:
	matrx_M.set(line_ind,0, mz11**2)
	matrx_M.set(line_ind,1, 2*mz11*mz12)
	matrx_M.set(line_ind,2, mz12**2)
	#----------------------------
	z_rms_2_vctr.set(line_ind,z_rms_exp**2)
	#----------------------------
	line_ind += 1
	
#---- error matrix   sigma_z in deg and sigma_z_rel - relative
sigma_z = 0.1
#sigma_z_rel = 0.05
print "======== Error of longitudinal rms size in deg from DTL scans ===="
print "z sigma[deg] = ",sigma_z
#print "z sigma relative[%] = ",sigma_z_rel*100.
print "==================================================================="
matrx_W = Matrix(n_cases,n_cases)
for ind0 in range(n_cases):
	for ind1 in range(n_cases):
		matrx_W.set(ind0,ind1,0.)
		if(ind0 == ind1):
			#matrx_W.set(ind0,ind1,1./(2*z_rms_2_vctr.get(ind0)*sigma_z_rel)**2)
			matrx_W.set(ind0,ind1,1./(2*math.sqrt(z_rms_2_vctr.get(ind0))*sigma_z)**2)
			
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

print "================Zero Iteration for Long. Twiss at BNCH02 entrance================"
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
print "================================================================================="
#=================================================================
# The first step is finished. We found the long. Twiss with the transport
# matrices generated by the design Twiss, and used these matrices to solve
# the LSQ problem for the initial Twiss by using the real length data.
#=================================================================

#=================================================================
# Now we will solve the self-consistent problem: by fitting we will
# find the initial Twiss reproducing the real length data.
# In this case the transport matrices will be generated on fly.
#=================================================================
(alphaZ,betaZ,emittZ) = (-phase_ekin_corr_res/emittance,phase2_rms_res/emittance,emittance)
(alphaZ,betaZ,emittZ) = (alphaZ,(1000./z_to_phase_coeff)*betaZ,emittZ/(z_to_phase_coeff*1000))
print "Zero Iter. Twiss: (alphaZ, betaZ[mm/MeV], emittZ[mm*MeV] = (%5.4f  %12.5g  %12.5g )"%(alphaZ,betaZ,emittZ*1.0e+6)

twissX = TwissContainer(alphaX,betaX,emittX)
twissY = TwissContainer(alphaY,betaY,emittY)
twissZ = TwissContainer(alphaZ,betaZ,emittZ)

class FitFunction:
	def __init__(self, bunch_gen, twissX, twissY, twissZ , E0TL_Z_RMS_arr):
		self.E0TL_Z_RMS_arr = E0TL_Z_RMS_arr
		self.bunch_gen = bunch_gen
		self.twissX = twissX
		self.twissY = twissY
		self.twissZ = twissZ
		self.min_twissZ = None
		self.diff2_min = 1.0e+46
		self.n_iter = 0
		
	def getDiff2(self,params_arr, print_info = True):
		"""
		params_arr = [alphaZ,betaZ,emittZ]
		"""
		(alphaZ,betaZ,emittZ) = params_arr[0:3]
		self.twissZ = TwissContainer(alphaZ,betaZ,emittZ)
		self.bunch_gen.twiss = (self.twissX,self.twissY,self.twissZ)
		#bunch_in = self.bunch_gen.getBunch(nParticles = 10000, distributorClass = WaterBagDist3D)
		bunch_in = self.bunch_gen.getBunch(nParticles = 10000, distributorClass = GaussDist3D)
		accLattice.trackDesignBunch(bunch_in,None,None,bnch02_ind)
		self.n_iter += 1
		diff2 = 0.
		for [E0TL,z_rms_exp] in self.E0TL_Z_RMS_arr:
			bnch02.setParam("E0TL",E0TL)
			bnch02.setParam("E0L",E0TL)
			bunch = Bunch()
			bunch_in.copyBunchTo(bunch)
			accLattice.trackBunch(bunch,None,None,bnch02_ind)
			twiss_analysis.analyzeBunch(bunch)
			z_rms = math.sqrt(twiss_analysis.getTwiss(2)[1]*twiss_analysis.getTwiss(2)[3])*1000.
			z_to_phase_coeff = bunch_gen.getZtoPhaseCoeff(bunch)
			z_rms_deg = z_to_phase_coeff*z_rms/1000.0	
			diff2 += (z_rms_deg-z_rms_exp)**2
		diff2 /= len(self.E0TL_Z_RMS_arr)
		if(diff2 < self.diff2_min):
			if(print_info == True):
				print "==== new min iter =",self.n_iter," diff = ",math.sqrt(diff2)
				print "  (alphaZ,betaZ,emittZ) = ",(alphaZ,betaZ,emittZ*1.0e+6)			
			self.diff2_min = diff2
			self.min_twissZ = self.twissZ
		return diff2
		
	def printDataComparison(self):
		self.twissZ = self.min_twissZ
		self.bunch_gen.twiss = (self.twissX,self.twissY,self.twissZ)
		#bunch_in = self.bunch_gen.getBunch(nParticles = 10000, distributorClass = WaterBagDist3D)
		bunch_in = self.bunch_gen.getBunch(nParticles = 10000, distributorClass = GaussDist3D)
		accLattice.trackDesignBunch(bunch_in,None,None,bnch02_ind)
		self.n_iter += 1
		diff2 = 0.
		for [E0TL,z_rms_exp] in self.E0TL_Z_RMS_arr:
			bnch02.setParam("E0TL",E0TL)
			bnch02.setParam("E0L",E0TL)
			bunch = Bunch()
			bunch_in.copyBunchTo(bunch)
			accLattice.trackBunch(bunch,None,None,bnch02_ind)
			twiss_analysis.analyzeBunch(bunch)
			z_rms = math.sqrt(twiss_analysis.getTwiss(2)[1]*twiss_analysis.getTwiss(2)[3])*1000.
			z_to_phase_coeff = bunch_gen.getZtoPhaseCoeff(bunch)
			z_rms_deg = z_to_phase_coeff*z_rms/1000.0
			print "E0TL[keV] = %6.2f "%(E0TL*1.0e+6),"  z_rms model / experiment [deg] = %6.2f / %6.2f "%(z_rms_deg,z_rms_exp)
	

fitFunctionObj = 	FitFunction(bunch_gen, twissX, twissY, twissZ , E0TL_Z_RMS_arr)

fitFunction = fitFunctionObj.getDiff2

(alphaZ,betaZ,emittZ) =  (-1.5207884851831515, 200.170606305362, 0.03665142894049984*1.0e-6)

fit_params_arr = []
fit_params_arr.append(alphaZ)
fit_params_arr.append(betaZ)
fit_params_arr.append(emittZ)

fit_params_step_arr = []
fit_params_step_arr.append(alphaZ*0.05)
fit_params_step_arr.append(betaZ*0.05)
fit_params_step_arr.append(emittZ*0.05)

print "INITIAL diff=",math.sqrt(fitFunctionObj.getDiff2(fit_params_arr,False))
print "INITIAL params =",fit_params_arr

"""
from Simplex import Simplex
simplex = Simplex(fitFunction,fit_params_arr,fit_params_step_arr)

(fit_params_arr, err2, nIter) =  simplex.minimize(epsilon = 0.000000001, maxiters = 100, monitor = 0)
print "total iteration =",nIter," err=",math.sqrt(err2)
"""

print "Min TwissZ (alpha,beta,gamma,emitt) = ",fitFunctionObj.min_twissZ.getAlphaBetaGammaEmitt()

fitFunctionObj.printDataComparison()

#==============================================================================
print "===========Let's do the matrix analysis again!         =========="
print "===========now with the long. Twiss found from fittting=========="
#==============================================================================

bunch_in = bunch_gen.getBunch(nParticles = 10000, distributorClass = GaussDist3D)

twiss_analysis.analyzeBunch(bunch_in)

(alphaZ,betaZ,gammaZ,emittZ) = fitFunctionObj.min_twissZ.getAlphaBetaGammaEmitt()

z_rms = math.sqrt(twiss_analysis.getTwiss(2)[1]*twiss_analysis.getTwiss(2)[3])*1000.
gammaZ = (1.0+alphaZ**2)/betaZ
dE_rms = math.sqrt(gammaZ*emittZ*1.0e+6)*1000.
z_rms_deg = z_to_phase_coeff*z_rms/1000.0
print "========== Bunch RMS sizes at BNCH02 entrance:"
print "(alphaZ, betaZ[mm/MeV], emittZ[mm*MeV] = (%5.4f  %12.5g  %12.5g )"%(alphaZ,betaZ,emittZ*1.0e+6)
print "(alphaZ, betaZ[deg/MeV], emittZ[deg*keV] = (%5.4f  %5.4f  %12.5g )"%(alphaZ,betaZ*z_to_phase_coeff/1000.,emittZ*z_to_phase_coeff*1000.*1000.)
print " z_rms[deg] = %5.2f "%(z_rms_deg)
print " corr_z_dE  = %5.2f "%(-alphaZ*emittZ*z_to_phase_coeff*1000.*1000.)
print " dE_rms[keV]= %5.2f "%(dE_rms)

#---- memorize initial coordinates of particles in PartAttributes
copyCoordsToInitCoordsAttr(bunch_in)

accLattice.trackDesignBunch(bunch_in,None,None,bnch02_ind)

res_sizes_trMtrx_arr = []

for [E0TL,z_rms_exp] in E0TL_Z_RMS_arr:
	trMtrx = Matrix(7,7)
	bnch02.setParam("E0TL",E0TL)
	bnch02.setParam("E0L",E0TL)
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
	#----------------------------------------------------------
	det_z = trMtrx.get(0+4,0+4)*trMtrx.get(1+4,1+4) - trMtrx.get(1+4,0+4)*trMtrx.get(0+4,1+4)	
	res_sizes_trMtrx_arr.append([E0TL,z_rms_deg,z_rms_exp,mz11,mz12])
	print "E0TL[keV] =  %5.1f (z_rms,z_rms_exp) =  %5.3f  %5.1f    det(TrMtrxZ) = %5.4f "%(E0TL*1.0e+6,z_rms_deg,z_rms_exp,det_z)
	#-------------------



#--------------------------------------------------------------
#----------- LSQ method for Longitudinal Twiss parameters
#--------------------------------------------------------------
n_cases = len(res_sizes_trMtrx_arr)
matrx_M = Matrix(n_cases,3)
z_rms_2_vctr = PhaseVector(n_cases)

line_ind = 0
for [E0TL,z_rms_deg,z_rms_exp,mz11,mz12] in res_sizes_trMtrx_arr:
	matrx_M.set(line_ind,0, mz11**2)
	matrx_M.set(line_ind,1, 2*mz11*mz12)
	matrx_M.set(line_ind,2, mz12**2)
	#----------------------------
	z_rms_2_vctr.set(line_ind,z_rms_exp**2)
	#----------------------------
	line_ind += 1
	
#---- error matrix   sigma_z in deg and sigma_z_rel - relative
sigma_z = 0.1
#sigma_z_rel = 0.05
print "======== Error of longitudinal rms size in deg from DTL scans ===="
print "z sigma[deg] = ",sigma_z
#print "z sigma relative[%] = ",sigma_z_rel*100.
print "==================================================================="
matrx_W = Matrix(n_cases,n_cases)
for ind0 in range(n_cases):
	for ind1 in range(n_cases):
		matrx_W.set(ind0,ind1,0.)
		if(ind0 == ind1):
			#matrx_W.set(ind0,ind1,1./(2*z_rms_2_vctr.get(ind0)*sigma_z_rel)**2)
			matrx_W.set(ind0,ind1,1./(2*math.sqrt(z_rms_2_vctr.get(ind0))*sigma_z)**2)
			
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

print "================Zero Iteration for Long. Twiss at BNCH02 entrance================"
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
print "================================================================================="

#===========================================================================
#  Now we will track bunch from MEBT02 to the start of MEBT 
#           by using backward tracking
# 1. Step: Reverse lattice and track reverse bunch from BNCH02(not included) to end
# 2. Reverse resulted bunch. It is initial bunch at the entrance of the MEBT
# 3. Reverse lattice back to original order
# 4. Track initial bunch from start to end of MEBT for experimental values 
#    of the E0TL for BNCH02
#===========================================================================
accLattice.reverseOrder()

bnch02 = accLattice.getNodeForName("LI_MEBT1:BNCH02:Rg01")
bnch02_ind = accLattice.getNodeIndex(bnch02)

def BunchTransformerFunc(bunch):
	""" 
	This function will reverse all xp, yp, z coordinates of the bunch.
	We have to change the sign of the z because the tail will be the head
	of the bunch, but the sign of dE will not change because of the same reason. 
	"""
	nParts = bunch.getSize()
	for i in range(nParts):
		(xp,yp,z,dE) = (bunch.xp(i),bunch.yp(i),bunch.z(i),bunch.dE(i))
		bunch.xp(i,-xp)
		bunch.yp(i,-yp)
		bunch.z(i,-z)
		#--- dE should not change the sign
		#bunch.dE(i,-dE)

bunch = Bunch()
bunch_in.copyBunchTo(bunch)

BunchTransformerFunc(bunch)
bunch.getSyncParticle().time(0.)

accLattice.trackDesignBunch(bunch,None,None,bnch02_ind+1)
accLattice.trackBunch(bunch,None,None,bnch02_ind+1)

BunchTransformerFunc(bunch)
bunch.getSyncParticle().time(0.)

bunch_in = bunch


twiss_analysis.analyzeBunch(bunch)
(alphaX,betaX,emittX) = (twiss_analysis.getTwiss(0)[0],twiss_analysis.getTwiss(0)[1],twiss_analysis.getTwiss(0)[3]*1.0e+6)
(alphaY,betaY,emittY) = (twiss_analysis.getTwiss(1)[0],twiss_analysis.getTwiss(1)[1],twiss_analysis.getTwiss(1)[3]*1.0e+6)
(alphaZ,betaZ,emittZ) = (twiss_analysis.getTwiss(2)[0],twiss_analysis.getTwiss(2)[1],twiss_analysis.getTwiss(2)[3]*1.0e+6)
print "PyORBIT initial Twiss at the MEBT entrance:"
s = ""
s += "(alphaX,betaX,emittX) = (%+6.4f , %+6.4f , %+6.4f)  \n "%(alphaX,betaX,emittX)
s += "(alphaY,betaY,emittY) = (%+6.4f , %+6.4f , %+6.4f)  \n "%(alphaY,betaY,emittY)
s += "(alphaZ,betaZ,emittZ) = (%+6.4f , %+6.4f , %+6.4f)  \n "%(alphaZ,betaZ,emittZ)

accLattice.reverseOrder()

bnch02 = accLattice.getNodeForName("LI_MEBT1:BNCH02:Rg01")

accLattice.trackDesignBunch(bunch_in)

print "Tracking the initial bunch through the whole MEBT:"
for [E0TL,z_rms_exp] in E0TL_Z_RMS_arr:
	bnch02.setParam("E0TL",E0TL)
	bnch02.setParam("E0L",E0TL)
	bunch = Bunch()
	bunch_in.copyBunchTo(bunch)
	accLattice.trackBunch(bunch)	
	twiss_analysis.analyzeBunch(bunch)
	x_rms = math.sqrt(twiss_analysis.getTwiss(0)[1]*twiss_analysis.getTwiss(0)[3])*1000.
	y_rms = math.sqrt(twiss_analysis.getTwiss(1)[1]*twiss_analysis.getTwiss(1)[3])*1000.
	z_rms = math.sqrt(twiss_analysis.getTwiss(2)[1]*twiss_analysis.getTwiss(2)[3])*1000.
	z_to_phase_coeff = bunch_gen.getZtoPhaseCoeff(bunch)
	z_rms_deg = z_to_phase_coeff*z_rms/1000.0	
	#----------------------------------------------------------
	print "E0TL[keV] =  %5.1f (z_rms,z_rms_exp) =  %5.3f  %5.1f "%(E0TL*1.0e+6,z_rms_deg,z_rms_exp)
	#----------------------------------------------------------
	
#================================================================================
#  Calculate Twiss and sizes through the whole MEBT for the production conditions
#================================================================================
bnch02.setParam("E0TL",0.147204/1000.)
bunch = Bunch()
bunch_in.copyBunchTo(bunch)


#prepare to track through the lattice 
paramsDict = {"old_pos":-1.,"count":0,"pos_step":0.01}
actionContainer = AccActionsContainer("Bunch Tracking")

pos_start = 0.

twiss_analysis = BunchTwissAnalysis()

file_out = open("pyorbit_twiss_sizes_ekin_after_dtl_scan_analysis.dat","w")

s = " Node   position "
s += "   alphaX betaX emittX  normEmittX"
s += "   alphaY betaY emittY  normEmittY"
s += "   alphaZ betaZ emittZ  emittZphiMeV"
s += "   sizeX sizeY sizeZ_deg"
s += "   eKin Nparts "
file_out.write(s+"\n")
#file_out.write("pos sizeX sizeY sizeZdeg \n")
print " N node   position    sizeX  sizeY  sizeZdeg  eKin Nparts "

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
	gamma = bunch.getSyncParticle().gamma()
	beta = bunch.getSyncParticle().beta()
	twiss_analysis.analyzeBunch(bunch)
	x_rms = math.sqrt(twiss_analysis.getTwiss(0)[1]*twiss_analysis.getTwiss(0)[3])*1000.
	y_rms = math.sqrt(twiss_analysis.getTwiss(1)[1]*twiss_analysis.getTwiss(1)[3])*1000.
	z_rms = math.sqrt(twiss_analysis.getTwiss(2)[1]*twiss_analysis.getTwiss(2)[3])*1000.
	z_to_phase_coeff = bunch_gen.getZtoPhaseCoeff(bunch)
	z_rms_deg = z_to_phase_coeff*z_rms/1000.0
	nParts = bunch.getSizeGlobal()
	(alphaX,betaX,emittX) = (twiss_analysis.getTwiss(0)[0],twiss_analysis.getTwiss(0)[1],twiss_analysis.getTwiss(0)[3]*1.0e+6)
	(alphaY,betaY,emittY) = (twiss_analysis.getTwiss(1)[0],twiss_analysis.getTwiss(1)[1],twiss_analysis.getTwiss(1)[3]*1.0e+6)
	(alphaZ,betaZ,emittZ) = (twiss_analysis.getTwiss(2)[0],twiss_analysis.getTwiss(2)[1],twiss_analysis.getTwiss(2)[3]*1.0e+6)		 
	norm_emittX = emittX*gamma*beta
	norm_emittY = emittY*gamma*beta
	#---- phi_de_emittZ will be in [pi*deg*MeV]
	phi_de_emittZ = z_to_phase_coeff*emittZ	
	eKin = bunch.getSyncParticle().kinEnergy()*1.0e+3
	s = " %45s  %4.5f "%(node.getName(),pos+pos_start)
	s += "   %+6.4f  %+6.4f  %+6.4f  %+6.4f   "%(alphaX,betaX,emittX,norm_emittX)
	s += "   %+6.4f  %+6.4f  %+6.4f  %+6.4f   "%(alphaY,betaY,emittY,norm_emittY)
	s += "   %+6.4f  %+6.4f  %+6.4f  %+6.4f   "%(alphaZ,betaZ,emittZ,phi_de_emittZ)
	s += "   %5.3f  %5.3f  %5.3f "%(x_rms,y_rms,z_rms_deg)
	s += "  %10.6f   %8d "%(eKin,nParts)
	file_out.write(s +"\n")
	file_out.flush()
	#s_samll =  "%4.5f "%(pos+pos_start)+"  %5.3f  %5.3f   %5.3f "%(x_rms,y_rms,z_rms_deg)
	#file_out.write(s_samll +"\n")
	s_prt = " %5d  %35s  %4.5f "%(paramsDict["count"],node.getName(),pos+pos_start)
	s_prt += "  %5.3f  %5.3f   %5.3f "%(x_rms,y_rms,z_rms_deg)
	s_prt += "  %10.6f   %8d "%(eKin,nParts)
	print s_prt	
	
def action_exit(paramsDict):
	action_entrance(paramsDict)
	
actionContainer.addAction(action_entrance, AccActionsContainer.ENTRANCE)
actionContainer.addAction(action_exit, AccActionsContainer.EXIT)

#---- This is actual tracking of the bunch
time_start = time.clock()

accLattice.trackBunch(bunch, paramsDict = paramsDict, actionContainer = actionContainer)

time_exec = time.clock() - time_start
print "time[sec]=",time_exec

file_out.close()


