#! /usr/bin/env python

"""
This is a library with classes and functions to generate the 
PyORBIT linac lattice by using Data Adaptor elements from J-PARC XAL
XML file.
"""

import sys
import time
import math

# import the XmlDataAdaptor XML parser
from orbit.utils.xml import XmlDataAdaptor

#---- phase wrappers from orbit.utils
from orbit.utils import phaseNearTargetPhase, phaseNearTargetPhaseDeg

class JPARC_Linac_Lattice_XAL_Generator:
	def __init__(self,xal_acc_seqs_init_da):
		self.xal_acc_seqs_da = self.filter_sequences(xal_acc_seqs_init_da)
		self.jparc_lattice_da = XmlDataAdaptor("JPARC_LINAC")
		self.acc_seqs_da = []
		self.acc_node_gen = Acc_Node_JPARC_Generator()

	def filter_sequences(self,xal_acc_seqs_init_da):
		"""
		filter out the original sequences (Data Adaptors)
		"""
		do_not_use_names = []
		do_not_use_names.append("RCS")
		do_not_use_names.append("LI_IS")
		do_not_use_names.append("LI_LEBT")
		do_not_use_names.append("LI_RFQ")		
		
		def filter_seq_name(seq_name):
			res = False
			for bad_name in do_not_use_names:
				res = res or (seq_name.find(bad_name) >= 0)
			return res
					
		xal_acc_seqs_da = []
		for xal_acc_seq_da in xal_acc_seqs_init_da:
			seq_name = xal_acc_seq_da.stringValue("id")
			#print "debug seq =",seq_name
			if(not filter_seq_name(seq_name)):
				xal_acc_seqs_da.append(xal_acc_seq_da)
		return xal_acc_seqs_da
		
	def makeLattice_da(self):
		self.jparc_lattice_da = XmlDataAdaptor("JPARC_LINAC")
		self.acc_seqs_da = []
		for xal_acc_seq_da in self.xal_acc_seqs_da:
			accSeq_da = self.makeSequence(xal_acc_seq_da)
			self.acc_seqs_da.append(accSeq_da)
			self.jparc_lattice_da.addChildAdaptor(accSeq_da)
		#---- perform the check of the elements' positions
		self.elementsPosCechk()
		return self.jparc_lattice_da
		
	def makeSequence(self,xal_acc_seq_da):
		accSeq_name = xal_acc_seq_da.stringValue("id")
		accSeq_da = XmlDataAdaptor(accSeq_name)
		accSeq_da.setValue("name",accSeq_name)
		accSeq_da.setValue("bpmFrequency",324.0E6)
		accSeq_length = xal_acc_seq_da.doubleValue("len")
		accSeq_da.setValue("length",accSeq_length)
		#=========================================
		cavities_da = []
		cav_name_gap_counter_dict = {}
		#=========================================
		nodes_da = []
		for xal_node_da in xal_acc_seq_da.childAdaptors("node"):
			node_da = self.acc_node_gen.getAccNode(xal_node_da)
			if(node_da != None): nodes_da.append(node_da)
		#=========================================
		#  let's calculate all RG not hidden in the sub-sequences
		xal_count_rg = 0
		xal_total_length = 0.
		for xal_node_da in xal_acc_seq_da.childAdaptors("node"):
			if(xal_node_da.hasAttribute("len")):
				if(xal_node_da.stringValue("type") != "DCV" and xal_node_da.stringValue("type") != "DCH"):
					xal_total_length += xal_node_da.doubleValue("len")
			if(xal_node_da.stringValue("type") == "RG"):
				xal_count_rg += 1
		#======= add RG that are in additional sub-sequences
		for xal_acc_seq_in_da in xal_acc_seq_da.childAdaptors("sequence"):
			node_da = self.acc_node_gen.getNodeRG_forSubSeq(cavities_da,cav_name_gap_counter_dict,xal_acc_seq_in_da)
			"""
			parent_name = xal_acc_seq_da.stringValue("id")
			child_name = xal_acc_seq_in_da.stringValue("id")
			child_type = xal_acc_seq_in_da.stringValue("type")
			child_len = xal_acc_seq_in_da.stringValue("len")
			child_nodes_da = xal_acc_seq_in_da.childAdaptors("node")
			child_nodes_n = len(child_nodes_da)
			st = "parent seq="+parent_name+" child seq="+child_name+" type="+child_type+" len="+child_len+"  n_childs="+str(child_nodes_n)
			print st	
			"""
			xal_count_rg += 1
			nodes_da.append(node_da)
		#===================================
		# add RG - RF Gaps and create additional cavities
		if(xal_acc_seq_da.hasAttribute("type")):
			if(xal_acc_seq_da.stringValue("type") == "DTLTank" or xal_acc_seq_da.stringValue("type") == "CCL"):
				for xal_node_da in xal_acc_seq_da.childAdaptors("node"):
					node_da = self.acc_node_gen.getNodeRG_forCavity(cavities_da,cav_name_gap_counter_dict,xal_acc_seq_da,xal_node_da)
					if(node_da != None): nodes_da.append(node_da)
		#===================================
		def positionComp(node1_da,node2_da):
			if(node1_da.doubleValue("pos") > node2_da.doubleValue("pos")):
				return 1
			else:
				if(node1_da.doubleValue("pos") == node2_da.doubleValue("pos")):
					return 0
			return -1
		#===================================
		nodes_da.sort(positionComp)
		orbit_count_rg = 0
		total_length = 0.
		for node_da in nodes_da:
			total_length += node_da.doubleValue("length")
			#print "debug node=",node_da.stringValue("name")," L=",node_da.doubleValue("length"), " total L=",total_length			
			if(node_da.stringValue("type") == "RFGAP"): orbit_count_rg += 1
			accSeq_da.addChildAdaptor(node_da)
		st = "debug difference seq= %10s "%accSeq_name
		st += " xal_rg_n = %2d"%xal_count_rg
		st += " orbit_rfgap_n = %2d"%orbit_count_rg
		st += " l_diff= %+8.6f"%(xal_total_length-total_length)
		st += " L = %8.4f"%accSeq_da.doubleValue("length")
		print st
		#   set up the cavities positions as average of all rf gaps
		#   and modes for these gaps
		self.acc_node_gen.setCavitiesPositions(cavities_da,nodes_da)
		#===================================
		cavities_in_seq_da = accSeq_da.createChild("Cavities")
		for cavity_da in cavities_da:
			cavities_in_seq_da.addChildAdaptor(cavity_da)
		#===================================
		return accSeq_da
		
	def elementsPosCechk(self):
		"""
		It will check that the positions of all accElements are inside the sequence.
		It also will move some elements to the right sequence
		"""
		eps = 0.0000001
		#---------move LI_S16B:FCT00 and LI_S16B:SCT00 to LI_MEBT2
		seq0 = None
		seq1 = None
		for acc_seq_da in self.acc_seqs_da:
			if(acc_seq_da.stringValue("name") == "LI_S16B"):
				seq0 = acc_seq_da
			if(acc_seq_da.stringValue("name") == "LI_MEBT2"):
				seq1 = acc_seq_da
		bad_nodes_da = []
		for node_da in seq0.childAdaptors("accElement"):
			if(node_da.stringValue("name") == "LI_S16B:FCT00"):
				bad_nodes_da.append(node_da)
			if(node_da.stringValue("name") == "LI_S16B:SCT00"):
				bad_nodes_da.append(node_da)
		length = seq0.doubleValue("length")
		nodes_da = seq0.childAdaptors()
		for node_da in bad_nodes_da:
			nodes_da.remove(node_da)
			pos0 = node_da.doubleValue("pos")
			pos1 = pos0 - length
			node_da.setValue("pos",pos1)
			seq1.addChildAdaptor(node_da)
		nodes_da = seq1.childAdaptors()
		#---------------------------------------------------------
		#---- remove node=LI_L3BD100:BLMP66 from LI_L3BD100
		for acc_seq_da in self.acc_seqs_da:
			if(acc_seq_da.stringValue("name") == "LI_L3BD100"):
				nodes_da = acc_seq_da.childAdaptors()
				bad_node_da = None
				for node_da in nodes_da:
					if(node_da.hasAttribute("name") and node_da.stringValue("name") == "LI_L3BD100:BLMP66"):
						bad_node_da = node_da
				nodes_da.remove(bad_node_da)
		#----------------------------------------------------------
		#---- set pos to 1.9493572 for LI_ACS19B:FCT00 and LI_ACS19B:SCT00 
		for acc_seq_da in self.acc_seqs_da:
			if(acc_seq_da.stringValue("name") == "LI_ACS19B"):
				nodes_da = acc_seq_da.childAdaptors("accElement")
				for node_da in nodes_da:
					if(node_da.stringValue("name") == "LI_ACS19B:FCT00" or node_da.stringValue("name") == "LI_ACS19B:SCT00"):
						node_da.setValue("pos",1.9493572)
		#---------------------------------------------------------
		for acc_seq_da in self.acc_seqs_da:
			length = acc_seq_da.doubleValue("length")
			for node_da in acc_seq_da.childAdaptors("accElement"):
				pos = node_da.doubleValue("pos")
				if(pos > length + eps):
					st = "debug bad node position! Seq.="+acc_seq_da.stringValue("name")
					st += " L="+str(length)
					st += " node="+node_da.stringValue("name")
					st += " pos="+str(pos)
					print st		
		
		
class Acc_Node_JPARC_Generator:
	def __init__(self):
		pass
	
	def getAccNode(self,xal_node_da):
		"""
		This is node_da generator for the following types
		DCV     DCH
		PMQV    PMQH
		"""
		node_da = XmlDataAdaptor("accElement")
		node_name = xal_node_da.stringValue("id")
		node_da.setValue("name",node_name)
		pos = xal_node_da.doubleValue("pos")
		node_da.setValue("pos",pos)			
		#======================================
		type_name = xal_node_da.stringValue("type")
		#======================================
		# for type DV see the read me file
		if(type_name == "DCV" or type_name == "DCH" or type_name == "DV"):
			if(type_name == "DV"): type_name == "DCV"
			node_da.setValue("type",type_name)
			node_da.setValue("length",0.)
			params_da = node_da.createChild("parameters")
			params_da.setValue("B",0.)
			params_da.setValue("effLength",xal_node_da.doubleValue("len"))
			return node_da
		#======================================
		if(type_name== "PMQV" or type_name== "PMQH"):
			node_da.setValue("length",xal_node_da.doubleValue("len"))
			node_da.setValue("type","QUAD")
			xal_attr_da = xal_node_da.childAdaptors("attributes")[0]
			xal_magnet_da = xal_attr_da.childAdaptors("magnet")[0]
			xal_aperture_da = xal_attr_da.childAdaptors("aperture")[0]
			params_da = node_da.createChild("parameters")
			params_da.setValue("field",xal_magnet_da.doubleValue("dfltMagFld"))
			params_da.setValue("aperture",2*xal_aperture_da.doubleValue("x"))
			params_da.setValue("aprt_type",1)
			params_da.setValue("radIn",xal_magnet_da.doubleValue("radIn"))
			params_da.setValue("radOut",xal_magnet_da.doubleValue("radOut"))
			return node_da
		#======================================
		node_type_marker = (type_name == "Marker")
		node_type_marker = node_type_marker or (type_name == "BLM")
		node_type_marker = node_type_marker or (type_name == "FCT")
		node_type_marker = node_type_marker or (type_name == "BCM")
		node_type_marker = node_type_marker or (type_name == "WS")
		node_type_marker = node_type_marker or (type_name == "GV")
		if(node_type_marker):
			node_da.setValue("type","MARKER")
			node_da.setValue("length",0.)
			params_da = node_da.createChild("parameters")
			return node_da
		#======================================
		if(type_name== "DH"):
			node_da.setValue("length",xal_node_da.doubleValue("len"))
			node_da.setValue("type","BEND")
			xal_attr_da = xal_node_da.childAdaptors("attributes")[0]
			xal_magnet_da = xal_attr_da.childAdaptors("magnet")[0]
			bendAngle = xal_magnet_da.doubleValue("bendAngle")*math.pi/180.
			path_length = xal_magnet_da.doubleValue("pathLength")
			xal_aperture_da = xal_attr_da.childAdaptors("aperture")[0]
			aptr_x = 2*xal_aperture_da.doubleValue("x")
			aptr_y = 2*xal_aperture_da.doubleValue("y")
			#node_da.setValue("length",path_length)
			params_da = node_da.createChild("parameters")
			params_da.setValue("aprt_type",3)
			params_da.setValue("aperture_x",aptr_x)
			params_da.setValue("aperture_y",aptr_y)
			params_da.setValue("ea1",0.)
			params_da.setValue("ea2",0.)
			params_da.setValue("kls","")
			params_da.setValue("poles","")
			params_da.setValue("skews","")
			params_da.setValue("theta",bendAngle)
			return node_da
		return None
			
	def getNodeRG_forSubSeq(self,cavities_da,cav_name_gap_counter_dict,xal_acc_seq_in_da):
		"""
		It is analysis for a subsequence which is Buncher with one RG
		It returns RF gap data adaptor and update cavities_da and cav_name_gap_counter_dict
		"""
		pos = xal_acc_seq_in_da.doubleValue("pos")
		attrb_da = xal_acc_seq_in_da.childAdaptors("attributes")[0]
		xal_rfcavity_da = attrb_da.childAdaptors("rfcavity")[0]
		E0TL = xal_rfcavity_da.doubleValue("amp")/1.0e+9
		phase = xal_rfcavity_da.doubleValue("phase")
		freq = xal_rfcavity_da.doubleValue("freq")*1.0E+6
		rg_node_da = xal_acc_seq_in_da.childAdaptors("node")[0]
		node_name = rg_node_da.stringValue("id")
		aperture_da = attrb_da.childAdaptors("aperture")[0]
		aperture = 2*aperture_da.doubleValue("x")
		st_arr = node_name.split(":")
		rf_cavity_name = st_arr[0]+":"+st_arr[1].split("_")[0]
		#print "rf cav name=",rf_cavity_name," E0TL=",E0TL*1.0E+3," freq=",freq/1.0E+6
		gap_count = -1
		if(cav_name_gap_counter_dict.has_key(rf_cavity_name)):
			cav_name_gap_counter_dict[rf_cavity_name] += 1
		else:
			cav_name_gap_counter_dict[rf_cavity_name] = 1
			#========================================================
			cavity_da = XmlDataAdaptor("Cavity")
			cavity_da.setValue("ampl",1.0)
			cavity_da.setValue("frequency","%12.5e"%freq)
			cavity_da.setValue("name",rf_cavity_name)
			cavity_da.setValue("pos",pos)
			cavities_da.append(cavity_da)
			#========================================================
		gap_count = cav_name_gap_counter_dict[rf_cavity_name]			
		st_count = str(gap_count)
		if(len(st_count) != 2): st_count = "0"+st_count
		node_da = XmlDataAdaptor("accElement")
		node_da.setValue("name",rf_cavity_name+":Rg"+st_count)
		node_da.setValue("pos",pos)
		node_da.setValue("type","RFGAP")
		node_da.setValue("length",0.0)
		parameters_da = node_da.createChild("parameters")
		parameters_da.setValue("E0L",E0TL)
		parameters_da.setValue("E0TL",E0TL)
		parameters_da.setValue("EzFile","")
		parameters_da.setValue("aperture",aperture)
		parameters_da.setValue("aprt_type",1)
		parameters_da.setValue("cavity",rf_cavity_name)
		parameters_da.setValue("mode",0)
		parameters_da.setValue("phase",phase)
		ttfs_da = node_da.createChild("TTFs")
		ttfs_da.setValue("beta_max",1.0)
		ttfs_da.setValue("beta_min",0.0)
		poly_t_da = ttfs_da.createChild("polyT")
		poly_s_da = ttfs_da.createChild("polyS")
		poly_tp_da = ttfs_da.createChild("polyTP")
		poly_sp_da = ttfs_da.createChild("polySP")
		poly_t_da.setValue("order",0)
		poly_t_da.setValue("pcoefs","1.0 ")
		poly_tp_da.setValue("order",0)
		poly_tp_da.setValue("pcoefs","0 ")
		poly_s_da.setValue("order",0)
		poly_s_da.setValue("pcoefs","0 ")
		poly_sp_da.setValue("order",0)
		poly_sp_da.setValue("pcoefs","0 ")
		return node_da

	def getNodeRG_forCavity(self,cavities_da,cav_name_gap_counter_dict,xal_acc_seq_da,xal_node_da):	
		"""
		It returns RF gap data adaptor and update cavities_da and cav_name_gap_counter_dict if there is a
		new cavity.
		"""		
		type_name = xal_node_da.stringValue("type")
		if(type_name != "RG"): return None
		rf_cavity_name = xal_acc_seq_da.stringValue("id")
		attrb_da = xal_acc_seq_da.childAdaptors("attributes")[0]
		xal_rfcavity_da = attrb_da.childAdaptors("rfcavity")[0]
		E0TL = xal_rfcavity_da.doubleValue("amp")/1.0e+9
		phase = xal_rfcavity_da.doubleValue("phase")
		freq = xal_rfcavity_da.doubleValue("freq")*1.0E+6
		mode = xal_rfcavity_da.stringValue("structureMode")
		if(cav_name_gap_counter_dict.has_key(rf_cavity_name)):
			cav_name_gap_counter_dict[rf_cavity_name] += 1
		else:
			cav_name_gap_counter_dict[rf_cavity_name] = 1
			cavity_da = XmlDataAdaptor("Cavity")
			cavity_da.setValue("ampl",1.0)
			cavity_da.setValue("frequency","%12.5e"%freq)
			cavity_da.setValue("name",rf_cavity_name)
			cavity_da.setValue("pos",0.)
			cavities_da.append(cavity_da)
		#====== make node from xal node
		node_name = xal_node_da.stringValue("id")
		pos = xal_node_da.doubleValue("pos")		
		xal_attrb_da = xal_node_da.childAdaptors("attributes")[0]
		xal_rf_gap_da = xal_attrb_da.childAdaptors("rfgap")[0]
		amp_factor = xal_rf_gap_da.doubleValue("ampFactor")
		phase_factor = xal_rf_gap_da.doubleValue("phaseFactor")
		E0TL = E0TL*amp_factor
		phase += phase_factor	
		xal_aperture_da = xal_attrb_da.childAdaptors("aperture")[0]
		aperture = 2*xal_aperture_da.doubleValue("x")
		node_da = XmlDataAdaptor("accElement")
		node_da.setValue("name",node_name)
		node_da.setValue("pos",pos)
		node_da.setValue("type","RFGAP")
		node_da.setValue("length",0.0)
		parameters_da = node_da.createChild("parameters")
		parameters_da.setValue("E0L",E0TL)
		parameters_da.setValue("E0TL",E0TL)
		parameters_da.setValue("EzFile","")
		parameters_da.setValue("aperture",aperture)
		parameters_da.setValue("aprt_type",1)
		parameters_da.setValue("cavity",rf_cavity_name)
		parameters_da.setValue("mode",mode)
		parameters_da.setValue("phase",phase)
		ttfs_da = node_da.createChild("TTFs")
		ttfs_da.setValue("beta_max",1.0)
		ttfs_da.setValue("beta_min",0.0)
		poly_t_da = ttfs_da.createChild("polyT")
		poly_s_da = ttfs_da.createChild("polyS")
		poly_tp_da = ttfs_da.createChild("polyTP")
		poly_sp_da = ttfs_da.createChild("polySP")
		poly_t_da.setValue("order",0)
		poly_t_da.setValue("pcoefs","1.0 ")
		poly_tp_da.setValue("order",0)
		poly_tp_da.setValue("pcoefs","0 ")
		poly_s_da.setValue("order",0)
		poly_s_da.setValue("pcoefs","0 ")
		poly_sp_da.setValue("order",0)
		poly_sp_da.setValue("pcoefs","0 ")
		return node_da

	def setCavitiesPositions(self,cavities_da,nodes_da):
		#cav_rfgaps_dict[cav_name] = [node_gap_da0,...]
		cav_rfgaps_dict = {}
		for node_da in nodes_da:
			if(node_da.stringValue("type") == "RFGAP"):
				params_da = node_da.childAdaptors("parameters")[0]
				cav_name = params_da.stringValue("cavity")
				if(cav_rfgaps_dict.has_key(cav_name)):
					cav_rfgaps_dict[cav_name].append(node_da)
				else:
					cav_rfgaps_dict[cav_name] = [node_da,]
		#---- calculate cavity position as a avg position of the gaps
		cav_avg_pos_dict = {}
		for cav_name in cav_rfgaps_dict.keys():
			avg_pos = 0.
			for rfgap_da in cav_rfgaps_dict[cav_name]:
				avg_pos += rfgap_da.doubleValue("pos")
			avg_pos /= len(cav_rfgaps_dict[cav_name])
			cav_avg_pos_dict[cav_name] = avg_pos
		for cavity_da in cavities_da:
			cav_name = cavity_da.stringValue("name")
			pos = cav_avg_pos_dict[cav_name]
			cavity_da.setValue("pos",pos)
			#------- set modes for cavities with mode != 0
			rfgap_da = cav_rfgaps_dict[cav_name][0]
			params_da = rfgap_da.childAdaptors("parameters")[0]
			mode0 = params_da.intValue("mode")
			if(mode0 > 0):
				mode = 0
				for rfgap_da in cav_rfgaps_dict[cav_name]:
					params_da = rfgap_da.childAdaptors("parameters")[0]
					params_da.setValue("mode",mode)
					if(mode == 0):
						mode = 1
					else:
						mode = 0
		#---- redefine cavities' phases to the PyORBIT definitions
		for cavity_da in cavities_da:
			cav_name = cavity_da.stringValue("name")
			for rfgap_da in cav_rfgaps_dict[cav_name]:
				params_da = rfgap_da.childAdaptors("parameters")[0]
				phase = params_da.doubleValue("phase")
				phase = phaseNearTargetPhaseDeg(180.+phase,0.)
				params_da.setValue("phase",phase)

class JPARC_Linac_Lattice_Transformation:
	"""
	This class will reorder the acc. sequences according to the 
	specific order defined by user, not by XAL xml file.
	It also will combine the A and B parts of the SDTL and RCS 
	cavities.
	"""
	def __init__(self,jparc_lattice_da):
		self.jparc_lattice_da = jparc_lattice_da
		self.acc_seqs_init_da = jparc_lattice_da.childAdaptors()
		#----------------------------------------------------------------
		self.acc_seq_names = ["LI_MEBT1","LI_DTL1","LI_DTL2","LI_DTL3"]
		for ind in range(16):
			seq_name = "LI_S"+"%02d"%(ind+1)
			self.acc_seq_names.append(seq_name)
		self.acc_seq_names.append("LI_MEBT2")
		for ind in range(21):
			seq_name = "LI_ACS"+"%02d"%(ind+1)
			self.acc_seq_names.append(seq_name)
		self.acc_seq_names.append("LI_L3BT")
		self.acc_seq_names.append("LI_L3BD0")
		self.acc_seq_names.append("LI_L3BD100")
		self.acc_seq_names.append("LI_L3BD30")
		self.acc_seq_names.append("LI_L3BD90")
		self.acc_seq_names.append("LI_L3L3ANA")
		self.acc_seq_names.append("LI_BD0")
		self.acc_seq_names.append("LI_BD30")
		self.acc_seq_names.append("LI_BD100")
		self.acc_seq_names.append("LI_BD90")
		self.acc_seq_names.append("LI_L3ANA")
		#----------------------------------------------------------------
			
	def getTransformedLattice(self):
		jparc_lattice_da = XmlDataAdaptor("JPARC_LINAC")
		for seq_name in self.acc_seq_names:
			accSeq_da = self.getAccSeq(seq_name)
			jparc_lattice_da.addChildAdaptor(accSeq_da)
		return jparc_lattice_da
			
	def getAccSeq(self,seq_name):
		seq_da = None
		if(seq_name.find("LI_S") < 0 and seq_name.find("LI_ACS") < 0):
			for acc_seq_da in self.acc_seqs_init_da:
				if(acc_seq_da.stringValue("name") == seq_name):
					seq_da = acc_seq_da
		else:
			acc_seqs_ab_da = []
			for st in ["A","B"]:
				seq_name_tmp = seq_name+st
				for acc_seq_da in self.acc_seqs_init_da:
					if(acc_seq_da.stringValue("name") == seq_name_tmp):
						acc_seqs_ab_da.append(acc_seq_da)
			seq_da = self.combineTwoSeqs(seq_name,acc_seqs_ab_da[0],acc_seqs_ab_da[1])
		return seq_da
		
	def combineTwoSeqs(self,seq_name,acc_seq0_da,acc_seq1_da):
		accSeq_da = XmlDataAdaptor(seq_name)
		accSeq_da.setValue("name",seq_name)
		accSeq_da.setValue("bpmFrequency",324.0E6)
		length0 = acc_seq0_da.doubleValue("length")
		length1 = acc_seq1_da.doubleValue("length")
		accSeq_da.setValue("length",length0+length1)
		accElems0_da = acc_seq0_da.childAdaptors("accElement")
		accElems1_da = acc_seq1_da.childAdaptors("accElement")
		for accElem_da in accElems1_da:
			pos = accElem_da.doubleValue("pos") + length0
			accElem_da.setValue("pos",pos)
		accElems_da = accElems0_da + accElems1_da
		#-----------------------------------------
		def positionComp(node1_da,node2_da):
			if(node1_da.doubleValue("pos") > node2_da.doubleValue("pos")):
				return 1
			else:
				if(node1_da.doubleValue("pos") == node2_da.doubleValue("pos")):
					return 0
			return -1
		#===================================
		accElems_da.sort(positionComp)		
		#-----------------------------------------
		cav_pos_avg = 0.
		count = 0
		for accElem_da in accElems_da:
			if(accElem_da.stringValue("type") == "RFGAP"):
				cav_pos_avg += accElem_da.doubleValue("pos")
				params_da = accElem_da.childAdaptors("parameters")[0]
				params_da.setValue("cavity",seq_name)
				count += 1
		cav_pos_avg /= count
		cavs_da = acc_seq0_da.childAdaptors("Cavities")[0]
		cav_da = cavs_da.childAdaptors("Cavity")[0]
		cav_da.setValue("pos",cav_pos_avg)
		cav_da.setValue("name",seq_name)
		for accElem_da in accElems_da:
			accSeq_da.addChildAdaptor(accElem_da)
		cavs_da = accSeq_da.createChild("Cavities")
		cavs_da.addChildAdaptor(cav_da)
		return accSeq_da
			
			
def LI_MEBT2_RF_Gaps_Mode_Fix(jparc_lattice_da):
	"""
	This function will assign the mode parameters to the RF gaps in
	the BNCH01 and BNCH02 cavities.Each cavity has 5 and 5 RF gaps (total 10).
	Each group of 5 gaps has modes 0,1,0,1,0. The XAL XML file does not have
	this information.
	"""
	count = 0
	for accSeq_da in jparc_lattice_da.childAdaptors():
		if(accSeq_da.stringValue("name") == "LI_MEBT2"):
			for accElem_da in accSeq_da.childAdaptors("accElement"):
				if(accElem_da.stringValue("type") == "RFGAP"):
					params_da = accElem_da.childAdaptors("parameters")[0]
					ind = count % 5
					mode = 0
					if(ind == 1 or ind == 3): mode = 1
					params_da.setValue("mode",mode)
					count += 1
					
			


