#! /usr/bin/env python

"""
This script reads Trace3D file and prepares the quad-field dictionary.
"""

import sys
import time
import math

#---- phase wrappers from orbit.utils
from orbit.utils import phaseNearTargetPhase, phaseNearTargetPhaseDeg

class Trace3D_Record:
	def __init__(self,index,ln):
		self.index = index
		self.ln = ln
		self.name = ""
		self.typeId = -1
		self.params_arr = []
		self.parseLine()
		self.pos_start = 0.
		self.pos_end = 0.
		
	def parseLine(self):
		res_arr = self.ln.split("=")
		#print "debug res_arr =",res_arr
		self.name = res_arr[1].split(",")[0].strip("'")
		self.typeId = int(res_arr[2].split(",")[0])
		params_arr = res_arr[3].split(",")
		self.params_arr = []
		for par in params_arr:
			st = par.strip()
			if(len(st) > 0):
				self.params_arr.append(float(st))
		#print "debug name=",self.name
		#print "debug type=",self.typeId
		#print "debug pars=",self.params_arr
		#print "debug =============================="
		

def getLength(trace3d_record):
	tp = trace3d_record.typeId
	if(tp == 1): return trace3d_record.params_arr[0]
	if(tp == 3 or tp == 4): return trace3d_record.params_arr[1]
	if(tp == 8):
		alpha = trace3d_record.params_arr[0]*math.pi/180.
		rpho = trace3d_record.params_arr[1]
		return alpha*rpho
	return 0.

#------------------------------------
#          START OF THE SCRIPT
#------------------------------------


#fl_in = open("./jparc_trace3d_lattice/jparc-LI_RCS-50mA_DB2_0_20150108.t3d","r")
fl_in = open("./jparc_trace3d_lattice/jparc_matched2_40mA100.t3d","r")
lns_init = fl_in.readlines()
fl_in.close()

lns = []
for ln in lns_init:
	ln = ln.strip()
	if(ln.find("CMT") == 0):
		lns.append(ln)

print "debug n=",len(lns)

rec_arr = []
pos = 0.
for ind in range(len(lns)):
	rec = Trace3D_Record(ind,lns[ind])
	rec.pos_start = pos
	rec.pos_end = rec.pos_start + getLength(rec)
	pos = rec.pos_end
	rec_arr.append(rec)
	
quad_field_dict = {}
quad_length_dict = {}
quad_pos_start_stop_dict = {}
quad_names = []
types_arr = []
types_count_dict = {}
for rec in rec_arr:
	if(not rec.typeId in types_arr): 
		types_arr.append(rec.typeId)
		types_count_dict[rec.typeId] = 0
	types_count_dict[rec.typeId] += 1
	if(rec.typeId == 4 or rec.typeId == 3):
		if(not quad_field_dict.has_key(rec.name)):
			quad_field_dict[rec.name] = rec.params_arr[0]
			quad_length_dict[rec.name] = 0.
			quad_pos_start_stop_dict[rec.name] = [rec.pos_start,0.]
			quad_names.append(rec.name)
		quad_length_dict[rec.name] += rec.params_arr[1]
		quad_pos_start_stop_dict[rec.name][1] = rec.pos_end
		
print "debug N quads =",len(quad_names)
print "debug types =",types_arr
print "debug ==============================="
for typeId in types_arr:
	print "debug type=",typeId, " n elems.=",types_count_dict[typeId]

#=========Print Quads' fields[T/m], lengths[m], and positions[m]
#----- Trace3D has length in mm
fl_out = open("../jparc_linac_model/optics/jparc_matched2_40mA100_trace3d_quads.dat","w")
for quad_name in quad_names:
	field = quad_field_dict[quad_name]
	length = quad_length_dict[quad_name]/1000.
	pos = (quad_pos_start_stop_dict[quad_name][0]+quad_pos_start_stop_dict[quad_name][1])/2.
	pos = pos/1000.
	fl_out.write(quad_name + "  %+10.6f"%field+"  %10.6f"%length  + "  %12.6f"%pos  +"\n")
fl_out.close()

#=======Print RF Gaps E0TL in MeV, phases[deg]+180, position[m]
eKin = 3.0
fl_out = open("../jparc_linac_model/optics/jparc_matched2_40mA100_trace3d_rfgaps.dat","w")
for rec in rec_arr:
	if(rec.typeId == 10):
		rfgap_name = rec.name
		E0TL = rec.params_arr[0]
		phase = phaseNearTargetPhaseDeg(180.+rec.params_arr[1],0.)
		eKin_in = eKin
		eKin -= E0TL*math.cos(math.pi*phase/180.)
		eKin_out = eKin
		pos = rec.pos_start/1000.
		fl_out.write(rfgap_name + "  %11.9f"%E0TL+"  %+6.2f"%phase + " %10.6f"%pos + "     %8.6f"% eKin_in + "  %8.6f"% eKin_out  +"\n")
fl_out.close()		


print "total length =",rec_arr[len(rec_arr)-1].pos_end/1000


print "Stop."