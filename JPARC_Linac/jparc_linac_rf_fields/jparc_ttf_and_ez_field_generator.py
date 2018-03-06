#-------------------------------------------------------------------------
# This script reads the SuperFish files with the field on the axis of RF gap
# and generates tables of the TTF functions and their polynomial interpolations. 
# TTFs (T,S,Tp,Sp) are funcftions of the cappa variable = 2*pi*f/(c*beta).
# The Super Fish files are the data files with electric and magnetic fields 
# along the RF cavities.
#
# For gaps in the RF Cavities the average field should be the same.
# It means the realtive amplitudes of the gaps will be inversly proportional 
# to the reference length of the gap.
#--------------------------------------------------------------------------


import sys
import math
import string

from orbit_utils import Function
from orbit_utils import SplineCH
from orbit_utils import GaussLegendreIntegrator
from orbit.utils.fitting import PolynomialFit
from orbit_utils import Polynomial

from orbit.py_linac.rf_field_readers import RF_AxisFieldAnalysis

# import the XmlDataAdaptor XML parser
from orbit.utils.xml import XmlDataAdaptor


def getAxisFieldFunction(fl_name):
	"""
	This function will read the .FSO file and extract the axis field Ez(z). 
	"""
	fl_in = open(fl_name,"r")
	lns = fl_in.readlines()
	fl_in.close()
	func = Function()
	i_start = -1
	for i in range(len(lns)):
		if(lns[i].find("    Z(cm)      Ez(V/m)") >= 0):
			i_start = i+1
			break
	for i in range(i_start,len(lns)):
		res_arr = lns[i].split()
		if(len(res_arr) != 2):
			break	
		x = 0.01*float(res_arr[0])
		y = float(res_arr[1])
		func.add(x,y)
	#fix for the field at the 0 region
	func1 = Function()
	for i in range(1,func.getSize()):
		ind = func.getSize() - i
		x = -func.x(ind)
		y = func.y(ind)
		func1.add(x,y)
	for i in range(func.getSize()):		
		x = func.x(i)
		y = func.y(i)
		func1.add(x,y)
	return func1
	
def DumpSplineIntoFile(fl_name,spline):
	"""
	Dumps spline with Ez(z) function into the file.
	"""
	fl_out = open(fl_name,"w")
	for iz in range(spline.getSize()):
		x = spline.x(iz)
		y = spline.y(iz)
		if(iz == 0): y = 0.
		if(iz == spline.getSize() -1): y = 0.
		st = " %15.9g  %16.10g "%(x,y)
		fl_out.write(st+"\n")
	fl_out.close()	

def SetUpETL_in_RFGAP_Da(rf_gap_analysis, rf_gap_da, rf_gap_beta_avg_dict, axis_file_name, beta_min, beta_max):
	"""
	Modifies RF gap XmlDataAdaptor to include actual polynomial coefficients for TTF functions and 
	E0L = E0TL/Tttf. It also dump Ez(z) spline into the external file.
	"""
	spline = rf_gap_analysis.getNormilizedSpline()
	DumpSplineIntoFile("./axis_fields/"+axis_file_name,spline)
	rf_gap_name = rf_gap_da.stringValue("name")
	parameters_da = rf_gap_da.childAdaptors("parameters")[0]
	TTFs_da = rf_gap_da.childAdaptors("TTFs")[0]
	polyT_da = TTFs_da.childAdaptors("polyT")[0]
	polyTP_da = TTFs_da.childAdaptors("polyTP")[0]
	polyS_da = TTFs_da.childAdaptors("polyS")[0]
	polySP_da = TTFs_da.childAdaptors("polySP")[0]
	#--------------------------------------------------------
	parameters_da.setValue("EzFile",axis_file_name)
	beta_avg = rf_gap_beta_avg_dict[rf_gap_name]
	(T,TP,S,SP) = rf_gap_analysis.getTTPandSSP_Values(beta_avg)
	E0TL = parameters_da.doubleValue("E0TL")
	E0L = E0TL/T
	parameters_da.setValue("E0L",E0L)
	#--------------------------------------------------------
	TTFs_da.setValue("beta_min",beta_min)
	TTFs_da.setValue("beta_max",beta_max)
	#--------------------------------------------------------
	gap_polynoms_coef_arr = rf_gap_analysis.gap_polynoms_coef_arr
	[t_coef_err_arr,tp_coef_err_arr,s_coef_err_arr,sp_coef_err_arr] = gap_polynoms_coef_arr[0]
	#--------------------------------------------------------
	[coef_arr,err_arr] = t_coef_err_arr
	order = len(coef_arr)-1
	st = ""
	for coeff in coef_arr:
		st += " %16.9e "%coeff
	polyT_da.setValue("order",order)
	polyT_da.setValue("pcoefs",st)
	#----------------------------------
	[coef_arr,err_arr] = tp_coef_err_arr
	order = len(coef_arr)-1
	st = ""
	for coeff in coef_arr:
		st += " %16.9e "%coeff
	polyTP_da.setValue("order",order)
	polyTP_da.setValue("pcoefs",st)
	"""
	#==== all rf gaps so far are symmetrical, so S and SP := 0
	#----------------------------------
	[coef_arr,err_arr] = s_coef_err_arr
	order = len(coef_arr)-1
	st = ""
	for coeff in coef_arr:
		st += " %16.9e "%coeff
	polyS_da.setValue("order",order)
	polyS_da.setValue("pcoefs",st)
	#----------------------------------
	[coef_arr,err_arr] = sp_coef_err_arr
	order = len(coef_arr)-1
	st = ""
	for coeff in coef_arr:
		st += " %16.9e "%coeff
	polySP_da.setValue("order",order)
	polySP_da.setValue("pcoefs",st)
	"""
	#----------------------------------	

def gapNameTransformDTL_SDTL(gap_name):
	tmp = "LI_" + gap_name
	if(tmp.find("D") >= 0):
		tmp = tmp.replace("D","DTL")
	tmp = tmp.replace("G",":RG")
	return tmp	

#--------------------------------------------------------
#            START SCRIPT
#--------------------------------------------------------

#==== read xml file with all data
xml_file_name = "./lattice_xml/jparc_linac_init.xml"
acc_da = XmlDataAdaptor.adaptorForFile(xml_file_name)
rf_gap_da_dict = {}
for seq_da in acc_da.childAdaptors():
	for accElem_da in seq_da.childAdaptors("accElement"):
		if(accElem_da.hasAttribute("type") and accElem_da.stringValue("type") == "RFGAP"):
			rf_gap_da_dict[accElem_da.stringValue("name")] = accElem_da

"""
for rf_gap_name in rf_gap_da_dict.keys():
	rf_gap_da = rf_gap_da_dict[rf_gap_name]
	print "debug rf_gap =",rf_gap_da.stringValue("name")
"""

#----- read beta average values for the RF gaps from the file
#----- prepared after pyORBIT tracking with BaseRF_Gap model.
#----- these data will be used to get E0L parameter from E0TL
#----- E0L = E0TL/T(beta_avg) where T - symmetric TTF
rf_gap_beta_avg_dict = {}
fl_in = open("./base_rf_tracking_data/pyorbit_jparc_rf_gaps_ekin_in_out.dat","r")
lns = fl_in.readlines()[1:]
fl_in.close()
for ln in lns:
	res_arr = ln.split()
	if(len(res_arr) == 4):
		gap_name = res_arr[0]
		beta_avg = float(res_arr[3])
		rf_gap_beta_avg_dict[gap_name] = beta_avg
		#print "debug gap=",gap_name," beta=",beta_avg



# rf_freq - rf frequency in Hz
rf_freq = 324.0e+6

data_root_dir = "./JPARC_SuperFish_FIELD_DATA/dtl_sdtl/"

dtl1_rf_gap_names = []
for i in range(1,77):
	name = "D1G"+string.zfill(str(i),2)
	dtl1_rf_gap_names.append(name)
	
dtl2_rf_gap_names = []
for i in range(1,44):
	name = "D2G"+string.zfill(str(i),2)
	dtl2_rf_gap_names.append(name)

dtl3_rf_gap_names = []
for i in range(1,28):
	name = "D3G"+string.zfill(str(i),2)
	dtl3_rf_gap_names.append(name)

sdtl_rf_cav_names = []
for i in range(1,17):
	nameA = "S"+string.zfill(str(i),2)+"A"
	nameB = "S"+string.zfill(str(i),2)+"B"
	sdtl_rf_cav_names.append(nameA)	
	sdtl_rf_cav_names.append(nameB)	
	

#data_in_arr[[cav_directory_name,cav_name,rf_gap_names_arr,rf_freq, beta_min, beta_max]]
data_in_arr = []
data_in_arr.append(["DTL1","DTL1",dtl1_rf_gap_names,rf_freq, 0.07, 0.22])
data_in_arr.append(["DTL2","DTL2",dtl2_rf_gap_names,rf_freq, 0.18, 0.30])
data_in_arr.append(["DTL3","DTL3",dtl3_rf_gap_names,rf_freq, 0.27, 0.35])	

for sdtl_cav in sdtl_rf_cav_names:
	sdtl_cav_rf_gap_names = []
	for i in range(1,6):
		sdtl_cav_rf_gap_names.append(sdtl_cav+"G"+string.zfill(str(i),2))
	data_in_arr.append(["SDTL",sdtl_cav,sdtl_cav_rf_gap_names,rf_freq, 0.27, 0.65])	
	
n_poly_order = 4
n_table_points = 100
dir_name = "ttf_results_data"

for [cav_directory_name,cav_name,rf_gap_names_arr,rf_freq, beta_min, beta_max] in data_in_arr:
	file_cav_name = cav_name+".dat"
	for rf_gap_name in rf_gap_names_arr:
		file_gap_name = data_root_dir + cav_directory_name + "/" + rf_gap_name + ".SFO"
		func = getAxisFieldFunction(file_gap_name)
		splineGap = SplineCH()					
		splineGap.compile(func)
		rf_gap_analysis = RF_AxisFieldAnalysis(splineGap, zeroIsCenter = True)
		rf_gap_analysis.makeTransitTimeTables(beta_min,beta_max,n_table_points,rf_freq)
		gap_polynoms_coef_arr = rf_gap_analysis.makePlynomialFittings(n_poly_order)
		gap_polynoms_t_tp_s_sp_err_arr = rf_gap_analysis.gap_polynoms_t_tp_s_sp_err_arr
		(err_t,err_tp,err_s,err_sp) = gap_polynoms_t_tp_s_sp_err_arr[0]
		print "gap=",rf_gap_name," spline size = ",splineGap.getSize()	," err = %5.6f %5.6f %5.6f %5.6f "%(err_t,err_tp,err_s,err_sp)
		rf_gap_analysis.dumpTransitTimeTables("./"+dir_name+"/"+rf_gap_name+"_t_tp_s_sp.dat")
		rf_gap_analysis.dumpTTFandFitting("./"+dir_name+"/"+rf_gap_name+"_ttf_fiitings.dat")
		#--------------------------------------------
		real_gap_name = gapNameTransformDTL_SDTL(rf_gap_name)
		rf_gap_da = rf_gap_da_dict[real_gap_name]
		axis_file_name = rf_gap_name+".dat"
		SetUpETL_in_RFGAP_Da(rf_gap_analysis, rf_gap_da, rf_gap_beta_avg_dict, axis_file_name, beta_min, beta_max)

#--------------------------------------------------------------------------
# ACS Part of the linac
# For ACS part the Ez(z) fields 
# are defined by the relativistic beta
#--------------------------------------------------------------------------
# rf_freq - rf frequency in Hz
rf_freq = 972.0e+6

beta_min = 0.45
beta_max = 0.75

n_poly_order = 4

dir_acs = "./JPARC_SuperFish_FIELD_DATA/acs/"

fl_in = open(dir_acs+"ACSgeobeta.txt")
lns_tmp = fl_in.readlines()
fl_in.close()

#-------------cav_name_beta_arr = ["cav_name","beta"]
#------------- cav names in XML PyORBIT file starts with like LI_ACS01A
cav_name_beta_arr = []

lns = []
for ln in lns_tmp:
	if(ln.find("ACS") >= 0):
		res_tmp_arr = ln.split()
		res_arr = []
		for st in res_tmp_arr:
			st = st.strip(",")
			res_arr.append(st)
		if(len(res_arr) == 2):
			cav_name_beta_arr.append(["LI_"+res_arr[0],"%6.4f"%(float(res_arr[1]))])
		else:
			cav_name_beta_arr.append(["LI_"+res_arr[0],"%6.4f"%(float(res_arr[2]))])
			cav_name_beta_arr.append(["LI_"+res_arr[1],"%6.4f"%(float(res_arr[2]))])

#---- cav_name_beta_dict[cav_name] = beta_st where beta_st it is a string
cav_name_beta_dict = {}
for [cav_name,beta_st] in cav_name_beta_arr:
	cav_name_beta_dict[cav_name] = beta_st

#---- axis_field_dict[beta_st] = RF_AxisFieldAnalysis instance
axis_field_dict = {}

for [cav_name,beta_st] in cav_name_beta_arr:
	#print "debug cav=",cav_name,"  beta=",beta_st
	if(not axis_field_dict.has_key(beta_st)):
		fl_axis_name = dir_acs+"data/"+beta_st+".SFO"
		func = getAxisFieldFunction(fl_axis_name)
		splineGap = SplineCH()					
		splineGap.compile(func)
		rf_gap_analysis = RF_AxisFieldAnalysis(splineGap, zeroIsCenter = True)
		rf_gap_analysis.makeTransitTimeTables(beta_min,beta_max,n_table_points,rf_freq)
		gap_polynoms_coef_arr = rf_gap_analysis.makePlynomialFittings(n_poly_order)	
		gap_polynoms_t_tp_s_sp_err_arr = rf_gap_analysis.gap_polynoms_t_tp_s_sp_err_arr
		(err_t,err_tp,err_s,err_sp) = gap_polynoms_t_tp_s_sp_err_arr[0]
		print "gap=",beta_st," spline size = ",splineGap.getSize()	," err = %5.6f %5.6f %5.6f %5.6f "%(err_t,err_tp,err_s,err_sp)
		axis_field_dict[beta_st] = rf_gap_analysis
		
count = 0
for rf_gap_name in rf_gap_da_dict.keys():
	rf_gap_da = rf_gap_da_dict[rf_gap_name]
	rf_gap_name = rf_gap_da.stringValue("name")
	cav_name = rf_gap_name.split(":")[0]
	if(cav_name.find("ACS") < 0): continue
	#print "debug rf_gap_name=",rf_gap_name
	beta_st = cav_name_beta_dict[cav_name]
	rf_gap_analysis = axis_field_dict[beta_st]
	axis_file_name = "ACS_"+beta_st+".dat"
	SetUpETL_in_RFGAP_Da(rf_gap_analysis, rf_gap_da, rf_gap_beta_avg_dict, axis_file_name, beta_min, beta_max)
	count += 1
	
print "Total ACS gaps =",count


#------ write changed XML data file with the linac lattice into the new file
xml_file_name = "./lattice_xml/jparc_linac_new.xml"
acc_da.writeToFile(xml_file_name)

sys.exit(0)
