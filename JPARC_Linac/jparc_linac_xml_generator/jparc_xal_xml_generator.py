#! /usr/bin/env python

"""
This script reads the J-PARC XAL XML file with the lattice data.
"""

import sys
import time

# import the XmlDataAdaptor XML parser
from orbit.utils.xml import XmlDataAdaptor

from jparc_lattice_factory_lib import JPARC_Linac_Lattice_XAL_Generator
from jparc_lattice_factory_lib import JPARC_Linac_Lattice_Transformation
from jparc_lattice_factory_lib import LI_MEBT2_RF_Gaps_Mode_Fix

print "==============START======================="
#---- the XML file name with the structure
xml_file_name = "./jparc_xal_xml/jparc-LI_RCS-40mA_-2483_20160603.xdxf"
acc_da = XmlDataAdaptor.adaptorForFile(xml_file_name)
acc_seqs_init_da = acc_da.childAdaptors("sequence")
print "Acc seq n=",len(acc_seqs_init_da)

jparc_lattice_gen = JPARC_Linac_Lattice_XAL_Generator(acc_seqs_init_da)
lattice_da = jparc_lattice_gen.makeLattice_da()

#---- Transformation of the lattice: all RCS and ACS cavities 
#---- with A and B indeces will be combined
transformation = JPARC_Linac_Lattice_Transformation(lattice_da)
lattice_da = transformation.getTransformedLattice()

#---- Fix the mode parameters of RF gaps in the LI_MEBT2 sequence
LI_MEBT2_RF_Gaps_Mode_Fix(lattice_da)

lattice_da.writeToFile("../jparc_linac_lattice_xml/jparc_linac.xml")

print "Stop."
