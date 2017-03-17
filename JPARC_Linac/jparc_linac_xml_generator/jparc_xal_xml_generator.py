#! /usr/bin/env python

"""
This script reads the J-PARC XAL XML file with the lattice data.
"""

import sys
import time

# import the XmlDataAdaptor XML parser
from orbit.utils.xml import XmlDataAdaptor

from jparc_lattice_factory_lib import JPARC_Linac_Lattice_XAL_Generator


print "==============START======================="
#---- the XML file name with the structure
xml_file_name = "./jparc_xal_xml/jparc-LI_RCS-a400_5mA131214_S16.xdxf"
xml_file_name = "./jparc_xal_xml/jparc-LI_RCS-40mA_-2483_20160603.xdxf"
acc_da = XmlDataAdaptor.adaptorForFile(xml_file_name)
acc_seqs_init_da = acc_da.childAdaptors("sequence")
print "Acc seq n=",len(acc_seqs_init_da)

jparc_lattice_gen = JPARC_Linac_Lattice_XAL_Generator(acc_seqs_init_da)
lattice_da = jparc_lattice_gen.makeLattice_da()
lattice_da.writeToFile("../jparc_linac_lattice_xml/jparc_linac.xml")

print "Stop."
