The PyORBIT linac lattice file for J-PARC was generated from the XAL (J-PARC version)
xdxf file ./jparc_xal_xml/jparc-LI_RCS-40mA_-2483_20160603.xdxf

PyORBIT Script: jparc_xal_xml_generator.py
Library with classes: jparc_lattice_factory_lib.py
Output XML file: jparc_lattice.xml
Time: March 2017, J-PARC visit
Author: Andrei Shishlo

The XAL RF gap phases were redefined as 
phase = xal_phase + 180.
That is inheritance from the XAL initial approach to the RF gap model.


Only linac related sequences were translated from XAL. There were some 
abnormalities in the XAL xml file:

1. Some sequences have "sub-sequences" which are just rebunchers with one RF gap.
	The RF gaps from these rebunchers were extracted to the sequence level.
	
2. There are elements with a "GV" type. We considered them as "gate valve". They 
   will be MARKERs, and they will have 0 length in the PyORBIT XML.
   
3. There is one element with a "DV" type, which was treated as DCV with effLength = 0.4
   and the real length = 0.

    <node type="DV" id="LI_L3BT:PBVM01" pos="170.2624782" len="0.4">
      <attributes>
        <magnet polarity="+1" bendAngle="0.0" pathLength="0.4" len="0.4" dipoleEntrRotAngle="0.0" dipoleExitRotAngle="0.0" dfltMagFld="0.0" fieldFactor="1.0" fringeFlag="0" f1="0.0" magnetGap="0.0" fringeFactor1="0.45" fringeFactor2="2.8" fieldPathFlag="3" />
        <align x="0.0" y="0.0" z="0" pitch="0.0" yaw="0.0" roll="0.0" />
        <aperture shape="1" x="0.06" y="0.06" />
      </attributes>
      <ps main="RCS_IERM:PBVMPS01" />
      <channelsuite name="magnetsuite">
        <channel handle="fieldRB" signal="RCS_IERM:PBVMPS01:MON:B" settable="false" />
      </channelsuite>
    </node>
    
4. There are some elements with the position outside the sequence where they are suppose to be:
   debug bad node position! Seq.=LI_ACS19B L=2.80755 node=LI_ACS19B:FCT00 pos=7.5911432
   debug bad node position! Seq.=LI_ACS19B L=2.80755 node=LI_ACS19B:SCT00 pos=7.5911432
   debug bad node position! Seq.=LI_L3BD100 L=146.3360637 node=LI_L3BD100:BLMP66 pos=147.84634894
   debug bad node position! Seq.=LI_S16B L=2.5590817 node=LI_S16B:FCT00 pos=2.660617
   debug bad node position! Seq.=LI_S16B L=2.5590817 node=LI_S16B:SCT00 pos=2.660617
   
   We did following:
   moved LI_S16B:FCT00 and LI_S16B:SCT00 to LI_MEBT2
   removed LI_L3BD100:BLMP66 completely
   changed positions LI_ACS19B:FCT00 and LI_ACS19B:SCT00 to 1.9493572 instead of 7.5911432

5. LI_MEBT2 accelerator sequence RF gaps mode were wrong. The right modes for RF BNCH01 and BNCH02
   are 0,1,0,1,0 - 0,1,0,1,0. Each buncher has 10 RF gaps. The function LI_MEBT2_RF_Gaps_Mode_Fix was 
   introduced to fix this.