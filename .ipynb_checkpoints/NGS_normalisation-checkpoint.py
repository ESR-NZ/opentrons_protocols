from opentrons import protocol_api 
import pandas as pd
import numpy as np
import string
import glob
import os
import time
import re

## Find excel file 
list_of_xlsx_files = glob.glob('Example_data/*.xlsx') # will need path of where these are on the robot file system
list_of_NXT_files = [worksheet for worksheet in list_of_xlsx_files if worksheet.startswith("Nextera")]
latest_file = max(list_of_xlsx_files, key=os.path.getctime)
## Time of file creation
c_time = os.path.getctime(latest_file)
local_time = time.ctime(c_time) 
    
## load up excel sheet
lib_pre_data = pd.read_excel(latest_file, skiprows=6, sheet_name=1)  

## Get the columns of interests from excel data
nmol4 = lib_pre_data.iloc[:,5]
RBS4 = lib_pre_data.iloc[:,6]

well_pos = lib_pre_data.iloc[:,7]
## Fix up the leading 0's
well_pos = well_pos.apply(lambda well: ''.join([x.lstrip("0") for x in re.split(r'([A-H])', well)]))


nmol2 = lib_pre_data.iloc[:,11]
RBS2 = lib_pre_data.iloc[:,12]

add_lib = list(lib_pre_data.iloc[:,8])

## Look at conc data
low_conc_mask = nmol4 > 14 ## this is a volume (uLs). Dont use more than 14ul.

nmol4_vols = list(nmol4[~low_conc_mask]) ## may need to fix up the leading 0's from the excel sheet (ie. 'A01 -> A1')
RBS4_vols = list(RBS4[~low_conc_mask])
wells_4 = list(well_pos[~low_conc_mask])

nmol2_vols = list(nmol2[low_conc_mask])
RBS2_vols = list(RBS2[low_conc_mask])
wells_2 = list(well_pos[low_conc_mask])

## Combine RSB into one list so we can do it in one run 
RBS_vols_combined = RBS4_vols + RBS2_vols
RBS_wells_combined = wells_4 + wells_2

## Robot setup ## 
metadata = {
    'apiLevel': '2.3',
    'author': 'Storey'}

def run(protocol: protocol_api.ProtocolContext):
    
    sample_plate = protocol.load_labware('opentrons_96_aluminumblock_nest_wellplate_100ul', 4)
    norm_plate = protocol.load_labware('opentrons_96_aluminumblock_nest_wellplate_100ul', 5) 
    RBS = protocol.load_labware('usascientific_12_reservoir_22ml', 11)['A1']
    pool = protocol.load_labware('opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap',6)
    
    
    tiprack_20ul_1 = protocol.load_labware('opentrons_96_filtertiprack_20ul', 8)
    tiprack_20ul_2 = protocol.load_labware('opentrons_96_filtertiprack_20ul', 9)
    
    p20 = protocol.load_instrument('p20_single_gen2', 'left', tip_racks = [tiprack_20ul_1,tiprack_20ul_2])
   
   
    ## Transfer RBS
    p20.transfer(RBS_vols_combined, 
                 RBS,
                 [norm_plate.wells_by_name()[well_name] for well_name in RBS_wells_combined])
    

    ## Transfer DNA    
    p20.transfer(nmol2_vols, 
                 [sample_plate.wells_by_name()[well_name] for well_name in wells_2],
                 [norm_plate.wells_by_name()[well_name] for well_name in wells_2], new_tip = 'always')
    
    p20.transfer(nmol4_vols, 
                 [sample_plate.wells_by_name()[well_name] for well_name in wells_4],
                 [norm_plate.wells_by_name()[well_name] for well_name in wells_4],new_tip = 'always')

    
    ## Pool sampels
    p20.transfer(add_lib, 
                 [sample_plate.wells_by_name()[well_name] for well_name in well_pos],
                 pool['A1'])
    