from opentrons import protocol_api
import pandas as pd
import numpy as np
import string
import sys
import glob
import os
import time
#!{sys.executable} -m pip install xlrd

# Set the required concentration of final diluted solution per well (ng/ul)
final_conc = 2

## Get most recent uploaded input file from plate reader 
list_of_xlsx_files = glob.glob('*.xlsx') # will need path of where these are on the robot file system
latest_file = max(list_of_xlsx_files, key=os.path.getctime)
c_time = os.path.getctime(latest_file)
local_time = time.ctime(c_time) 
print("Input file created:", local_time) 
## Continue if correct
#input("Press Enter to continue...") ## This might be useful

## Get raw data from excel file (sheet = "End point")
optima_raw = pd.read_excel(latest_file, usecols="B:E", skiprows=14)  

## Calc standard curve equation
stnds_values = optima_raw.loc[:4,4] ## Make sure this is set up clearly in the SOP 
stnds_concs = [0, 1000, 100, 10, 1] ## Make sure this is set up clearly in the SOP 

## Standard curve equation 
f = np.polyfit(stnds_values, stnds_concs, deg=1)

## Calc the concentrations of each sample.
sample_concs = (optima_raw.loc[:,:3]*f[0]+f[1])/10
## Set to Nan of too low to be useful. This will 'count' what samples are to be processed 
sample_concs = sample_concs.applymap(lambda x: np.NaN if x <= (final_conc * 2) else x)

## Complementry functions to calulate the dilution volumes
## Assume a max PCR vol avaliable of 20ul 

def get_pcr_prod(sample_conc, final_conc = final_conc):
    ## Amount of PCR product to get for nomalisation. Skip if too low.
    if sample_conc > 200:
        ## Just get the min the pipette can transfer
        vol1 = 1
    elif 10 <= sample_conc <= 200:  
        vol1 = (final_conc/sample_conc) * 100
    elif final_conc <= sample_conc < 10:
        ## Just take all of it 
        vol1 = 20
    else: #if the conc is too low or 
        return np.NaN 
    
    return round(vol1, 1)


def dilute(sample_conc, final_conc = final_conc):
    ## Amount of H2O to add to get the conc to final_conc
    if sample_conc > 200:
        ## We know we always have 1ul for all samples over 200ng/ul
        vol2 = (sample_conc/final_conc) - 1
    elif 10 <= sample_conc <= 200:  
        vol2 = 100 - ((final_conc/sample_conc) * 100)
        
    elif final_conc <= sample_conc < 10:
        vol2 = ((sample_conc * 20) / final_conc) - 20 
        
    else: #if the conc is too low 
        return np.NaN 
    
    return round(vol2, 1)
    
## Make some dfs for the amounts to be transfered
add_pcr_df = sample_concs.applymap(lambda conc: get_pcr_prod(conc))
add_water_df = sample_concs.applymap(lambda conc: dilute(conc))

## A 96 well plate coordinate dataframe
plate_96 = pd.DataFrame({k:[letter + str(k) for letter in string.ascii_uppercase[:8]] for k in range(1,13)})

## Dilutant volume and pos conversion from dataframes to lists for robot pipettes. 
# Where to use the p20 to add dilutant (transposed across plate to empty columns)
p20_dilute_pos_nans =  plate_96.shift(-4,axis='columns')[add_water_df < 20].T.values.flatten() ## convert df to list for robot input
p20_dilute_pos = [pos for pos in p20_dilute_pos_nans if type(pos) == str] ## drop nan vlaues from list

#Vol list for p20 addition of dilutant
p20_dilute_vol_nans = add_water_df[add_water_df < 20].T.values.flatten()
p20_dilute_vol = [vol for vol in p20_dilute_vol_nans if vol > 0]

#Where to use the p300 to add dilutant (transposed across plate to empty columns)
p300_dilute_pos_nans = plate_96.shift(-4,axis='columns')[add_water_df > 20].T.values.flatten()
p300_dilute_pos = [pos for pos in p300_dilute_pos_nans if type(pos) == str]

#Vol list for p300 addition of dilutant
p300_dilute_vol_nans = add_water_df[add_water_df > 20].T.values.flatten() 
p300_dilute_vol = [vol for vol in p300_dilute_vol_nans if vol > 0] 

## PCR vols and pos lists for pipettes (always going to be 20ul or less so p20)
##list of volumes of PCR product to transfer
p20_pcr_vols = [vol for vol in add_pcr_df.T.values.flatten() if vol > 0] ## dataframe to list with no NaNs

## List of positions to transfer some pcr product from
p20_pcr_pos_from_nan = plate_96[add_pcr_df > 0].T.values.flatten() 
p20_pcr_pos_from = [pos for pos in p20_pcr_pos_from_nan if type(pos) == str] ##rm NaNs and converts to list
## List of positions to transfer some pcr product
p20_pcr_pos_to_nan = plate_96.shift(-4,axis='columns')[add_pcr_df > 0].T.values.flatten()
p20_pcr_pos_to = [pos for pos in p20_pcr_pos_to_nan if type(pos) == str]

## Wells for doing end repair reaction 
EP_wells_nan = plate_96.shift(-8,axis='columns')[add_water_df > 0].T.values.flatten()
EP_wells = list(EP_wells_nan[pd.notna(EP_wells_nan)])

num_samples = len(EP_wells)

## Do end repair rxn with same pipettes 
## Calculate a master mix for the number of smaples + 5%.
"""
Per RNX
-------
1) DNA amplicons (2ng/ul) - 5ul
2) H20 - 7.5 ul 
3) Ultra II End Prep Reaction Buffer - 1.75
4) Ultra II End Prep Enzyme Mix - 0.75
"""

def add_6pc(ul, num_samples):
    initial_amount = (ul * num_samples)
    plus_6pc = initial_amount + (initial_amount * 0.06)
    if plus_6pc < 1:
        plus_6pc = 1
    return plus_6pc
    

master_mix = {'H20': add_6pc(7.5, num_samples),
             'ER_Buffer': add_6pc(1.75, num_samples),
             "ER_Enzyme": add_6pc(0.75, num_samples),
             "Lig_master_mix": add_6pc(17.5, num_samples),
             "Lig_enhance": add_6pc(0.5, num_samples)} 


### Barcoding logic ###
## set up array of the barcodes postion in plate/rack 
bc_pos_array = np.array([letter + str(num) for letter in string.ascii_uppercase[:4] for num in range(1,7)])
## select which wells to pick barcoeds from
bc_to_use = list(bc_pos_array[pd.notna(add_water_df).T.values.flatten()])
## Map to the wells to put the barcodes in
bc_to_wells = zip(bc_to_use, EP_wells)


## Robot setup ## 
metadata = {
    'apiLevel': '2.3',
    'author': 'Kemp and Storey'}

def run(protocol: protocol_api.ProtocolContext):
    # Create labware
    
    ## Sample plate on tempdeck.
    #tempdeck = protocol.load_module('tempdeck', 10)
    
    #sample_plate = tempdeck.load_labware('corning_96_wellplate_360ul_flat')
    sample_plate = protocol.load_labware('corning_96_wellplate_360ul_flat', 10)

    ## Reagents and solutions
    dilutant = protocol.load_labware('usascientific_12_reservoir_22ml', 11)['A1']
    enzyme_rack = protocol.load_labware('opentrons_24_tuberack_generic_2ml_screwcap', 7)
    ER_buffer = enzyme_rack['A1']
    ER_enzyme = enzyme_rack['A2']
    ER_mastermix_Tube = enzyme_rack['A3']
    Lig_master_mix = enzyme_rack['A4']
    Lig_enhance = enzyme_rack['A5']
    Lig_tube = enzyme_rack['A6']
    barcoded_combined_1 = enzyme_rack['D1'] ## Place to collect the barcoded frags together
    
    ## Barcodes rack 
    barcodes = protocol.load_labware('opentrons_24_tuberack_generic_2ml_screwcap', 4)
    
    ## Tips
    tiprack_20ul_1 = protocol.load_labware('opentrons_96_filtertiprack_20ul', 8)
    tiprack_20ul_2 = protocol.load_labware('opentrons_96_filtertiprack_20ul', 9)
    
    tiprack_300ul_1 = protocol.load_labware('opentrons_96_tiprack_300ul', 6)

    
    
    #Pipettes
    p20 = protocol.load_instrument('p20_single_gen2', 'left', tip_racks = [tiprack_20ul_1,tiprack_20ul_2])
    p300 = protocol.load_instrument('p300_single_gen2', 'right', tip_racks = [tiprack_300ul_1])
    
    ## Transfer any volumes of H2O under 20ul with the p20
    if p20_dilute_pos: ## chck for empty list, pass if empty to save time and tip 
        p20.distribute(p20_dilute_vol, dilutant, [sample_plate.wells_by_name()[well_name] for well_name in p20_dilute_pos])
    
    ## Transfer any volumes of H2O for dilution over 20ul with the p300
    if p300_dilute_vol:
        p300.distribute(p300_dilute_vol, dilutant, [sample_plate.wells_by_name()[well_name] for well_name in p300_dilute_pos])
    
       
    ## Transfer the PCR products over
    p20.transfer(p20_pcr_vols, 
                 [sample_plate.wells_by_name()[well_name] for well_name in p20_pcr_pos_from], 
                 [sample_plate.wells_by_name()[well_name] for well_name in p20_pcr_pos_to], 
                 new_tip='always'
    
    
    ## Make the end repair master mix in mastermix_Tube 
    # Transfer the H20
    if num_samples < 3:
         p20.transfer(master_mix['H20'], dilutant, ER_mastermix_Tube)
    else:
        p300.transfer(master_mix['H20'], dilutant, ER_mastermix_Tube)
    
    # buffer
    if num_samples < 11:
        p20.transfer(master_mix['ER_Buffer'], ER_buffer, ER_mastermix_Tube, mix_after=(3, 20)) 
        ## mix amount needs looking at 5/num_samples? 
        
    else:
        p300.transfer(master_mix['ER_Buffer'], ER_buffer, ER_mastermix_Tube, mix_after=(3, 8 * num_samples))
    
    p20.transfer(master_mix['ER_Enzyme'], ER_enzyme, ER_mastermix_Tube, mix_after=(3, 20))
    
    ## Dispense, mix before aspirate
    if num_samples < 3:
        p20.distribute(10, ER_mastermix_Tube, 
                        [sample_plate.wells_by_name()[well_name] for well_name in EP_wells], 
                        mix_before=(3, 8 * num_samples), 
                        disposal_volume=(num_samples*10)*0.04) ## set to 4%
    else:
        p300.distribute(10, ER_mastermix_Tube, 
                        [sample_plate.wells_by_name()[well_name] for well_name in EP_wells], 
                        mix_before=(3, 8 * num_samples),
                        disposal_volume=(num_samples*10)*0.04) ## set to 4%
        
    #tempdeck.set_temperature(25)
    
    ## Transfer the normalised PCR over to the wells with EP_mastermix
    p20.transfer(5, [sample_plate.wells_by_name()[well_name] for well_name in p20_pcr_pos_to],
        [sample_plate.wells_by_name()[well_name] for well_name in EP_wells], new_tip='always')
    
    
    #protocol.delay(minutes=20) ## incubate rxn
    #tempdeck.set_temperature(65) ## stop rxn
    
    
    ## Barcode mastermix 
    ## Add the ligation master mix
    if num_samples < 2:
         p20.transfer(master_mix['Lig_master_mix'], Lig_master_mix, Lig_tube)
    else:
        p300.transfer(master_mix['Lig_master_mix'], Lig_master_mix, Lig_tube)
        
    ## Add the ligation enhancer 
    if num_samples < 3: ## Add and mix with p20 
        p20.transfer(master_mix['Lig_enhance'], Lig_enhance, Lig_tube, mix_after=(3, 15))
    else: ## Mix afer addition with p300
        p20.transfer(master_mix['Lig_enhance'], Lig_enhance, Lig_tube)
        
        ## check total vol and ensure mix vol is within pipette range
        if master_mix['Lig_master_mix'] +  master_mix['Lig_enhance'] > 200:
            mix_vol = 200
        else:
            mix_vol = num_samples * 15
        
        p300.pick_up_tip()
        p300.mix(3, mix_vol, Lig_tube)
        p300.drop_tip()
    
    
    #protocol.delay(minutes=2)
    #tempdeck.set_temperature(12)
    
    
    ## Dispense the appropriate barcodes to the right wells
    p20.transfer(2.5, 
                 [barcodes.wells_by_name()[well_name] for well_name in bc_to_use],
                 [sample_plate.wells_by_name()[well_name] for well_name in EP_wells], new_tip='always')
    
    
    #tempdeck.set_temperature(22)
    
    ## Add the ligase 
    p20.transfer(18, Lig_tube,
                 [sample_plate.wells_by_name()[well_name] for well_name in EP_wells], new_tip='always', mix_after=(3, 20))
    
    
    ## Incubate the ligation
#     protocol.delay(minutes=30)
#     tempdeck.set_temperature(70)
#     protocol.delay(minutes=10)
#     tempdeck.set_temperature(10)
#     protocol.delay(seconds=30)
#     tempdeck.deactivate() 
    
    
    ## Collect all the samples together
    p300.consolidate(35.5, 
                     [sample_plate.wells_by_name()[well_name] for well_name in EP_wells],
                     enzyme_rack['D1'])
    
    
    
