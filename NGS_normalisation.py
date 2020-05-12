from opentrons import protocol_api
import pandas as pd
import numpy as np
import string
import glob
import os
import time


list_of_xlsx_files = glob.glob('Example_data/*.xlsx') # will need path of where these are on the robot file system

list_of_NXT_files = [worksheet for worksheet in list_of_xlsx_files if worksheet.startswith("Nextera")]

latest_file = max(list_of_xlsx_files, key=os.path.getctime)

lib_pre_data = pd.read_excel(latest_file, usecols="B:E", skiprows=14, sheet_name=1)  

wells = lib_pre_data.








def run(protocol: protocol_api.ProtocolContext):

    c_time = os.path.getctime(latest_file)
    local_time = time.ctime(c_time) 
    print("Input file created:", local_time) 
    
    ## Continue if correct
    #input("Press Enter to continue...") ## This might be useful, dont know if it works on the robot. Could put it in the run() funtion as a pause

