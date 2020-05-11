from opentrons import protocol_api
import glob
import time
import sys

## Get most recent uploaded input file from plate reader 
list_of_xlsx_files = glob.glob('/root/plateReaderData/*.xlsx') # will need path of where these are on the robot file system

#latest_file = "/root/plateReaderData/pciogreen_pcr-20200414-xr.xlsx"
latest_file = max(list_of_xlsx_files, key=os.path.getctime)

c_time = os.path.getctime(latest_file)
local_time = time.ctime(c_time) 
print("Input file created:", local_time) 