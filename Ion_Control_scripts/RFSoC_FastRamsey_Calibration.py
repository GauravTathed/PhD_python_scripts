import time
import json, urllib.request
HOST = "pynq"
# HOST = "129.97.41.202"
PORT = 9009
URL  = f"http://{HOST}:{PORT}/upload_rows"
 
def upload_rows(rows, time_unit="us"):
    payload = {"rows": rows, "time_unit": time_unit}
    data = json.dumps(payload).encode()
    req = urllib.request.Request(URL, data=data, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read().decode())

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

import sys

import os
import datetime
import glob
import json
import numpy as np
import random
from scipy.optimize import curve_fit
sys.path.append(r'C:\Users\ions\Documents\IonControl\project_dir\QuditLab\config\Scripts')

from Functions_Data import *
from Functions_AWG import *
from Functions_RFSoC_RamseyCalibration import *
from Functions_Measurement import *
script_functions = (getGlobal, setGlobal, setScan, startScan, stopScan, getAllData, createTrace, closeTrace, plotPoint, scriptIsStopped, setEvaluation)

# get latest 545 and 623 freqs
f_offset = getGlobal('f_offset')
f_upper = getGlobal('f_upper')

check_3pt_coherence = False
ref_Ramsey_wait_time = 150
Ramsey_wait_time = 50
line_trigger = True
num_points = 31
line_trigger_waits = 1e6*np.arange(0/60/num_points,1/60+(1/60/num_points),(1/60/num_points)) # in us

# line_trigger_waits = np.concatenate((1e6*np.arange(15/60/num_points,1/60+(1/60/num_points),(1/60/num_points)),1e6*np.arange(0/60/num_points,1/60+(1/60/num_points),(1/60/num_points)),1e6*np.arange(0/60/num_points,1/60+(1/60/num_points),(1/60/num_points)),1e6*np.arange(0/60/num_points,1/60+(1/60/num_points),(1/60/num_points)),1e6*np.arange(0/60/num_points,1/60+(1/60/num_points),(1/60/num_points)),1e6*np.arange(0/60/num_points,1/60+(1/60/num_points),(1/60/num_points))))

#line_trigger_waits = [0,0,0,0,0,0,0,0,0]

# get latest pitimes for all delta-m transition types
pitime_n2 = getGlobal('pitime_n2') # [-2, 4, -4]
pitime_n1 = getGlobal('pitime_n1') # [-2, 3, -3]
pitime_0 = getGlobal('pitime_0') # [2, 4, 2]
pitime_p1 = getGlobal('pitime_p1') # [2, 4, 3]
pitime_p2 = getGlobal('pitime_p2') # [2, 4, 4]

ref_pitimes_input = [pitime_n2,pitime_n1,pitime_0,pitime_p1,pitime_p2]


f_offset_index = [0,2,0]
f_upper_index = [-1,4,-3]

for LTwait in line_trigger_waits:
            dt_string = datetime.datetime.now().strftime("%Y%m%d_%H%M")

            LTWaitTime = getGlobal('LineTriggerWait_Global')
            if not "%s"%LTWaitTime == "%s us"%0:
                setGlobal("LineTriggerWait_Global", 0, "us")
            
            f_offset_0 = run_fast_calibration(script_functions, 0, ref_Ramsey_wait_time,target_transition=f_offset_index,f_offset_input=None, f_upper_input=None,ref_pitimes=ref_pitimes_input,check_3pt_coherence=check_3pt_coherence)
            f_offset_dummy = getGlobal('f_offset')
            if not "%s"%f_offset_dummy == f_offset_0:
                setGlobal("f_offset", f_offset_0, "")

            f_upper_0 = run_fast_calibration(script_functions, 0, ref_Ramsey_wait_time,target_transition=f_upper_index,f_offset_input=None, f_upper_input=None,ref_pitimes=ref_pitimes_input,check_3pt_coherence=check_3pt_coherence)
            f_upper_dummy = getGlobal('f_upper')
            if not "%s"%f_upper_dummy == f_upper_0:
                setGlobal("f_upper", f_upper_0, "")


            LTWaitTime = getGlobal('LineTriggerWait_Global')
            if not "%s"%LTWaitTime == "%s us"%LTwait:
                setGlobal("LineTriggerWait_Global", LTwait, "us")
            
            f_offset_cal = run_fast_calibration(script_functions, 0, Ramsey_wait_time,target_transition=f_offset_index,f_offset_input=None, f_upper_input=None,ref_pitimes=ref_pitimes_input, write_to_IonControl = False)
            f_offset_dummy = getGlobal('f_offset_LT_varying')
            if not "%s"%f_offset_dummy == f_offset_cal:
                setGlobal("f_offset_LT_varying", f_offset_cal, "")

            f_upper_cal = run_fast_calibration(script_functions, 0, Ramsey_wait_time,target_transition=f_upper_index,f_offset_input=None, f_upper_input=None,ref_pitimes=ref_pitimes_input, write_to_IonControl = False)
            f_upper_dummy = getGlobal('f_upper_LT_varying')
            if not "%s"%f_upper_dummy == f_upper_cal:
                setGlobal("f_upper_LT_varying", f_upper_cal, "")
            
            output_file = fr'Z:\Lab Data\D52_Calibration_Ba137\Fast_Ramsey_calibrations_180Hz_PhaseSweep\properly_self_referenced_sweeps\ramsey_reference_transitions_LTWait{np.round(LTwait)}us_Shots250_RamseyWait{Ramsey_wait_time}us_{dt_string}.txt'
            with open(output_file,'w') as file:
                file.write(f"{LTwait}\n")
                file.write(f"[{f_offset_cal},{f_upper_cal}]\n")
                file.write(f"[{f_offset_index},{f_upper_index}]\n")
                file.write(f"[{f_offset_0},{f_upper_0}]\n")
