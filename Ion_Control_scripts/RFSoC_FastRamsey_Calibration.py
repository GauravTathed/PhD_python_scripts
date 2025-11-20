#D52_Qudit_heralded_Ramsey.py created 2024-09-26 15:04:20.886427

import json
import urllib.request
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
 
import sys
import os
import datetime
import glob
import json
import numpy as np
import shutil
from scipy.optimize import curve_fit
sys.path.append(r'C:\Users\ions\Documents\IonControl\project_dir\QuditLab\config\Scripts')
from Functions_Data import *
from Functions_AWG import *
from Functions_RFSoC_RamseyCalibration import *
from Functions_Measurement import *
script_functions = (getGlobal, setGlobal, setScan, startScan, stopScan, getAllData, createTrace, closeTrace, plotPoint, scriptIsStopped, setEvaluation)
dt_string = datetime.datetime.now().strftime("%Y%m%d_%H%M")
year = datetime.datetime.now().strftime("%Y")
month = datetime.datetime.now().strftime("%m")
day = datetime.datetime.now().strftime("%d")
pattern = "Z:\Lab Data\Qudit_Ramsey_raw_data\Raw_data\qudit_ramsey_scan_*"
def insert(x,p):
    return x[:int((len(x))/2)] + [p] +x[int((len(x))/2):]
Side_band_cooling_reps = 0


#initial_state = [[0,2,2]]

do_calibrations_0 = False
Measure_U1_only = False
full_phase_scan = True

#ramsey_wait_times = np.arange(0,5100,500)
repeat_experiments = 1
ramsey_wait_times = [1]
f_offset = float(getGlobal('f_offset'))
f_upper = float(getGlobal('f_upper'))

pitime_n2 = float(getGlobal('pitime_n2')) # [-2, 4, -4]
pitime_n1 = float(getGlobal('pitime_n1')) # [-2, 3, -3]
pitime_0 = float(getGlobal('pitime_0')) # [2, 4, 2]
pitime_p1 = float(getGlobal('pitime_p1')) # [2, 4, 3]
pitime_p2 = float(getGlobal('pitime_p2')) # [2, 4, 4]


#waitTime_list  = [1,10,300,400,500,600,700,800,900,100,350,450,550,650,750,850,950,200]
#waitTime_list = [1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500]
waitTime_list = [0]
wait_trans_piTime = list(Get_1762_PiTimes([[0,0,0]],pitime_n2,pitime_n1,pitime_0,pitime_p1,pitime_p2))[0]
print(wait_trans_piTime)
for wtr in waitTime_list:
    num_repeat = 0

    
    #d=2 
    #, (np.sin((400/1013.667)*(np.pi/2)))**2, , [-2,4,-3]
    initial_state = [[0, 3, 2]]

    pulse_train_U1 = [[0, 3, 2], [0, 2, 2], [0,0,0]]
    fractions_U1 = [0.5, 1,(np.sin((wtr/wait_trans_piTime)*(np.pi/2)))**2]
    simulated_phase_mask_U1 = [0, 0, 0]
    fixed_phase_mask_U1 = [1, 0, 0]

    pulse_train_U2 = [[0, 2, 2], [0, 3, 2]]
    fractions_U2 = [1, 0.5]
    simulated_phase_mask_U2 = [0, 1]
    fixed_phase_mask_U2 = [1, 0]

    s12_state_shelvings = []

    probe_trans = [[0, 3, 2], [0, 2, 2]]



############################### Just make sure s12 state is dealt with (this is shared for all dimensions)
    pulse_train_U2 = pulse_train_U2 + s12_state_shelvings
    fractions_U2 = fractions_U2 #+ s12_state_fractions
    fixed_phase_mask_U2 = fixed_phase_mask_U2 #+ s12_state_fixed_phases
    simulated_phase_mask_U2 = simulated_phase_mask_U2 #+ s12_state_simulated_phases

############################### Finished define the pulse sequence for a given dimension

    dimension = int(np.max(np.round(1/np.array(fractions_U1))))
    dimension_string = str(dimension)

    list_of_inits = [[[-1,3,-2],[0,4,-2],[1,4,0],[2,4,2]],
    [[-2,3,-1],[0,3,0],[1,4,0],[2,4,3]],
    [[-2,2,-1],[-1,3,0],[1,3,2],[2,4,2]],
    [[-2,3,-1],[-1,3,0],[0,3,2],[2,4,3]],
    [[-2,2,-1],[-1,3,0],[0,3,2],[1,3,3]]
    ]

    pulse_train = pulse_train_U1 + pulse_train_U2
    simulated_phase_mask = simulated_phase_mask_U1 + simulated_phase_mask_U2
    fixed_phase_mask = fixed_phase_mask_U1 + fixed_phase_mask_U2

    phase_shifts = zip(simulated_phase_mask,fixed_phase_mask)
    detunings = np.zeros(len(probe_trans))
    s12_level = initial_state[0][0]

    periodicity = 100 #us
    start_time = 0
    if full_phase_scan:
        # for full phase scans
        stop_time = 100
        start_time = 0
        time_step = 30*(1/dimension)
    else:
        # take just two points for a contrast measurement
        if dimension%2==0:
            time_step = (dimension)*100/(2*dimension)
            stop_time = 1+(dimension)*100/(2*dimension)
        else:
            time_step = (dimension+1)*100/(2*dimension)
            stop_time = 1+(dimension+1)*100/(2*dimension)

    F1PumpTime = 1.5 #us
    F1PumpReps = 60
    InitReps = 0
    threshold = 10

    pulse_program = "Qudit_ramsey_experiment_bused_ket0_in_D52"


    def findPiTime_withfit_plotPD(ramsey_wait,stop_time,start_time, time_step, threshold, pulse_program, script_functions):
        getGlobal, setGlobal, setScan, startScan, stopScan, getAllData, createTrace, closeTrace, plotPoint, scriptIsStopped, setEvaluation = script_functions
        createTrace('Qudit Ramsey, 0 state', 'Scan Data', xLabel=f'Pulse Time (us)')
        
        if Measure_U1_only:
            U1_only_bool_list = [False, True]
        else:
            U1_only_bool_list = [False]
        
        pulse_time = start_time

        file_names_list_all = []

        fluor_ave = []
        pulse_time_list = []
        file_names_list = []
        fluor_ave_U1 = []
        pulse_time_list_U1 = []
        file_names_list_U1 = []
        while pulse_time <= stop_time:
            for U1_only in U1_only_bool_list:
                if U1_only:
                    fractions = fractions_U1 + list(np.zeros(np.size(fractions_U2)))
                    do_calibrations = False
                else:
                    fractions = fractions_U1 + fractions_U2
                    do_calibrations = do_calibrations_0

                if scriptIsStopped():

                    break
        
                print('1')
                setEvaluation('Eval2')
                if pulse_time == 0:
                #if 1 == 1:
                    if do_calibrations:
                        delta = 0
                        ramsey_real_wait_time = 100
                        freq_offset = run_fast_calibration(script_functions,delta,ramsey_real_wait_time,target_transition=[0,2,0],f_offset_input=None,f_upper_input=None,ref_pitimes=None)
                        freq_upper = run_fast_calibration(script_functions,delta,ramsey_real_wait_time,target_transition=[-1,4,-3],f_offset_input=None,f_upper_input=None,ref_pitimes=None, check_3pt_coherence = True)
                setEvaluation('Eval3')
        
                f_offset = float(getGlobal('f_offset'))
                f_upper = float(getGlobal('f_upper'))
        
                pitime_n2 = float(getGlobal('pitime_n2')) # [-2, 4, -4]
                pitime_n1 = float(getGlobal('pitime_n1')) # [-2, 3, -3]
                pitime_0 = float(getGlobal('pitime_0')) # [2, 4, 2]
                pitime_p1 = float(getGlobal('pitime_p1')) # [2, 4, 3]
                pitime_p2 = float(getGlobal('pitime_p2')) # [2, 4, 4]
        
        
                init_trans = list_of_inits[s12_level+2]
                init_freqs_array = list(Get_1762_EOM_Freqs_an1an2(init_trans,f_offset,f_upper))
                init_times_array = list(Get_1762_PiTimes(init_trans,pitime_n2,pitime_n1,pitime_0,pitime_p1,pitime_p2))
                init_pulse_time = sum(init_times_array) 
        
                set_freq = list(Get_1762_EOM_Freqs_an1an2(probe_trans,f_offset,f_upper))
        
                initial_state_frequency = list(Get_1762_EOM_Freqs_an1an2(initial_state,f_offset,f_upper))
        
                pulse_train_frequencies = list(Get_1762_EOM_Freqs_an1an2(pulse_train,f_offset,f_upper))
        
                for idx, delta in enumerate(detunings):
                    set_freq[idx] = set_freq[idx] + delta
        
                for delta in detunings:
                    pulse_train_frequencies[0] = pulse_train_frequencies[0] + delta
        
                pi_time_initial_state = list(Get_1762_PiTimes(initial_state,pitime_n2,pitime_n1,pitime_0,pitime_p1,pitime_p2))
        
                pi_time_pulse_train = list(Get_1762_PiTimes(pulse_train,pitime_n2,pitime_n1,pitime_0,pitime_p1,pitime_p2))
        
                pi_times_probe_trans = list(Get_1762_PiTimes(probe_trans,pitime_n2,pitime_n1,pitime_0,pitime_p1,pitime_p2))
        
        
        
                print(pulse_train)
        
                F1_Pump_Time = getGlobal('F1_PumpTime')
                if not "%s"%F1_Pump_Time == "%s us"%F1PumpTime:
                    setGlobal("F1_PumpTime", F1PumpTime, "us")
        
                Init_reps = getGlobal('InitialisationReps')
                if not "%s"%Init_reps == InitReps:
                    setGlobal("InitialisationReps", InitReps, "")
        
                F1Pump_reps = getGlobal('F1_PumpReps')
                if not "%s"%F1Pump_reps == F1PumpReps:
                    setGlobal("F1_PumpReps", F1PumpReps, "")
        
                Init_PulseTime = getGlobal('Init_Shelving_Pulse_Time')
                if not "%s"%Init_PulseTime == init_pulse_time:
                    setGlobal("Init_Shelving_Pulse_Time", init_pulse_time, "us")
        
                SB_cooling_reps = getGlobal("Sideband_Cooling_Reps")
                if not "%s"%SB_cooling_reps == Side_band_cooling_reps:
                    setGlobal("Sideband_Cooling_Reps", Side_band_cooling_reps, "")
        
                corrected_pulse_train_times = []
        
                for i,_ in enumerate(fixed_phase_mask):
                    if np.isnan(_):
                        corrected_pulse_train_times.append(ramsey_wait)
                    else:
                        frac = fractions[i]
                        corrected_pulse_train_times.append(2*pi_time_pulse_train[i]*np.arcsin(np.sqrt(frac))/np.pi)
        
                print(corrected_pulse_train_times)
        
                Ramsey_wait_time_dummy = getGlobal('Ramsey_Wait_Time')
                if not "%s"%Ramsey_wait_time_dummy == "%s us"%sum(corrected_pulse_train_times):
                    setGlobal("Ramsey_Wait_Time", sum(corrected_pulse_train_times), "us")
        
                half_pi_time_dummy = getGlobal('Shelving_Pulse_Time')
                if not "%s"%half_pi_time_dummy == "%s us"%sum(corrected_pulse_train_times):
                    setGlobal("Shelving_Pulse_Time", sum(corrected_pulse_train_times), "us")
                print('2')
        
        # Phase Calculations ===================================================================
        
                phases = []
                for i,_ in enumerate(fixed_phase_mask):
                    if np.isnan(_):
                        phases.append(0)
                    else:
                        pi_phase = fixed_phase_mask[i]
                        real_phase = simulated_phase_mask[i]
                        x = (2 * np.pi * pulse_time * real_phase / periodicity)
                        y = pi_phase*np.pi
                        print(i,x,y)
                        phases.append(x + y)
                    
                full_phases = phases
                print(full_phases)
        #=======================================================================================


                init_row1 = [1,init_freqs_array, [0]*len(init_freqs_array),init_times_array, [0]*len(init_freqs_array),0]
                dummy_row2 = [2, [800],[0], [0.1], [0], 0]

                herald_row3 = [3, initial_state_frequency, [0], pi_time_initial_state, [0], 1]
                dummy_row4 = [4, [800],[0], [0.1], [0], 0]
                
                ramsey_row5 = [5, pulse_train_frequencies, full_phases, corrected_pulse_train_times, [1]*len(pulse_train_frequencies), 1]
                dummy_row6 = [6, [800],[0], [0.1], [0], 0]   
                
                rows = [init_row1, dummy_row2, herald_row3, dummy_row4, ramsey_row5, dummy_row6]
                
                row_number = 6
                for idx, freq in enumerate(set_freq):
                    row_number = row_number + 1
                    readout_row = [row_number, [freq], [0], pi_times_probe_trans[idx], [0], 1]  
                    rows.append(readout_row)
                    row_number = row_number + 1
                    dummy_row = [row_number, [800],[0], [0.1], [0], 0]
                    rows.append(dummy_row)
                    print(row_number,[freq],[pi_times_probe_trans[idx]])


        
  
        
                setScan(pulse_program)
                startScan(globalOverrides=list(), wait=True)
                stopScan()
        
                data = getAllData()
                herald_data = data['PMT Index 0'][1]
                print(herald_data)
                ket1_data = data['PMT Index 1'][1]
                print(ket1_data)
                ket2_data = data['PMT Index 2'][1]
                print(ket2_data)
        
                arrays = []
                for idx,herald_outcome in enumerate(herald_data):
                    if herald_outcome < threshold and ket1_data[idx] < threshold:
                        arrays.append(ket2_data[idx])
                mean_value = np.mean(np.array(arrays)<threshold)

#--------------------------------------------------------Raw_data_file_saving-----------------------------------------------------------------------------------------------

                year = datetime.datetime.now().strftime("%Y")
                month = datetime.datetime.now().strftime("%m")
                day = datetime.datetime.now().strftime("%d")

                file_path_today = f'C:\\Users\\ions\\Documents\\IonControl\\project_dir\\QuditLab\\{year}\\{year}_{month}\\{year}_{month}_{day}\\'
                destination_today = f'Z:\\Lab Data\\Qudit_Ramsey_raw_data\\Raw_data_copied\\{year}\\{year}_{month}\\{year}_{month}_{day}\\'
                pattern = f'C:\\Users\\ions\\Documents\\IonControl\\project_dir\\QuditLab\\{year}\\{year}_{month}\\{year}_{month}_{day}\\qudit_ramsey_scan_*'
                soruce_file_path = f'C:\\Users\\ions\\Documents\\IonControl\\project_dir\\QuditLab\\{year}\\{year}_{month}\\{year}_{month}_{day}\\qudit_ramsey_scan'

                def copy_file(src_file, dest_path):
                 
                    if not os.path.isfile(src_file):
                        print(f"Source file does not exist: {src_file}")
                        return
                 
                    dest_folder = os.path.dirname(dest_path)
                 
                    if not os.path.exists(dest_folder):
                        os.makedirs(dest_folder)
                        print(f"Created destination folder: {dest_folder}")
                 
                    shutil.copy2(src_file, dest_path)
                    print(f"File copied from {src_file} to {dest_path}")

                #if not file_names_list_all:
                matching_files = glob.glob(pattern)
                matching_files = sorted(matching_files, key=lambda t: os.stat(t).st_mtime)
                file_path = matching_files[-1]
                chunks = file_path.split('\\')
                print(chunks[:-1])        
                fname = chunks[-1]
                file_names_list_all.append(destination_today+fname)
                file_names_list.append(destination_today+fname)
                source_file = file_path_today + fname
                copied_file = destination_today + fname
                copy_file(source_file,copied_file)

                if np.isnan(mean_value):
                    mean_value = 1
                PD = mean_value
                print("PD is", PD)
                
                if U1_only:
                    fluor_ave_U1.append(PD)
                else:
                    fluor_ave.append(PD)
                
                if U1_only:
                    pulse_time_list_U1.append(pulse_time)
                else:
                    pulse_time_list.append(pulse_time)
                
                if not U1_only:
                    plotPoint(pulse_time, PD,'Qudit Ramsey, 0 state', plotStyle=2)
        
                freq_string = str(round(set_freq[0],3)).replace('.','p')
                if U1_only:
                    combined_data = zip(pulse_time_list_U1,fluor_ave_U1)
                else:
                    combined_data = zip(pulse_time_list,fluor_ave)
                
                # for bussed qudit Ramsey contrasts
                filename = f'Z:\\Lab Data\\Qudit_Ramsey_raw_data\\Raw_data_PD\\\Ramsey_qudit_WaitTime{ramsey_wait}us_d={dimension_string}_Cal_{do_calibrations_0}_U1only_{U1_only}_{repeat_ind}_{dt_string}.txt'
                # for bussed qubit Ramsey measurements
                #filename = f'Z:\\Lab Data\\Qudit_Ramsey_raw_data\\Raw_data_PD\\\Ramsey_experiment_{pulse_train_U1[:-1]}_{wtr}_us_{num_repeat}_{dt_string}.txt'
                with open(filename,'w') as file:
                    for x,y in combined_data:
                        file.write(f"{x},{y}\n")
                    file.write(f"{pulse_train}\n")
                    file.write(f"{[f_offset,f_upper]}\n")
                    file.write(f"{pi_time_pulse_train}\n")
                    file.write(f"{corrected_pulse_train_times}\n")
                    file.write(f"{fractions_U1}\n")
                    file.write(f"{fractions_U2}\n")
                    file.write(f"{probe_trans}\n")
                    file.write(f"{phases}\n")
                    if U1_only:
                        file.write(f"{file_names_list_U1}\n")
                    else:
                        file.write(f"{file_names_list}\n")
            pulse_time = pulse_time + time_step
        closeTrace('Qudit Ramsey, 0 state')


    # looping over different wait times
    for repeat_ind in range(repeat_experiments):
        for wait_time in ramsey_wait_times:
            print('this worked')
            setEvaluation('Eval3')
            pi_time = findPiTime_withfit_plotPD(wait_time,stop_time,start_time, time_step, threshold, pulse_program, script_functions)
            print(pi_time)
            setEvaluation('Eval2')
            #with open(output_file_pitimes,'a') as outfile:
            #    outfile.write(f'{pi_time}\n')

