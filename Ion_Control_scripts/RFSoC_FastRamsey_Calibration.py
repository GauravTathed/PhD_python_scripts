import sys
sys.path.append(r'C:\Users\ions\Documents\IonControl\project_dir\QuditLab\config\Scripts')

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
       
import time
from Functions_Data import *
from Functions_AWG import *
from Functions_Measurement import *

transition_strengths = np.loadtxt(
    'Z:\\Lab Data\\Phase_and_freq_correction_180Hz\\Transition_strengths_4p216.txt', delimiter=','
)

def compute_pi_times(pi_t, ref_transitions = np.array([
        transition_strengths[23, 0], transition_strengths[14, 0],
        transition_strengths[17, 4], transition_strengths[16, 4], transition_strengths[15, 4]
    ])):
    transition_strengths = np.loadtxt(
        'Z:\\Lab Data\\Phase_and_freq_correction_180Hz\\Transition_strengths_4p216.txt', delimiter=','
    )
    transition_strengths[transition_strengths == 0] = np.nan

    # pi_t = np.array([19.470, 35.554, 41.166, 30.108, 39.326])
    strengths = ref_transitions

    factors = np.array(pi_t) * strengths
    Fs = [1, 2, 3, 4]
    row_labels = [[i, i - j] for i in Fs for j in range(2 * i + 1)]
    col_labels = [-2, -1, 0, 1, 2]

    pi_times = np.zeros((24, 5))
    for i in range(np.shape(transition_strengths)[0]):
        for j in range(np.shape(transition_strengths)[1]):
            if not np.isnan(transition_strengths[i, j]):
                delta_m = (row_labels[i][1] - col_labels[j]) + 2
                pi_times[i, j] = factors[delta_m] / transition_strengths[i, j]
    
    return pi_times

def get_pi_times(transitions, pi_t,ref_transitions = np.array([
        transition_strengths[23, 0], transition_strengths[14, 0],
        transition_strengths[17, 4], transition_strengths[16, 4], transition_strengths[15, 4]
    ])):
    matrix = compute_pi_times(pi_t, ref_transitions)
    pi_times_list = []
    col_labels = [-2, -1, 0, 1, 2]
    for transition in transitions:
        row_label = [transition[1],transition[2]]
        Fs = [1,2,3,4]
        states = []
        for i in Fs:
            for j in range(2*i+1):
                mF = i-j
                states.append([i,mF])
    
        row_labels = states
    
        col_label = transition[0]
    
        # Find the index of the row label
        row_index = next((i for i, label in enumerate(row_labels) if label == row_label), None)
        # Find the index of the column label
        col_index = col_labels.index(col_label)
        
        if row_index is not None and col_index in range(len(col_labels)):
            pi_times_list.append(matrix[row_index, col_index])
        else:
            pi_times_list.append(np.nan)

    return pi_times_list


def run_fast_calibration(script_functions,delta,real_wait_time,target_transition=[0,2,0],f_offset_input=None,f_upper_input=None,ref_pitimes=None,diff_pass=0.3,check_3pt_coherence=False):
    getGlobal, setGlobal, setScan, startScan, stopScan, getAllData, createTrace, closeTrace, plotPoint, scriptIsStopped, setEvaluation = script_functions

    import time
    import json, urllib.request
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

    dt_string = datetime.datetime.now().strftime("%Y%m%d_%H%M")
    year = datetime.datetime.now().strftime("%Y")
    month = datetime.datetime.now().strftime("%m")
    day = datetime.datetime.now().strftime("%d")
    pattern = "Z:\\Lab Data\\D52_Calibration_Ba137\\Fast_calibration_ion_control_raw_data\\fast_calibration_*"
    
    Side_band_cooling_reps = 0
    real_wait_time_bool = True
    f_offset_index = [0,2,0]
    f_upper_index = [-1,4,-3]
    pitime_ref_index = [[-2,4,-4],[-2,3,-3],[2,4,2],[2,4,3],[2,4,4]]

    if f_offset_input and f_upper_input:
        [f_offset,f_upper] = [f_offset_input,f_upper_input]
    else:
        f_offset = float(getGlobal('f_offset'))
        f_upper = float(getGlobal('f_upper'))
        print('getting ion control info',f_upper,f_offset)

    if ref_pitimes:
        pass
    else:
        pitime_n2 = float(getGlobal('pitime_n2'))
        pitime_n1 = float(getGlobal('pitime_n1'))
        pitime_0 = float(getGlobal('pitime_0'))
        pitime_p1 = float(getGlobal('pitime_p1'))
        pitime_p2 = float(getGlobal('pitime_p2'))
        ref_pitimes = [pitime_n2,pitime_n1,pitime_0,pitime_p1,pitime_p2]
        
    F1PumpTime = 0.5 #us
    F1PumpReps = 50
    InitReps = 0
    fs = 4e9
    threshold = 9
    
    list_of_inits = [[[-1,3,-2],[0,4,-2],[1,4,0],[2,4,2]],
                    [[-2,3,-1],[0,3,0],[1,4,0],[2,4,3]],
                    [[-2,2,-1],[-1,3,0],[1,3,2],[2,4,2]],
                    [[-2,3,-1],[-1,3,0],[0,3,2],[2,4,3]],
                    [[-2,2,-1],[-1,3,0],[0,3,2],[1,3,3]]
                    ]

    probe_transitions = [[target_transition]]
    for probe_trans in probe_transitions:

        
        s12_level = probe_trans[0][0]
        
        periodicity = 100 #us

        if check_3pt_coherence:
        
            wait_times = [periodicity/4,periodicity/2,periodicity*3/4]
            # Phase Calculations ===================================================================

            phases = [np.pi/2, np.pi, np.pi*3/2]
            #phases = [0.0*np.pi,0.1*np.pi,0.2*np.pi,0.3*np.pi,0.4*np.pi,0.5*np.pi,0.6*np.pi,0.7*np.pi,0.8*np.pi,0.9*np.pi,1.0*np.pi]
        else:
            wait_times = [periodicity/4,periodicity*3/4]
            # Phase Calculations ===================================================================

            phases = [np.pi/2, np.pi*3/2]

        #=======================================================================================

        init_trans = list_of_inits[s12_level+2]
        init_freqs_array = list(Get_1762_EOM_Freqs_an1an2(init_trans,f_offset,f_upper))
        init_times_array = list(get_pi_times(init_trans,ref_pitimes))
        init_pulse_time = sum(init_times_array) 
        
        set_freq = list(Get_1762_EOM_Freqs_an1an2(probe_trans,f_offset,f_upper))
        set_freq = [set_freq[0] + delta]

        freq_og = set_freq[0]

        print('set freq', set_freq)

        f_target_dummy = getGlobal('f_target')
        if not "%s"%f_target_dummy == set_freq[0]:
            setGlobal("f_target", set_freq[0], "")
               
      
        pi_times= list(get_pi_times(probe_trans,ref_pitimes))
        times = []
        for i in range(len(pi_times)):
            times.append(2*pi_times[i]*np.arcsin(np.sqrt(1/(len(pi_times)+1-i)))/np.pi)
        
        
        print(times)
        full_freqs = set_freq + [0] + list(np.flip(set_freq))
        print(full_freqs)
        
        LT_dummy = getGlobal('LineTriggerBoolean')
        if not "%s"%LT_dummy == 1:
            setGlobal("LineTriggerBoolean", 1, "")

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

        Hareld_pulse_time_dummy = getGlobal('Hareld_pulse_time')
        if not "%s"%Hareld_pulse_time_dummy == pi_times[0]:
                setGlobal("Hareld_pulse_time", pi_times[0], "us")
        
        pulse_program = "Qudit_ramsey_experiment_fast_calibration"
        
        def perform_fast_Ramsey_cal(f_offset,f_upper,real_wait_time, threshold, pulse_program, script_functions):
            getGlobal, setGlobal, setScan, startScan, stopScan, getAllData, createTrace, closeTrace, plotPoint, scriptIsStopped, setEvaluation = script_functions
            createTrace('Fast Calibration Points', 'Script Data', xLabel=f'Pulse Time (us)')
            wait_time_num = 0
            fluor_ave = []
            pulse_time_list = []
            file_names_list = []
            
            if check_3pt_coherence:
                ramsey_points=3
            else:
                ramsey_points=2
                
            while wait_time_num<ramsey_points:
                if scriptIsStopped():
                    break
                
                set_freq = float(getGlobal('f_target'))
                print('set freq', set_freq)

                full_freqs = [float(set_freq),800,float(set_freq)]
       
                print('1')
        
        #Pulse Times ==============================================================
                T_wait = real_wait_time
                full_times = list(times) + [T_wait] + list(np.flip(times))
        #Phases =================================================================
                full_phases = [0,0,phases[wait_time_num]]
                print(full_phases)
        #==============================================================
                Ramsey_wait_time_dummy = getGlobal('Ramsey_Wait_Time')
                if not "%s"%Ramsey_wait_time_dummy == "%s us"%sum(full_times):
                    setGlobal("Ramsey_Wait_Time", sum(full_times), "us")
        
                half_pi_time_dummy = getGlobal('Shelving_Pulse_Time')
                if not "%s"%half_pi_time_dummy == "%s us"%sum(full_times):
                    setGlobal("Shelving_Pulse_Time", sum(full_times), "us")
                print('2')


                if scriptIsStopped():
                    break
 
                init_row1 = [1, init_freqs_array, [0]*len(init_freqs_array), init_times_array, [0]*len(init_times_array), 0]    
                dummy_row2 = [2, [800], [0], [0.01], [0], 0]
                herald_row3 = [3, [set_freq], [0]*len([set_freq]), [pi_times[0]], [0]*len([set_freq]), 1]
                dummy_row4 = [4, [800], [0], [0.01], [0], 0]
                ramsey_row5 = [5, full_freqs, full_phases, full_times, [1,1,1], 1]
                dummy_row6 = [6, [800], [0], [0.01], [0], 0]


                rows = [init_row1, 
                        dummy_row2, 
                        herald_row3, 
                        dummy_row4, 
                        ramsey_row5, 
                        dummy_row6]
                print(rows)

                resp = upload_rows(rows, time_unit="us")
                print(json.dumps(resp, indent=2))

                seq_num = 6
                Table_length_dummy =  getGlobal('Table_length')
                if not "%s"%Table_length_dummy == seq_num:
                    setGlobal("Table_length", seq_num, "")
        
                
        
                #total_time = sum([piover2_time,pulse_time,piover2_time])+30
                #PulseTime_Dummy = getGlobal('Shelving_Pulse_Time').magnitude
                #if not "%s"%PulseTime_Dummy == "%s us"%total_time:
                #    setGlobal("Shelving_Pulse_Time", total_time, "us")
                #time.sleep(0.2)
        
                setScan(pulse_program)
                startScan(globalOverrides=list(), wait=True)
                stopScan()

                #new way of getting the right PMT counts for thresholding/heralding and plotting
                data = getAllData()
                herald_data = data['PMT Index 0'][1]
                print(herald_data)
                ket1_data = data['PMT Index 1'][1]
                print(ket1_data)


                arrays = []
                for idx,herald_outcome in enumerate(herald_data):
                    if herald_outcome < threshold:
                        arrays.append(ket1_data[idx])
                mean_value = np.mean(np.array(arrays)<threshold)
        
                #matching_files = glob.glob(pattern)
                #matching_files = sorted(matching_files, key=lambda t: os.stat(t).st_mtime)
                #file_path = matching_files[-1]
                #chunks = file_path.split('\\')
                #print(chunks[-1])        
                #fname = chunks[-1]
                #file_names_list.append(fname)
                #print("=========================================")
        
                #arrays = []
                #with open(file_path, 'r') as file:
                #    print(file_path)
                #    for line in file:
                #        data = json.loads(line)
                #        if data[0]["0"][0] < threshold:
                #            arrays.append(data[0]["0"][1])
                #mean_value = np.mean(np.array(arrays)<threshold)
                #print("=========================================")
        
                if np.isnan(mean_value):
                    mean_value = 1
                PD = mean_value
                print("PD is", PD)
                fluor_ave.append(PD)
                wait_time_ramsey = wait_times[wait_time_num]
                pulse_time_list.append(wait_times[wait_time_num])
        
                wait_time_num = wait_time_num + 1

                plotPoint(wait_time_ramsey, PD,'Fast Calibration Points', plotStyle=2)
                
            closeTrace('Fast Calibration Points')

            def unitary_product_with_ket_zero(pi_time, T,T_wait,T2_time, delta, t_wait_T):
                omega = np.pi / (pi_time * 1e-6)
                t_pi_over_2 = (pi_time * 1e-6) / 2
                T = T*1e-6
                T_wait_s = T_wait*1e-6
                def U1(omega, t):
                    H = np.array([[0, omega / 2],
                                  [omega / 2, delta]], dtype=complex)
                    U = expm(-1j * H * t)
                    return U

                def U_wait(omega, t):
                    H = np.array([[0, 0],
                                  [0, delta]], dtype=complex)
                    U = expm(-1j * H * t)
                    return U

                def U2(omega, t):
                    H = np.array([[0, omega * np.exp(-1j * 2 * np.pi * t_wait_T) / 2],
                                  [omega * np.exp(1j * 2 * np.pi * t_wait_T) / 2, delta]], dtype=complex)
                    U = expm(-1j * H * t)
                    return U

                U1 = U1(omega, t_pi_over_2)
                U_wait = U_wait(omega, T_wait_s)
                U2 = U2(omega, t_pi_over_2)

                product = U2 @ U_wait @ U1

                ket_zero = np.array([1, 0], dtype=complex)
                result = product @ ket_zero
                state_pop = np.abs(result)**2
				
                #T2_star = T2_time
                #decay_factor = np.exp(- T_wait**2 / T2_star**2)
                #state_pop *= decay_factor
				
                return state_pop

            def call_unitary_product_twice(pi_time, T,T_wait,T2_time, delta):
                result1 = unitary_product_with_ket_zero(pi_time, T,T_wait,T2_time, delta, 1/4)
                result2 = unitary_product_with_ket_zero(pi_time, T,T_wait,T2_time, delta, 3/4)
                return result1[0], result2[0]

            def fitting_function(delta, pi_time, T,T_wait,T2_time, target_results):
                computed_result1, computed_result2 = call_unitary_product_twice(pi_time, T,T_wait,T2_time, delta[0])
                target_result1, target_result2 = target_results
                error = (computed_result1 - target_result1)**2 + (computed_result2 - target_result2)**2
                return error
            
            pi_time = pi_times[0]
            initial_delta_guess = 0
            result1 = fluor_ave[0]
            
            if check_3pt_coherence:
                result2 = fluor_ave[2]
            else:
                result2 = fluor_ave[1]
            
            target_results = (result1 , result2)
            T = periodicity
            T2_time = None
            result = minimize(fitting_function, initial_delta_guess, args=(pi_time, T,T_wait,T2_time, target_results),method='Nelder-Mead')
            optimized_delta = result.x[0]/[2*np.pi]
            print(optimized_delta)

            freq_string = str(round(freq_og,6)).replace('.','p')
            freq_string_prev = str(round(set_freq,6)).replace('.','p')
            combined_data = zip(pulse_time_list,fluor_ave)
            
			
            upload_issues = result1==0 or result2==0 or result1==1 or result2==1
            ion_issues = result1==1 and result2==1
            if upload_issues or ion_issues:
                new_set_frequency = [set_freq]
            else:
                new_set_frequency = [np.abs(set_freq) - optimized_delta*1e-6]


                if abs(result1-result2)<0.45:
                    filename = f'Z:\\Lab Data\\D52_Calibration_Ba137\\Fast_calibration_data\\Fast_calibration_TS_{freq_string}_{T_wait}_{dt_string}.txt'
                    with open(filename,'w') as file:
                        for x,y in combined_data:
                            file.write(f"{x},{y}\n")
                        file.write(f"[{freq_og}],[{new_set_frequency[0][0]}]\n")
                        file.write(f"{pi_times},{ref_pitimes}\n")
                        file.write(f"{phases}\n")
                        file.write(f"{file_names_list}\n")

                    filename = f'Z:\\Lab Data\\D52_Calibration_Ba137\\Fast_calibration_fits\\Fast_calibration_TS_{freq_string_prev}_{T_wait}_{dt_string}.txt'
                    with open(filename,'w') as file:
                        for x,y in combined_data:
                            file.write(f"{x},{y}\n")
                        file.write(f"[{set_freq}],[{new_set_frequency[0][0]}]\n")
                        file.write(f"{pi_times},{ref_pitimes}\n")
                        file.write(f"{phases}\n")
                        file.write(f"{file_names_list}\n")

                else:
                    #dt_string_current = datetime.datetime.now().strftime("%Y%m%d_%H%M")
                    filename = f'Z:\\Lab Data\\D52_Calibration_Ba137\\Fast_calibration_fits\\Fast_calibration_TS_{freq_string_prev}_{T_wait}_{dt_string}.txt'
                    with open(filename,'w') as file:
                        for x,y in combined_data:
                            file.write(f"{x},{y}\n")
                        file.write(f"[{set_freq}],[{new_set_frequency[0][0]}]\n")
                        file.write(f"{pi_times},{ref_pitimes}\n")
                        file.write(f"{phases}\n")
                        file.write(f"{file_names_list}\n")

            return new_set_frequency[0],result1,result2
        
        setEvaluation('Eval3')
		
        iterator = 0.5
        mean_count_diff = diff_pass+1
        
        while mean_count_diff > diff_pass:
            if iterator>3:
                # print('Had upload or ion issues repeatedly.')
                # print('Nonsense to stop script', sdfasd)
                raise Exception("Stopping script due to upload or ion issues.")
                break   
            else:
                pass
                    
            new_set_frequency,fluor_ave0,fluor_ave1 = perform_fast_Ramsey_cal(f_offset,f_upper,real_wait_time, threshold, pulse_program, script_functions,)
            
            upload_issues = fluor_ave0==1 or fluor_ave1==1
            ion_issues = fluor_ave0==1 and fluor_ave1==1
            if ion_issues:
                iterator = iterator + 1
            else:
                iterator = abs(fluor_ave0-fluor_ave1)
                mean_count_diff = abs(fluor_ave0-fluor_ave1)
                f_target_dummy = getGlobal('f_target')
                if not "%s"%f_target_dummy == new_set_frequency:
                    setGlobal("f_target", new_set_frequency, "")
                            
        if target_transition==[0,2,0]:
            f_offset_dummy = getGlobal('f_offset')
            if not "%s"%f_offset_dummy == new_set_frequency:
                    setGlobal("f_offset", new_set_frequency, "")

        elif target_transition==[-1,4,-3]:
            f_upper_dummy = getGlobal('f_upper')
            if not "%s"%f_upper_dummy == new_set_frequency:
                    setGlobal("f_upper", new_set_frequency, "")
        else:
            pass
		
        print(new_set_frequency)
        setEvaluation('Eval3')
        return new_set_frequency[0]


#script_functions = (getGlobal, setGlobal, setScan, startScan, stopScan, getAllData, createTrace, closeTrace, plotPoint, scriptIsStopped, setEvaluation)
new_set_frequency = run_fast_calibration(script_functions, 0, 25,target_transition=[0,2,0],f_offset_input=None, f_upper_input=None,ref_pitimes=None,diff_pass=0.2,check_3pt_coherence=True)
'''deltas = np.linspace(-0.001,0.001,5)
for delta in deltas:
    #delta = 0
    new_set_frequency = run_fast_calibration(script_functions, 0, 100,target_transition=[-1,4,-3],f_offset_input=None, f_upper_input=None,ref_pitimes=None)
    f_upper_input1=new_set_frequency
print(new_set_frequency)'''