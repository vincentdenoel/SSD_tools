import csv
import pandas as pd
import pandas as pd
import datetime
import itertools
import itertools
import numpy as np
from scipy.signal import butter, filtfilt
from scipy.interpolate import interp1d

def importGCDC(filename, headerOnly=False):
    # Define structure for storing information
    info = {}


    try:
        # Read file header (first 30 lines)
        nlines = 30
        lines = []
        with open(filename, 'r') as file:
            for i, line in enumerate(file):
                lines.append(line.strip())
                if i >= nlines:
                    break           
            if not lines:
                print("File is empty or does not contain enough lines.")
                return None, None
            info['serial'] = (lines[1].split(',')[4]).replace(' ', '')
            info['source'] = filename
            unit_type = get_unit_type_by_serial(info['serial'])
            
            if unit_type == 'USB Accelerometer X2-5':
                # TO DO !!
                info['deviceType'] = 'USB Accelerometer X2-5'
                info['sampling'] = float(lines[4].split(',')[1])
                info['']
            
            elif unit_type == 'Yellow Multifunction Extended Life (MEL) Data Logger':
                info['deviceType'] = 'Yellow Multifunction Extended Life (MEL) Data Logger'
                info['sampling'] = float(lines[4].split(',')[2])
                info['sensitivity'] = float(lines[4].split(',')[9]) # LSB/g
                info['start-date'] = lines[2].split(',')[1]
                info['start-time'] = lines[2].split(',')[2]
                info['uptime'] = float(lines[3].split(',')[1])
                info['battery'] = float(lines[3].split(',')[4])
                info['eol'] = float(lines[3].split(',')[7])
                info['deadband'] = float(lines[5].split(',')[1]) # counts
                info['deadband-timeout'] = float(lines[6].split(',')[1]) # seconds
            
                # find a line containing "Time, Ax, Ay, Az"
                skiplines = 1
                for line in lines:
                    skiplines += 1
                    if 'Time, Ax, Ay, Az' in line:
                        break

            elif unit_type == 'Black Multifunction Extended Life (MEL) Data Logger':
                info['deviceType'] = 'Black Multifunction Extended Life (MEL) Data Logger'
                info['sampling'] = float(lines[6].split(',')[2]) # Hz
                info['sampling-gps'] = float(lines[7].split(',')[2]) # Hz                
                info['sensitivity'] = float(lines[6].split(',')[9]) # LSB/g
                info['start-date'] = lines[2].split(',')[1]
                info['start-time'] = lines[2].split(',')[2]
                info['uptime'] = float(lines[3].split(',')[1])
                info['battery'] = float(lines[3].split(',')[4])
                info['eol'] = float(lines[3].split(',')[7])
                info['deadband'] = float(lines[4].split(',')[1]) # counts
                info['deadband-timeout'] = float(lines[5].split(',')[1]) # seconds

                # find a line containing "Time, Ax, Ay, Az, T"
                skiplines = 1
                for line in lines:
                    skiplines += 1
                    if 'Time, Ax, Ay, Az, T' in line:
                        break

    except FileNotFoundError:
        print("Inexistent file")
        return None, None

    if headerOnly:
        return info, None, None, None, None

    # Read data from file
    print('\nReading file: {} ...\n  on device (serial): {}\n  Start date & time: {}\n'.format(
        filename, info['serial'], info['start'].strftime('%a %d-%b-%Y %H:%M:%S')))

    with open(filename, 'r') as file:
        reader = csv.reader(file)
        for i in range(8):
            next(reader)  # Skip first 7 lines (header)

        # Read data into arrays
        data = list(reader)
        time = [float(row[0]) for row in data[:-1]]
        ax = [float(row[1]) for row in data[:-1]]
        ay = [float(row[2]) for row in data[:-1]]
        az = [float(row[3]) for row in data[:-1]]

    # Resample for better accuracy
    b, a = butter(8, 0.90)
    if len(df['acc-x']) > 24:
        ax = filtfilt(b, a, ax)
        ay = filtfilt(b, a, ay)
        az = filtfilt(b, a, az)

    # check median diff(time) and find large gaps in the time series
    difftime = np.diff(time)
    if len(difftime) > 0:
        median_diff = np.median(difftime)
        if np.max(np.abs(difftime)) / median_diff > 3:
            # find where the gaps are
            gaps = np.where(np.abs(np.diff(time)) > median_diff * 3)[0]
            # shift all values behind the gap by the median difference
            for gap in gaps:
                print('Warning: it seems there was a gap in the data. Correction is applied.')
                time[int(gap)+1:] += (median_diff - difftime[int(gap)])
    
    # find negative dt (glitches), replace them with the median difference
    difftime = np.diff(time)
    neg_diffs = np.where(difftime < 0.1*median_diff)[0]
    if len(neg_diffs) > 0:
        print('Warning: it seems there were glitches in the data. Correction is applied.')
        for neg_diff in neg_diffs:
            time[int(neg_diff)+1:] += (median_diff - difftime[int(neg_diff)])
    
    # Resample the data to a regular time grid
    dt = 1 / info['sampling']       
    newtime = np.arange(time[0], time[-1] + dt, dt)


    ax = np.interp(newtime, time, ax)
    ay = np.interp(newtime, time, ay)
    az = np.interp(newtime, time, az)

    # Shift start time if required
    if start_date is not None:
        start_datetime = pd.to_datetime(start_date)
        start_timestamp = start_datetime.timestamp()
        time_shift = start_timestamp - time[0]
        newtime += time_shift

    # Convert the first column from float to datetime
    df_resampled = pd.DataFrame({
        'time': pd.to_datetime(newtime, unit='s'),
        'ax': ax,
        'ay': ay,
        'az': az
    })

    # Set the time column as the index
    df_resampled.set_index('time', inplace=True)
    print('Done.\n')

    return info, df_resampled

# Test the function
# info, df = importGCDC("your_file.csv")

def get_unit_type_by_serial(serial_number):
    for unit in units:
        if unit['sn'] == serial_number:
            return unit['type']
    raise ValueError("Unknown serial number")

import numpy as np
import itertools

def detect_and_correct_glitches(ax, ay, az, jump_factor=10, improvement_factor=5, force_gravity='z'):
    """
    Efficient detection and correction of 1-line axis-swap glitches in GCDC accelerometer data.
    Works directly on numpy arrays (ax, ay, az).
    """
    charact_val = np.mean(np.sqrt(ax**2 + ay**2 + az**2))
    acc = np.column_stack((ax, ay, az))
    
    # Check the gravity is in the right direction for the first sample
    tol = 0.5 # tolerance (in g) to decide if the axis is aligned with gravity
    if force_gravity in ['x', 'y', 'z']:
        axis_index = {'x': 0, 'y': 1, 'z': 2}[force_gravity]
        if np.abs(np.abs(acc[0, axis_index]) - 1.0) > tol: # acceleration along axis_index is too far from 1g
            # find if any other index is close to 1g
            other_indices = [i for i in range(3) if i != axis_index]
            if  np.abs(np.abs(acc[0, other_indices[0]]) - 1.0) < tol:
                p = np.zeros(3)
                p[axis_index] = other_indices[0]
                p[other_indices[1]] = other_indices[1]
                p[other_indices[0]] = axis_index
                acc = acc[:, p]
            elif np.abs(np.abs(acc[0, other_indices[1]]) - 1.0) < tol:
                p = np.zeros(3).astype(int)
                p[axis_index] = other_indices[1]
                p[other_indices[0]] = other_indices[0]
                p[other_indices[1]] = axis_index
                acc = acc[:, p]

    # Compute vectorized jump magnitudes
    diffs = np.diff(acc, axis=0)
    jump_norms = np.linalg.norm(diffs, axis=1)
    median_jump = np.median(jump_norms)
    #threshold = jump_factor * median_jump
    threshold = 0.1 * charact_val

    # Candidate glitch indices (large jumps only)
    glitch_idx = np.where(jump_norms > threshold)[0]
    if len(glitch_idx) == 0:
        print("✅ No potential glitches detected.")
        return acc[:, 0], acc[:, 1], acc[:, 2]

    print(f"🔍 Found {len(glitch_idx)} potential glitches out of {len(acc)} samples.")


    # Define transformations
    perms = list(itertools.permutations(range(3)))  # (0,1,2), (0,2,1), ...
    signs = list(itertools.product([-1, 1], repeat=3))
    
    for i in glitch_idx:
        jump = np.linalg.norm(acc[i+1] - acc[i])

        best_jump = jump
        best_transform = None

        for p in perms:
            for s in signs:
                transformed = np.array([
                    s[0] * acc[i+1, p[0]],
                    s[1] * acc[i+1, p[1]],
                    s[2] * acc[i+1, p[2]],
                ])
                new_jump = np.linalg.norm(acc[i] - transformed)
                if new_jump < best_jump:
                    best_jump = new_jump
                    best_transform = (p, s)

        if best_transform and jump / best_jump > improvement_factor:
            p, s = best_transform
            # correct all remaining data points
            acc[i+1:] = acc[i+1:, p] * s
            """
            for j in range(i+1, len(acc)):
                acc[j] = [
                    s[0] * acc[j, p[0]],
                    s[1] * acc[j, p[1]],
                    s[2] * acc[j, p[2]],
                ]"""

            print(f"\033[93m⚠️  SIGN GLITCH FIXED at line {i:6d} \033[0m")

    return acc[:, 0], acc[:, 1], acc[:, 2]
