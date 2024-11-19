import csv
import datetime
import numpy as np
from scipy.signal import butter, filtfilt
from scipy.interpolate import interp1d

def importGCDC(filename, headerOnly=False):
    # Define structure for storing information
    info = {
        'deviceType': None,
        'serial': None,
        'date': None,
        'start': None,
        'uptime': None,
        'battery': None,
        'sampling': None,
        'deadband': None,
        'source': None
    }

    # Read file header
    try:
        with open(filename, 'r') as file:
            # Read lines from file, and remove the last character (newline)
            lines = file.readlines()
            # remove return carriage
            for i in range(len(lines)):
                lines[i] = lines[i][:-1]

            # Extract information from header
            info['deviceType'] = lines[0].split(',')[3][1:]
            info['serial'] = lines[1].split(',')[4][1:].replace(' ', '')
            info['date'] = datetime.datetime.strptime(lines[2].split(',')[1][1:], '%Y-%m-%d').date()
            info['start'] = datetime.datetime.strptime(lines[2].split(',')[1][1:] + '-' + lines[2].split(',')[2][1:], '%Y-%m-%d-%H:%M:%S.%f')
            info['uptime'] = float(lines[3].split(',')[1])
            info['battery'] = float(lines[3].split(',')[4])

            units = {
                'SN': [
                    'SN:CCDC1002DA64743',
                    'SN:CCDC1002FD3B45E',
                    'SN:CCDC100287EA8D2',
                    'SN:CCDC1002617A2D8',
                    'SN:CCDC10027233E8E',
                    'SN:CCDC1002DF160DF',
                    'SN:CCDC1002A9BEA64',
                    'SN:CCDC1002B202393',
                    'SN:CCDC100247C9CE5',
                    'SN:CCDC1002B999430',
                    'SN:CCDC1002B561590',
                    'SN:CCDC10211F02908',
                    'SN:CCDC102146158B8',
                    'SN:CCDC10219779497',
                    'SN:CCDC1021AF3B8A3',
                    'SN:CCDC10211C00096',
                    'SN:CCDC1021DFA9251'
                ],
                'type': [
                    'USB Accelerometer X2-5',
                    'USB Accelerometer X2-5',
                    'USB Accelerometer X2-5',
                    'USB Accelerometer X2-5',
                    'USB Accelerometer X2-5',
                    'USB Accelerometer X2-5',
                    'USB Accelerometer X2-5',
                    'USB Accelerometer X2-5',
                    'USB Accelerometer X2-5',
                    'USB Accelerometer X2-5',
                    'USB Accelerometer X2-5',
                    'Multifunction Extended Life (MEL) Data Logger',
                    'Multifunction Extended Life (MEL) Data Logger',
                    'Multifunction Extended Life (MEL) Data Logger',
                    'Multifunction Extended Life (MEL) Data Logger',
                    'Multifunction Extended Life (MEL) Data Logger',
                    'Multifunction Extended Life (MEL) Data Logger'
                ]
            }

            if info['serial'] in units['SN'][:2]:
                info['sampling'] = float(lines[4].split(',')[1])
            elif info['serial'] in units['SN'][2:]:
                info['sampling'] = float(lines[4].split(',')[2])
            else:
                raise ValueError("Unknown serial number")

            info['deadband'] = float(lines[5].split(',')[1])
            info['source'] = filename

    except FileNotFoundError:
        print("Inexistent file")
        return None, None, None, None, None

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
    if len(ax) > 24:
        ax = filtfilt(b, a, ax)
        ay = filtfilt(b, a, ay)
        az = filtfilt(b, a, az)

    newtime = np.arange(time[0], time[-1] + 1 / info['sampling'], 1 / info['sampling'])

    interpolate = interp1d(time, ax, bounds_error=False, fill_value="extrapolate")
    ax = interpolate(newtime)
    interpolate = interp1d(time, ay, bounds_error=False, fill_value="extrapolate")
    ay = interpolate(newtime)
    interpolate = interp1d(time, az, bounds_error=False, fill_value="extrapolate")
    az = interpolate(newtime)
    
    time = newtime

    print('Done.\n')

    return info, time, np.column_stack((ax, ay, az))

# Test the function
# info, time, acc = importGCDC("your_file.csv")
