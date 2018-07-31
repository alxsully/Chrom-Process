import sys, os
import tkinter
from tkinter import filedialog
import numpy as np
import matplotlib
matplotlib.use('TkAgg') 
from matplotlib import pyplot as plt
from collections import OrderedDict
import numpy as np
import pandas as pd
import time

class GCMS_Params():
    def __init__(self):
        self.param_indeces = OrderedDict([
            ('CHANNEL'              ,       ' '),
            ('MACHINE_TIMESTAMP'    ,       'Time [msec]'),
            ('PUFF_TIMESTAMP'       ,       'Time [msec]'),
            ('PUFF_COUNT'           ,       'Puffs [#]'),
            ('CHARGE_ENABLE'        ,       ' '),
            ('PUFF_ENABLE'          ,       ' '),
            ('FLUSH_ENABLE'         ,       ' '),
            ('VAC_ENABLE'           ,       ' '),
            ('FLOW_RATE_RAW'        ,       'Flow Rate [mLPM]'),
            ('FLOW_RATE_FILT'       ,       'Flow Rate [mLPM]'),
            ('FLOW_RATE_SET'        ,       'Flow Rate [mLPM]'),
            ('LIGHT_CELL_SET'       ,       'Voltage [mV]'),
            ('LIGHT_CELL_RAW'       ,       'ADC [ticks]')
        ])



# root = tkinter.Tk()
# root.withdraw()
# time.sleep(1)
# os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "python" to true' ''')
# FILE_PATHS = filedialog.askopenfilenames(title="Select all data files", initialdir=(os.path.expanduser(os.getcwd())))
# root.destroy()
# os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "Terminal" to true' ''')
# file_names = []
# print('\n')
# print('Loaded file {}'.format(FILE_PATHS[0].split('/')[-1]))
# print('\n')

# f = open(FILE_PATHS[0], errors='ignore')
# print(f.readline())
# f.close()

FILE_PATHS = ['Nicotine_1ppm_180712201814.csv']
FILE_PATHS += ['Nicotine_2,5ppm_180712202451.csv']
FILE_PATHS += ['Nicotine_5ppm_180712203135.csv']
FILE_PATHS += ['Nicotine_7,5ppm_180712203822.csv']
FILE_PATHS += ['Nicotine_10ppm_180712204456.csv']
TIME_TOLERANCE = 0.025
ABUNDANCE_OFFSET = 0.0
raw_data = {}
Quinoline = {}
Nicotine = {}
area_ratio = []

# Read file and store into dataframe
for i in range(len(FILE_PATHS)):
	raw_data[i] = pd.read_csv(
	        FILE_PATHS[i], 
	        sep=';',
	        index_col=False,
	        engine='c',
	        error_bad_lines=False)

	# Calculate time stamp differential in milliseconds (for integral calculation)
	raw_data[i]['t_diff'] = np.append(np.diff(raw_data[i]['RT(milliseconds)']), [0])

	# Create dataframe subset of raw_data for Quinoline and Nicotine
	Quinoline[i] = raw_data[i][['76', '102', '129']].copy()
	Nicotine[i] = raw_data[i][['84', '133', '162']].copy()

	# Calculate the mean of each compound across the three m/z ratio columns
	Quinoline[i]['mean'] = Quinoline[i].mean(axis=1)
	Nicotine[i]['mean'] = Nicotine[i].mean(axis=1)

	# Find the maximum value (peak) and its index of the mean
	max_value_mean_Q = Quinoline[i]['mean'].max(axis=0)
	max_value_mean_N = Nicotine[i]['mean'].max(axis=0)
	index_max_N = Nicotine[i].loc[Nicotine[i]['mean'] == max_value_mean_N].index
	index_max_Q = Quinoline[i].loc[Quinoline[i]['mean'] == max_value_mean_Q].index

	# Notmalize the mean to 100%
	Quinoline[i]['normalized_mean'] = 100 * (Quinoline[i]['mean'] / max_value_mean_Q)
	Nicotine[i]['normalized_mean'] = 100 * (Nicotine[i]['mean'] / max_value_mean_N)

	# Find the 
	value_time_high = raw_data[i].loc[index_max_N - 1, 'RT(minutes) - NOT USED BY IMPORT'].values * (1 + TIME_TOLERANCE)
	value_time_low = raw_data[i].loc[index_max_N - 1, 'RT(minutes) - NOT USED BY IMPORT'].values * (1 - TIME_TOLERANCE)
	nicotine_range = Nicotine[i].loc[(raw_data[i]['RT(minutes) - NOT USED BY IMPORT'] >= value_time_low[0]) & (raw_data[i]['RT(minutes) - NOT USED BY IMPORT'] <= value_time_high[0]), 'mean'] - (ABUNDANCE_OFFSET * max_value_mean_N)
	tdiff_range = raw_data[i].loc[(raw_data[i]['RT(minutes) - NOT USED BY IMPORT'] >= value_time_low[0]) & (raw_data[i]['RT(minutes) - NOT USED BY IMPORT'] <= value_time_high[0]), 't_diff']
	area_nicotine = (nicotine_range * tdiff_range).sum()
	print(area_nicotine)

	value_time_high = raw_data[i].loc[index_max_Q + 1, 'RT(minutes) - NOT USED BY IMPORT'].values * (1 + TIME_TOLERANCE)
	value_time_low = raw_data[i].loc[index_max_Q + 1, 'RT(minutes) - NOT USED BY IMPORT'].values * (1 - TIME_TOLERANCE)
	# quinoline_range = Quinoline.loc[(raw_data['RT(minutes) - NOT USED BY IMPORT'] >= value_time_low[0]) & (raw_data['RT(minutes) - NOT USED BY IMPORT'] <= value_time_high[0]), 'mean'] - (ABUNDANCE_OFFSET * max_value_mean_Q)
	quinoline_range = Quinoline[i].loc[(raw_data[i]['RT(minutes) - NOT USED BY IMPORT'] >= value_time_low[0]) & (raw_data[i]['RT(minutes) - NOT USED BY IMPORT'] <= value_time_high[0]), 'mean']
	tdiff_range = raw_data[i].loc[(raw_data[i]['RT(minutes) - NOT USED BY IMPORT'] >= value_time_low[0]) & (raw_data[i]['RT(minutes) - NOT USED BY IMPORT'] <= value_time_high[0]), 't_diff']
	print (quinoline_range)
	print (tdiff_range)
	area_quinoline = (quinoline_range * tdiff_range).sum()
	print(area_quinoline)

	area_ratio += [area_nicotine / area_quinoline]

print(area_ratio)

x = [1,2.5,5,7.5,10]
y = area_ratio
(a, b) = np.polyfit(x, y, 1)
print('y = {:0.05f} x + {:0.05f}'.format(a, b))
# # print(max_value_mean_Q, max_value_max_Q, max_value_min_Q)

fig1, axes1 = plt.subplots(1,1,sharex=True)
fig2, axes2 = plt.subplots(1,1,sharex=True)
axes1.scatter(x, y, marker='s')
axes1.plot(np.unique(x), np.poly1d((a,b))(np.unique(x)), linewidth=2.0, linestyle='--', color='r')
axes1.set(xlabel='ppm', ylabel='Area Ratio')
fig1.suptitle('Nicotine\ny = {:0.05f} x + {:0.05f}'.format(a, b))
# axes.plot(raw_data['RT(minutes) - NOT USED BY IMPORT'], Quinoline['normalized_mean'])
# axes.plot(raw_data['RT(minutes) - NOT USED BY IMPORT'], Nicotine['normalized_mean'])
axes2.plot(raw_data[0]['RT(minutes) - NOT USED BY IMPORT'], Quinoline[0]['mean'])
axes2.plot(raw_data[0]['RT(minutes) - NOT USED BY IMPORT'], Nicotine[0]['mean'])


plt.show()

