import pandas as pd
import numpy as np
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt
import scipy as sp
#direct questions to AlexSolivais@gmail.com

#imports chromatogram.csv generated using arbitrary xcalibur method, returns only columns for
#selected masses(default is SIM_Nicotine selected masses)
def import_SIM_chrom(file, m_z = ['76','84','102','129','133','162']):
	rt_label = ['RT(minutes) - NOT USED BY IMPORT']
	import_list = rt_label + m_z
	return pd.read_csv(file, sep=";", usecols=import_list) 

#renames RT column
def rename_RT(df):
	df.rename(
	columns = {'RT(minutes) - NOT USED BY IMPORT': 'RT'}, inplace = True)
	return df

#sums intensities across various masses and adds column 'tot_I'
#specific to masses from SIM_Nicotine method
def sum_I(df):
	df['tot_I'] = np.sum(df, axis=1)
	df.drop(columns = ['76','84','102','129','133','162'], inplace=True)
	return df

#combines above three functions
def chrom_import(file):
	nic = sum_I(rename_RT(import_SIM_chrom(file)))	
	return nic

#Peak Find Function
'''create moving baseline, where lag = length of moving baseline and threshold = std devs above baseline 
influence variable lowers weighting for prior signal hits
baseline avg stored in avgFilter , std dev stored in stdFilter
hits when an intensity (y) > threshold*std dev + baseline
outputs bool array where True values indicate signals in peak
also, runs the baseline backwards to capture the back half of peaks 
draws pretty heavily on https://bit.ly/2NNSKEL'''
def baseline(df, lag = 15, rear_lag = 20, threshold = 1.9, rear_thresh=3, influence = 0.3):
	y = np.array(df['tot_I']) #y is numpy array of intensities
	filterY = np.array(y)
	signals = [False] * len(y) # blank array of bools for mask of dataframe
	avgFilter =[0]*len(y) #moving average over lag
	stdFilter = [0]*len(y) # std deviation over lag
	avgFilter[lag - 1] = np.mean(y[0:lag])#calculate avg at start point
	stdFilter[lag - 1] = np.std(y[0:lag])#calculate std dev at start poing

	for i in range(lag,len(y)): #checks for peaks moving forward through index
		if (y[i] - avgFilter[i-1] > threshold * stdFilter[i-1]):
			signals[i] = True
			signals[i-1] = True
			filterY[i] = y[i]*influence + filterY[1-i]*(1 - influence)
			avgFilter[i] = np.mean(filterY[(i-lag):i])
			stdFilter[i] = np.std(filterY[(i-lag):i])

		else:
			filterY[i] = y[i]
			avgFilter[i] = np.mean(filterY[(i-lag):i])
			stdFilter[i] = np.std(filterY[(i-lag):i])
		#print(avgFilter[i],stdFilter[i],(y[i]-avgFilter[i-1]))'''

	for i in range((len(y)-rear_lag),-1,-1): #checks for peaks moving backwards through index

		if (y[i] - avgFilter[i+1] > rear_thresh * stdFilter[i+1]):

			signals[i] = True
			filterY[i] = y[i]*influence + filterY[i+1]*(1 - influence)
			avgFilter[i] = np.mean(filterY[i:(i+rear_lag)])
			stdFilter[i] = np.std(filterY[i:(i+rear_lag)])

		else:
			filterY[i] = y[i]
			avgFilter[i] = np.mean(filterY[i:(i+rear_lag)])
			stdFilter[i] = np.std(filterY[i:(i+rear_lag)])

	return signals

# attempt at editing baseline function that would signal end of peak once signal returned to pre-peak baseline
# It doesn't work because total baseline noise increase after peaks
def baseline2(df, lag = 10, threshold = 2, rear_thresh=3, influence = 0.5):
	y = np.array(df['tot_I']) #y is numpy array of intensities
	filterY = np.array(y)
	signals = [False] * len(y) # blank array of bools for mask of dataframe
	avgFilter =[0]*len(y) #moving average over lag
	stdFilter = [0]*len(y) # std deviation over lag
	avgFilter[lag - 1] = np.mean(y[0:lag])
	stdFilter[lag - 1] = np.std(y[0:lag])
	avgFilter[len(y)-lag] = np.mean(y[(len(y)-lag):len(y)])
	stdFilter[len(y)-lag] = np.std(y[(len(y)-lag):len(y)])
	peak_detect = False
	peak_end = 0

	for i in range(lag,len(y)):
		delta = y[i] - y[i-1]

		if (peak_detect == True) and (delta < 0):
			signals[i] = True
			print(delta)
			peak_end = i
			

		elif (peak_detect == True) and (delta > 0):
			signals[i] = False
			print(delta)
			peak_detect = False
			avgFilter[i] = np.mean(filterY[peak_end:(peak_end+lag)])
			stdFilter[i] = np.std(filterY[peak_end:(peak_end+lag)])

		elif (y[i] - avgFilter[i-1] > threshold * stdFilter[i-1]):
			signals[i] = True
			signals[i-1] = True
			filterY[i] = y[i]*influence + filterY[1-i]*(1 - influence)
			avgFilter[i] = np.mean(filterY[(i-lag):i])
			stdFilter[i] = np.std(filterY[(i-lag):i])
			peak_detect = True
			peak_start = i

		else:
			filterY[i] = y[i]
			signals[i] = False
			avgFilter[i] = np.mean(filterY[(i-lag):i])
			stdFilter[i] = np.std(filterY[(i-lag):i])
		#print(avgFilter[i],stdFilter[i],(y[i]-avgFilter[i-1]))'''

	for x in range(lag, len(y)-lag):
		''' this was an attempt to find peaks based on the slope of a moving baseline,
		with the assumption that a negative slope would indicate the start of the peak

		Once I learn how git works I'll use actual version control instead of just commenting out 
		big chunks of code

		calculate slope of line over rear_lag
		y_coords = np.arange(i,(i+rear_lag))
		x_coords = y[i:(i+rear_lag)]
		(slope, intercept) = np.polyfit(x_coords,y_coords,1)

		
		if slope < -0.001:
			print("It's a hit!")
			print(slope)'''
		i = len(y)-x
		if (y[i] - avgFilter[i+1] > rear_thresh * stdFilter[i+1]):
			signals[i] = True
			filterY[i] = y[i]*influence + filterY[i+1]*(1 - influence)
			avgFilter[i] = np.mean(filterY[i:(i+lag)])
			stdFilter[i] = np.std(filterY[i:(i+lag)])

		else:
			filterY[i] = y[i]
			signals[i] = False
			avgFilter[i] = np.mean(filterY[i:(i+lag)])
			stdFilter[i] = np.std(filterY[i:(i+lag)])

	return signals

#takes criteria from baseline and returns bool array that hits on peaks > peak_width_min
def data_clean(mask, peak_width_min=3):
	clean_mask = [False] * len(mask)
	consec = 0 #stores consecutive hits in mask
	for i in range(1,len(mask)):
		try:#This try block catches errors thrown if the last signal is marked as a peak. If this happens, mask[i+1] throws IndexError
			if (mask[i] == True) and (mask[i+1] == True):
				clean_mask[i] = True
				consec += 1
				#print("A")
			elif (consec >= 1) and (consec <= peak_width_min):#eliminates peaks less than peak width min
				for x in range((i-consec),i+1):
					clean_mask[x] = False
				#print("b")
				consec = 0
			elif (clean_mask[i-1] == True) and (mask[i] == True): #prevents erasure of last signal in peak
				clean_mask[i] = True
				consec = 0
				#print("C")
			else:
				clean_mask[i] = False
				#print("D")
		except IndexError:
			print("peak detected at end")
			break
	return clean_mask

#splits peaks based on bool mask from data_clean
#there is probably a faster/cleaner way to do this in pandas, but I can't figure it out
def peak_split(msk,df):
	mask = np.array(msk)#array from bool mask
	peaks = {}
	#indices function taken from https://bit.ly/2LuWqi0
	indices = np.nonzero(mask[1:] != mask[:-1])[0]+1 #gives index values where bool mask switches from True to False or vice versa
	y = 0
	for i,x in enumerate(indices):
		peaks[i] = df[y:x]
		y = x
	peaks[len(indices)] = df[y:] #handles segment after the last peak.
	return peaks #peak[0,2,4,ect] are noise, odd numbered peaks are actual peaks

#integrates peak. Because t differential is identical, I'm just ignoring the x component 
#find total area underneath peak, then subtracts area between 0 and approximate baseline
def integrate_peak(df):
	intens= np.array(df['tot_I'])
	last = len(intens) - 1
	baseline_trap = len(intens) * ((intens[0] + intens[last])/2) #total area below peak. (0 to peak start) 
	area = 0
	for i in range(0,last): #find total area of of peak, from 0 to peak height
		area += ((intens[i] + intens[i+1])/2)
	return area-baseline_trap

#takes peaks and return ratio of second peak/first peak
def area_ratio(peaks):
	return integrate_peak(peaks[3])/integrate_peak(peaks[1])

#
#takes coefficients from np.polyfit, and the x and y coordinates fed to polyfit. Calculates r squared
def r_squared(coeffs,x,y):
	results = 0
	p = np.poly1d(coeffs)
	# fit values, and mean
	yhat = p(x)
	ybar = np.sum(y)/len(y)
	ssreg = np.sum((yhat-ybar)**2)#regression sum of squares
	sstot = np.sum((y - ybar)**2)#Total Sum of Squares
	results = ssreg / sstot
	return results

#Takes coefficients from repsonse curve, plus a chromatogram. Returns the area ratio (Nic/Quin)
def unknown_conc(coeffs, unk):
	unk_peaks = {}
	area_ratio_unk = 0
	unk_peaks = peak_split(data_clean(baseline(unk)), unk)
	area_ratio_unk = area_ratio(unk_peaks)
	curve = np.poly1d(coeffs)
	return  curve(area_ratio_unk)

#takes a list of peaks from peak_split and a retention time (rt), 
#returns the peak that contains the specificed retention time
def peak_select(peaks, rt):
	for i in range(len(peaks)):
		if rt in peaks[i].RT.round(decimals=2).values:
			return peaks[i]

#finds the width of a peak at half its height
def half_height(peak,chrom):
	left_points = np.zeros((2,2), dtype=float)
	right_points = np.zeros((2,2), dtype=float)
	baseline_points = np.zeros((2,2), dtype=float)
	baseline_points[0] = peak.iloc[0].values#RT and Intesity and for first point in peak
	baseline_points[1] = peak.iloc[len(peak.index)-1].values#RT and Intesity and for last point in peak
	baseline = np.polyfit(baseline_points[:,0],baseline_points[:,1],1)
	base_curve = np.poly1d(baseline)
	height_tot = peak['tot_I'].idxmax(axis = 0)#index for the highest signal in the peak
	peak_max = peak.loc[height_tot].values#RT and Intesity for highest point
	height = peak_max[1] - base_curve(peak_max[0])#Height = intensity at height minus calculated baseline
	half_height = height/2
	#splits peak at peak. left curve is before peak, right_curve is points after peak
	left_curve = peak.loc[0:height_tot]
	right_curve = peak.loc[height_tot:]
	#splits left_curve into peaks less than half height and greater than half height
	left_points_less_half = left_curve[(left_curve['tot_I'] < half_height)]
	left_points_greater_half = left_curve[(left_curve['tot_I'] > half_height)]
	#splits right_curve into peaks less than half height and greater than half height
	right_points_less_half = right_curve[(right_curve['tot_I'] < half_height)]
	right_points_greater_half = right_curve[(right_curve['tot_I'] > half_height)]
	#assigns points(RT,tot_I) to array that will be used to calculate eqn of line
	left_points[0] = left_points_less_half.iloc[-1].values
	left_points[1] = left_points_greater_half.iloc[0].values
	right_points[0] = right_points_less_half.iloc[0].values
	right_points[1] = right_points_greater_half.iloc[-1].values
	left_line = np.polyfit(left_points[:,1],left_points[:,0],1)
	left_curve = np.poly1d(left_line)
	right_line = np.polyfit(right_points[:,1],right_points[:,0],1)
	right_curve = np.poly1d(right_line)
	half_height_width = right_curve(half_height) - left_curve(half_height)
	return half_height_width


Nic_exp_RT = 2.49
Quinoline_exp_RT = 1.89	
t_diff_min = 0.0102#0.612 seconds		

C = chrom_import('Unknowns/C.csv')
coeffs = np.array([1.51147293, 0.86384225])
unknown_conc(coeffs, C)
nicotine_peak_df = peak_select(peak_split(data_clean(baseline(C)), C), Quinoline_exp_RT)
half_height(nicotine_peak_df,C)

'''FILE_PATHS = ['standards_data/Nicotine_1ppm_180712201814.csv']
FILE_PATHS += ['standards_data/Nicotine_2,5ppm_180712202451.csv']
FILE_PATHS += ['standards_data/Nicotine_5ppm_180712203135.csv']
FILE_PATHS += ['standards_data/Nicotine_7,5ppm_180712203822.csv']
FILE_PATHS += ['standards_data/Nicotine_10ppm_180712204456.csv']
area_ratios = [0,0,0,0,0]
chroms = {}
peaks = {}

for i in range(len(FILE_PATHS)):
	chroms[i] = chrom_import(FILE_PATHS[i])
	peaks[i] = peak_split(data_clean(baseline(chroms[i])),chroms[i])
	ax = chroms[i].plot(x ='RT',y='tot_I')
	print(area_ratio(peaks[i]))
	area_ratios[i] = area_ratio(peaks[i])
	peaks[i][1].plot(ax=ax, x ='RT',y='tot_I', style='ro')
	peaks[i][3].plot(ax=ax, x ='RT',y='tot_I', style='ro')
	plt.show()

x = area_ratios
y = [0.9904,2.476,4.952,7.428,9.904]
(a, b) = np.polyfit(x, y, 1)
coeffs = np.polyfit(x,y, 1)
r2 = r_squared(coeffs, x, y)
print('y = {:0.05f} x + {:0.05f} : R Squared = {:0.05f}'.format(a, b, r2))

fig1, axes1 = plt.subplots(1,1,sharex=True)
axes1.scatter(x, y, marker='s')
axes1.plot(np.unique(x), np.poly1d((a,b))(np.unique(x)), linewidth=2.0, linestyle='--', color='r')
axes1.set(ylabel='ppm', xlabel='Area Ratio')
fig1.suptitle('Nicotine\ny = {:0.05f} x + {:0.05f} -- R Squared = {:0.05f}'.format(a, b, r2))
plt.show()

UNKNOWNS = ['Unknowns/A.csv']
UNKNOWNS += ['Unknowns/B.csv']
UNKNOWNS += ['Unknowns/C.csv']
UNKNOWNS += ['Unknowns/D.csv']

for file in UNKNOWNS:
	unknown = chrom_import(file)
	unknown_conc(coeffs, unknown)'''

