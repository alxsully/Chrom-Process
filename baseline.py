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

#combines above three modules
def chrom_import(file):
	nic = sum_I(rename_RT(import_SIM_chrom(file)))	
	return nic

#Peak Find Module
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
def data_clean(mask, peak_width_min=5):
	clean_mask = [False] * len(mask)
	consec = 0 #stores consecutive hits in mask
	for i in range(1,len(mask)):
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

#I just stole this from stack overflow.
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

#def unknow_conc(unk):




FILE_PATHS = ['standards_data/Nicotine_1ppm_180712201814.csv']
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
	#ax = chroms[i].plot(x ='RT',y='tot_I')
	print(area_ratio(peaks[i]))
	area_ratios[i] = area_ratio(peaks[i])
	#peaks[i][1].plot(ax=ax, x ='RT',y='tot_I', style='ro')
	#peaks[i][3].plot(ax=ax, x ='RT',y='tot_I', style='ro')
	#plt.show()

y = [0.9904,2.476,4.952,7.428,9.904]
x = area_ratios
(a, b) = np.polyfit(x, y, 1)
coeffs = np.polyfit(x,y, 1)
r2 = r_squared(coeffs, x, y)
print('y = {:0.05f} x + {:0.05f}'.format(a, b))

fig1, axes1 = plt.subplots(1,1,sharex=True)
axes1.scatter(x, y, marker='s')
axes1.plot(np.unique(x), np.poly1d((a,b))(np.unique(x)), linewidth=2.0, linestyle='--', color='r')
axes1.set(ylabel='ppm', xlabel='Area Ratio')
fig1.suptitle('Nicotine\ny = {:0.05f} x + {:0.05f} -- R Squared = {:0.05f}'.format(a, b, r2))
plt.show()


'''nic = chrom_import("Nicotine_5ppm.csv")
crit = baseline(nic)
clean_crit = data_clean(crit)
peaks = peak_split(clean_crit, nic)
ax = nic.plot(x ='RT',y='tot_I')
peaks[1].plot(ax=ax, x ='RT',y='tot_I', style='ro')
peaks[3].plot(ax=ax, x ='RT',y='tot_I', style='ro')
plt.show()'''

