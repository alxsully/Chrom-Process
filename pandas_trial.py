import pandas as pd
import numpy as np
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt
import scipy as sp

#imports chromatogram.csv generated using arbitrary xcalibur method, returns only columns for
#selected masses(default is SIM_Nicotine selected masses)
def import_SIM_chrom(file, m_z = ['76','84','102','129','133','162']):
	rt_label = ['RT(minutes) - NOT USED BY IMPORT']
	import_list = rt_label + m_z
	return pd.read_csv(file, sep=";", usecols=import_list) 

#sets RT(minutes) column to index, renames to 'RT' and drops RT(milliseconds)
def index_RT(df):
	df.set_index(
	'RT(minutes) - NOT USED BY IMPORT', inplace=True)
	df.index.names = ['RT']
	return df

#sums intensities of various masses and adds column 'tot_I'
def sum_I(df):
	df['tot_I'] = np.sum(df, axis=1)
	df.drop(columns = ['76','84','102','129','133','162'], inplace=True)
	return df

def nearest_rt(df, rt):
	idx = (np.abs(df.index.values - rt)).argmin()
	#return df.iloc[idx,6]
	return idx

def peak_max(df):
	i_array = df['tot_I'].values
	idx = np.argmax(i_array)
	print (idx)
	return df.iloc[idx,0]

def peak_max_idx(df):
	i_array = df['tot_I'].values
	idx = np.argmax(i_array)		
	return idx

def near_val(df, val):
	return (np.abs(df['tot_I'].values - val)).argmin()

def baseline(df):
	front = df.iloc[nearest_rt(df,2.4):nearest_rt(df,2.44),6].values
	back = df.iloc[nearest_rt(df,2.65):nearest_rt(df,3),6].values
	baseline = np.concatenate((front,back),axis = 0)
	base_i = np.mean(baseline)
	return base_i

def baseline_std(df):
	front = df.iloc[nearest_rt(df,2.4):nearest_rt(df,2.44),6].values
	back = df.iloc[nearest_rt(df,2.65):nearest_rt(df,3),6].values
	baseline = np.concatenate((front,back),axis = 0)
	base_std = np.std(baseline)
	return base_std

def plot_baseline(df):
	front = df.iloc[nearest_rt(df,2.4):nearest_rt(df,2.44),6].values
	front_idx = df.iloc[nearest_rt(df,2.4):nearest_rt(df,2.44),6].index
	back = df.iloc[nearest_rt(df,2.65):nearest_rt(df,3),6].values
	back_idx = df.iloc[nearest_rt(df,2.65):nearest_rt(df,3),6].index
	baseline = np.concatenate((front,back),axis = 0)
	baseline_idx = np.concatenate((front_idx,back_idx), axis = 0)
	plt.plot(baseline_idx,baseline, 'ro')


def height(df):
	return peak_max(df)-baseline(df)

def df_slice(df):
	criteria = (df.index > 2.4) & (df.index < 3)
	print(criteria)
	return df[criteria]




nic = {}
FILE_PATHS = ['Nicotine_1ppm_180712201814.csv']
FILE_PATHS += ['Nicotine_2,5ppm_180712202451.csv']
FILE_PATHS += ['Nicotine_5ppm_180712203135.csv']
FILE_PATHS += ['Nicotine_7,5ppm_180712203822.csv']
FILE_PATHS += ['Nicotine_10ppm_180712204456.csv']
for i in range(len(FILE_PATHS)):
	nic[i] = index_RT(import_SIM_chrom(FILE_PATHS[i]))	
	nic[i] = sum_I(nic[i])

nic_peak = df_slice(nic[2])
peak_max(nic[2])


'''coefs = poly.polyfit(nic_peak.index, nic_peak['tot_I'],200)
ffit = poly.polyval(nic_peak.index,coefs)
plt.plot(nic_peak.index,ffit)
plt.show()'''


nic[2].plot(y='tot_I', use_index=True)
#nic_peak.plot.scatter(nic_peak[peak_max_idx(nic_peak)].index,nic_peak[peak_max_idx(nic_peak),0],'ro')
#print nic_peak.iloc[nearest_rt(nic_peak,2.65):nearest_rt(nic_peak,2.7),0].index

plt.show()




