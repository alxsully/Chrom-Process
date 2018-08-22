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

def nearest_rt(df, rt):
	idx = (np.abs(df.index.values - rt)).argmin()
	return df.iloc[idx,1]

blank_files = []
blanks = {}
rt_I = np.zeros(20)

for i in range(1,21):
	blank_files += ['blanksdata/Nicotine_SN_Cal_Blank{:0=2}.csv'.format(i)]

for i in range(len(blank_files)):
	blanks[i] = chrom_import(blank_files[i])
	rt_I[i] = nearest_rt(blanks[i],2.49)

print(np.mean(rt_I),np.std(rt_I))
	