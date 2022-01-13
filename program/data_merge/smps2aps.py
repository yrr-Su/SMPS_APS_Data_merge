

import numpy as n
from os import listdir, mkdir
from os.path import join as pth, exists, dirname, realpath
from pandas import date_range, concat, DataFrame
from datetime import datetime as dtm
from datetime import timedelta as dtmdt
import pickle as pkl


# bugs box
"""




# """

__all__ = [
		'merge_data',
		'save_data',
		'get_data',
	]

# read pkl file
# merge smps and aps
class merge_func:
	
	nam = None

	## initial setting
	## input : file path, 
	## 		   start time,
	## 		   final time,
	## 		   reset switch
	## 
	## because the pickle file will be generated after read raw data first time,
	## if want to re-read the rawdata, please set 'reset=True'

	def __init__(self,path_data,start_time,end_time,reset=False,):

		print(f'\n{self.nam}')
		print('='*65)
		print(f"Merge smps and aps, from mobility diameter to aerodynamic diameter")

		## class parameter
		self.index 	  = lambda _freq: date_range(start_time,end_time,freq=_freq)
		self.__time   = (start_time,end_time)

		self.path  = path_data
		self.reset = reset

		self.pkl_nam = f"smps2aps_{start_time.strftime('%Y%m%d')}-{end_time.strftime('%Y%m%d')}.pkl"
		
		print(f" from {start_time.strftime('%Y-%m-%d %X')} to {end_time.strftime('%Y-%m-%d %X')}")
		print('='*65)
		print(f"{dtm.now().strftime('%m/%d %X')}")




	def __pre_process(self,):
		pass
		## discard missing data 


		## quality control
		## SMPS overlap region remove the no value profile
		## if overlap region have no value(average value equal to 0), set as nan
		## overlap region is from 0.523um/523nm





	def __overlap_fitting(self,):
		pass
		## overlap fitting
		## lowest APS bins in overlap region (starting after bin 2)
		## return shift infomation




		## power law fit to SMPS num conc at upper bins to log curve
		## y = Ax^B, A = e**coefa, B=coefb, x = logx, y = logy
		## ref : http://mathworld.wolfram.com/LeastSquaresFittingPowerLaw.html



		## coefficient function




	def __shift_data_process(self,):
		pass
		## return data deal with shift infomation




	
	def __merge_data(self,):
		pass
		## return 






	def merge_data(self,ave_time='1h',smps_fit_lowbound=340.,smps_overlap_lowbound=523.):
		pass



		## pre-process data

		self.__pre_process()

		## shift infomation, calculate by powerlaw fitting
		shift = self.__overlap_fitting()
		

		## process data by shift infomation
		smps, aps = self.__shift_data_process(shift)


		self.fout = fout

		return fout



	def save_data(self,path_save=None):
		pass
		path_save = path_save if path_save is not None else self.path




	def get_data(self,**kwarg):
		pass
		
		default_par = {start_time : self.start_time,
					   end_time   : self.end_time,
					   path_data  : self.path,
					   ave_time   : '1h',
					}
		default_par.update(kwarg)
		



