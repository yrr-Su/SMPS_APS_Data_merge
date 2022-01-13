

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

	def __init__(self,_path,_sta,_fin,reset=False,):

		print(f'\n{self.nam}')
		print('='*65)
		print(f"Merge smps and aps, from mobility diameter to aerodynamic diameter")

		## class parameter
		self.index 	  = lambda _freq: date_range(_sta,_fin,freq=_freq)
		self.__time   = (_sta,_fin)

		self.path  = _path
		self.reset = reset

		self.pkl_nam = f"smps2aps_{_sta.strftime('%Y%m%d')}-{_fin.strftime('%Y%m%d')}.pkl"
		
		print(f" from {_sta.strftime('%Y-%m-%d %X')} to {_fin.strftime('%Y-%m-%d %X')}")
		print('='*65)
		print(f"{dtm.now().strftime('%m/%d %X')}")




	def __pre_process():
		
		## discard missing data 


		## quality control
		## SMPS overlap region remove the no value profile
		## if overlap region have no value(average value equal to 0), set as nan
		## overlap region is from 0.523um/523nm





	def __overlap_fitting():
		
		## overlap fitting
		## lowest APS bins in overlap region (starting after bin 2)
		## return shift infomation




		## power law fit to SMPS num conc at upper bins to log curve
		## y = Ax^B, A = e**coefa, B=coefb, x = logx, y = logy
		## ref : http://mathworld.wolfram.com/LeastSquaresFittingPowerLaw.html



		## coefficient function
















	def merge_data(ave_time='1h',smps_fit_lowbound=340.,smps_overlap_lowbound=523.):
		



		## pre-process data

		self.__pre_process()

		## shift infomation, calculate by powerlaw fitting
