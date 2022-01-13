

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

		self.path  = Path(path_data)
		self.reset = reset

		self.pkl_nam = f"smps2aps_{start_time.strftime('%Y%m%d')}-{end_time.strftime('%Y%m%d')}.pkl"
		
		print(f" from {start_time.strftime('%Y-%m-%d %X')} to {end_time.strftime('%Y-%m-%d %X')}")
		print('='*65)
		print(f"{dtm.now().strftime('%m/%d %X')}")


	def __read_data(self,):
		## read raw data
		with open(self.path/'smps_aps_raw.pkl','rb') as f:
			_dt = pkl.load(f)
		_aps, _smps = _dt['aps'], _dt['smps']

		return _aps['data'], _aps['total'], _smps['data'], _smps['total']

	def __test_read_data(self,):
		## read raw data
		with open(self.path/'smps_aps_raw.pkl','rb') as f:
			_dt = pkl.load(f)
		_aps, _smps = _dt['aps'].loc[dtm(2021,11,8,18):dtm(2021,11,8,18,54)], _dt['smps'].loc[dtm(2021,11,8,18):dtm(2021,11,8,18,54)]

		return _aps['data'], _aps['total'], _smps['data'], _smps['total']

	def __pre_process(self,_aps,_aps_t,_smps,_smps_t,_aps_hb,_smps_lb):
		## quality control

		## discard missing data(data equal to 0)
		## bins larger than aps smallest bin
		## bins smaller than aps fitting bin limit
		_smps_overlap = _smps[_smps.keys()[_smps.keys()>_smps_lb]].copy()
		_smps.loc[_smps_overlap.mean(axis=1)<=0,_smps_overlap.keys()] = n.nan

		_aps_overlap = _aps[_aps.keys()[_aps.keys()<_aps_hb]].copy()
		_aps.loc[_aps_overlap.mean(axis=1)<=0,_aps_overlap.keys()] = n.nan

		## discard outlier by 1-hr total conc. profile
		## data over 1.5 std will set to be nan
		def _quality_ctrl(_data,_total):
			def _filter(_df):
				_df = _df.total
				_std, _mean = _df.std().copy(), _df.mean().copy()
				_df.loc[(_df>_mean+1.5*_std)|(_df<_mean-1.5*_std)] = n.nan
				return _df

			_total['grp_time'] = _total.index.strftime('%Y%m%d %H00')

			_total = _total.groupby('grp_time',group_keys=False).apply(_filter).copy()
			_total_condi = _total[_total.isna()].index

			_data.loc[_total_condi] = n.nan

			return _data

		return  _quality_ctrl(_aps,_aps_t), _quality_ctrl(_smps,_smps_t)


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





	## aps_fit_highbound : the diameter I choose randomly
	def merge_data(self,ave_time='1h',aps_fit_highbound=1382.,smps_overlap_lowbound=523.):

		## read raw data
		# aps, aps_total, smps, smps_total = self.__read_data()
		aps, aps_total, smps, smps_total = self.__test_read_data()

		## pre-process data
		aps, smps = self.__pre_process(aps,aps_total,smps,smps_total,aps_fit_highbound,smps_overlap_lowbound)

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
		



