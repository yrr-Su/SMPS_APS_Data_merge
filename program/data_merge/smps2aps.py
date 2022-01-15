

import numpy as n
from os import listdir, mkdir
from os.path import join as pth, exists, dirname, realpath
from pandas import date_range, concat, DataFrame
from datetime import datetime as dtm
from datetime import timedelta as dtmdt
from pathlib import PurePath as Path
import pickle as pkl
np = n

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
class merger:

	## initial setting
	## input : file path, 
	## 		   start time,
	## 		   final time,
	## 		   reset switch
	## 
	## because the pickle file will be generated after read raw data first time,
	## if want to re-read the rawdata, please set 'reset=True'

	def __init__(self,path_data,start_time,end_time,reset=False,):

		print(f'\nSMPS and APS data merge')
		print('='*65)
		print(f"Merge smps and aps, from mobility diameter to aerodynamic diameter")

		## class parameter
		self.index 	  = date_range(start_time,end_time,freq='6T')
		self.__time   = (start_time,end_time)

		self.path  = Path(path_data)
		self.reset = reset

		self.pkl_nam = f"smps2aps_{start_time.strftime('%Y%m%d')}-{end_time.strftime('%Y%m%d')}.pkl"
		
		print(f" from {start_time.strftime('%Y-%m-%d %X')} to {end_time.strftime('%Y-%m-%d %X')}")
		print('='*65)
		print(f"{dtm.now().strftime('%m/%d %X')}")

	## Read raw data(has been processed by rawdata_process)
	## return : aps, 
	## 			total conc. of aps, 
	## 			smps, 
	## 			total conc. smps
	def __read_data(self,):
		## read raw data
		with open(self.path/'smps_aps_raw.pkl','rb') as f:
			_dt = pkl.load(f)
		_aps, _smps = _dt['aps'].reindex(self.index), _dt['smps'].reindex(self.index),

		return _aps[_aps.keys()[:-1]], _aps['total'].to_frame(), _smps[_smps.keys()[:-1]], _smps['total'].to_frame()

	## Quality control
	## return : aps after QC,
	## 			smps after QC,
	def __pre_process(self,_aps,_aps_t,_smps,_smps_t,_aps_hb,_smps_lb):

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

	## Overlap fitting
	## Create a fitting func. by smps data
	## return : shift factor
	def __overlap_fitting(self,_aps,_smps,_aps_hb,_smps_lb):

		## overlap fitting
		## parmeter
		_dt_size = len(_aps)

		## overlap diameter data
		_aps  = _aps[_aps.keys()[_aps.keys()<_aps_hb]].copy()
		_smps = _smps[_smps.keys()[_smps.keys()>_smps_lb]].copy()

		## power law fit to SMPS num conc at upper bins to log curve
		## y = Ax^B, A = e**coefa, B = coefb, x = logx, y = logy
		## ref : http://mathworld.wolfram.com/LeastSquaresFittingPowerLaw.html

		## coefficient A, B
		_aps_qc_cond = ((_aps!=0)&n.isfinite(_aps))
		_aps_qc = _aps.where(_aps_qc_cond)

		_size = _aps_qc_cond.sum(axis=1)
		_size = _size.where(_size!=0.).copy()

		_logx, _logy = n.log(_aps_qc.keys()._data.astype(float)), n.log(_aps_qc)
		_x, _y, _xy, _xx = _logx.sum(), _logy.sum(axis=1), (_logx*_logy).sum(axis=1), (_logx**2).sum()

		_coeB = ((_size*_xy-_x*_y)/(_size*_xx-_x**2.))
		_coeA = n.exp((_y-_coeB*_x)/_size).values.reshape(_dt_size,-1)
		_coeB = _coeB.values.reshape(_dt_size,-1)

		## rebuild shift smps data by coe. A, B
		## x_shift = (y_ori/A)**(1/B)
		## y_shift = A*x_shift**B
		## y_shift_min = sum((y_shift-ybase)**2).argmin() (calculate by each column)
		_smps_shift_x = (_smps/_coeA)**(1/_coeB)
		_smps_shift_y = (_smps_shift_x*_coeA)**_coeB

		_temp = DataFrame()
		_smps_base = n.log(_smps)
		for _idx, (_key, _df) in enumerate(_smps_shift_y.iteritems()):
			_df = _df.values.reshape(_dt_size,-1)
			_temp[_idx] = ((_df-_smps_base)**2).sum(axis=1)
		_smps_shift_ymin = _temp.idxmin(axis=1)

		## find the shift factor which contribute miniumum y shift
		## x_shift_factor = x_ori/x_shift
		_smps_shift_factor = (_smps_shift_x.keys()._data.astype(float)/_smps_shift_x).values

		return _smps_shift_factor[range(_dt_size),_smps_shift_ymin.values].copy()

	
	## Remove big shift data
	## Average datga in average time
	## Return : aps, smps (mean)
	def __shift_data_process(self,_aps,_smps,_shift,_ave):

		_big_shift = (_shift>10.)|(_shift<-2.)
		_smps.loc[_big_shift] = n.nan
		_aps.loc[_big_shift]  = n.nan


		return _aps.resample('1h').mean(), _smps.resample('1h').mean()


	## Create merge data
	## Return : merge data, (density, not yet, fk out)
	def __merge_data(self,):



		_smps_bins = n.full(smps.shape,smps.keys()._data)

		_aps_bins = aps.keys()[aps.keys()>aps[_aps_fit_region].keys()[-1]]._data
		_aps_bins_shift = n.full((aps.index.size,_aps_bins.size),_aps_bins)/shift_x.reshape(-1,1)

		merge_conc = []	
		def merge_bin_conc(_smps_bin,_aps_bin,_smps,_aps):

			if n.isnan(_smps).all():
				merge_conc.append(n.nan)
				return n.nan

			_condi = _smps_bin<=_aps_bin[0]

			merge_conc.append(n.nan if (~_condi.any()) else n.hstack((_smps[n.where(_condi)[0]],_aps)))
			return n.nan if (~_condi.any()) else n.hstack((_smps_bin[_condi],_aps_bin))

		merge_bin = list(map(merge_bin_conc,_smps_bins,_aps_bins_shift,smps.values,aps[_aps_bins].values))





	## aps_fit_highbound : the diameter I choose randomly
	def merge_data(self,ave_time='1h',aps_fit_highbound=1382.,smps_overlap_lowbound=523.):

		## read raw data
		aps, aps_total, smps, smps_total = self.__read_data()

		## pre-process data
		aps, smps = self.__pre_process(aps,aps_total,smps,smps_total,aps_fit_highbound,smps_overlap_lowbound)

		## shift infomation, calculate by powerlaw fitting
		shift = self.__overlap_fitting(aps,smps,aps_fit_highbound,smps_overlap_lowbound)

		## process data by shift infomation, and average data
		aps, smps = self.__shift_data_process(aps,smps,shift,ave_time)

		## merge aps and smps
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
		



