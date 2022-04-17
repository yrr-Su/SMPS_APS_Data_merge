

import numpy as n
from os import listdir, mkdir
from os.path import join as pth, exists, dirname, realpath
from pandas import date_range, concat, DataFrame
from datetime import datetime as dtm
from datetime import timedelta as dtmdt
from pathlib import Path
import pickle as pkl
from scipy.interpolate import interp1d
np = n

# bugs box
"""
raw_index : default 6T, but it is not suit for process data



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

	def __init__(self,path_data,reset=False,
				 input_QCdata=False,QCdata_freq=None,**kwarg):
		print(f'\nSMPS and APS data merge')
		print('='*65)
		print(f"Merge smps and aps, from mobility diameter to aerodynamic diameter")

		## class parameter
		self.path  = Path(path_data)
		self.reset = reset
		self.input_QCdata = input_QCdata
		self.data_freq  = QCdata_freq or '1h'
	
		print('='*65)
		print(f"{dtm.now().strftime('%m/%d %X')}")

	## Read raw data(has been processed by rawdata_process)
	## return : aps, 
	## 			total conc. of aps, 
	## 			smps, 
	## 			total conc. smps
	def __read_data(self,):
		print(f"\t{dtm.now().strftime('%m/%d %X')} : \033[96mreading raw data\033[0m")

		## read raw data
		with open(self.path/'smps_aps_raw.pkl','rb') as f:
			_dt = pkl.load(f)
			_aps, _smps = _dt['aps'], _dt['smps']

			## add env parameter
			self.out_index = _aps.index.copy()
			_st, _ed = self.out_index[[0,-1]]

			self.csv_nam = lambda _ : f"smps2aps_{_}_{_st.strftime('%Y%m%d%H00')}-{_st.strftime('%Y%m%d%H00')}.csv"
			self.pkl_nam = f"smps2aps_{_ed.strftime('%Y%m%d%H00')}-{_ed.strftime('%Y%m%d%H00')}.pkl"
	
		if self.input_QCdata: 
			return _aps, None, _smps, None

		return _aps[_aps.keys()[:-1]], _aps['total'].to_frame(), _smps[_smps.keys()[:-1]], _smps['total'].to_frame()

	## Quality control
	## return : aps after QC(1-hr ave.),
	## 			smps after QC(1-hr ave.),
	def __pre_process(self,_aps,_aps_t,_smps,_smps_t,_aps_hb,_smps_lb):
		
		## processed data, has been QC
		if self.input_QCdata: 
			return _aps, _smps

		print(f"\t{dtm.now().strftime('%m/%d %X')} : \033[96mpre-process data\033[0m")

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
				if (_df.isna().sum()) != (_df.__len__()):
					_std, _mean = _df.std(), _df.mean()
					_df.loc[(_df>_mean+1.5*_std)|(_df<_mean-1.5*_std)] = n.nan
				return _df

			_total['grp_time'] = _total.index.strftime('%Y%m%d %H00')

			_total = _total.groupby('grp_time',group_keys=False).apply(_filter).copy()
			_total_condi = _total[_total.isna()].index

			_data.loc[_total_condi] = n.nan

			return _data.resample(self.data_freq).mean()
		return  _quality_ctrl(_aps,_aps_t), _quality_ctrl(_smps,_smps_t)

	## Overlap fitting 
	## Create a fitting func. by smps data
	## return : shift factor
	def __overlap_fitting(self,_aps,_smps,_aps_hb,_smps_lb):
		print(f"\t{dtm.now().strftime('%m/%d %X')} : \033[96moverlap range fitting\033[0m")

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
		_coeA = n.exp((_y-_coeB*_x)/_size).values.reshape(-1,1)
		_coeB = _coeB.values.reshape(-1,1)

		## rebuild shift smps data by coe. A, B
		## x_shift = (y_ori/A)**(1/B)
		## y_shift = A*x_shift**B
		## y_shift_min = sum((y_shift-ybase)**2).argmin() (calculate by each column)
		_smps_shift_x = (_smps/_coeA)**(1/_coeB)
		_smps_shift_x = _smps_shift_x.where(~n.isinf(_smps_shift_x))
		_smps_shift_y = (_smps_shift_x*_coeA)**_coeB

		_temp = DataFrame()
		_smps_base = n.log(_smps)
		for _idx, (_key, _df) in enumerate(_smps_shift_y.iteritems()):
			_df = _df.values.reshape(-1,1)
			_temp[_idx] = ((_df-_smps_base)**2).sum(axis=1)
		_smps_shift_ymin = _temp.idxmin(axis=1)

		## find the shift factor which contribute miniumum y shift
		## x_shift_factor = x_aps/x_smps
		_smps_shift_factor = (_smps_shift_x.keys()._data.astype(float)/_smps_shift_x).values

		return _smps_shift_factor[range(_dt_size),_smps_shift_ymin.values].copy().reshape(-1,1)

	
	## Remove big shift data
	## Return : aps, smps, shift (without big shift data)
	def __shift_data_process(self,_aps,_smps,_shift):
		print(f"\t{dtm.now().strftime('%m/%d %X')} : \033[96mshift-data quality control\033[0m")
		_big_shift = (_shift>10.)|(_shift<-2.)|(n.isnan(_shift))
		return _aps.loc[~_big_shift], _smps.loc[~_big_shift], _shift[~_big_shift].reshape(-1,1)


	## Create merge data
	##  shift all smps bin and remove the aps bin which smaller than the latest old smps bin
	## Return : merge bins, merge data, density
	def __merge_data(self,_aps,_smps,_shift):
		print(f"\t{dtm.now().strftime('%m/%d %X')} : \033[96mcreate merge data\033[0m")

		## get merge data
		_smps_bin	  = _smps.keys()._data

		_aps_bin  	  = _aps.keys()._data.astype(float)
		_aps_bin	  = _aps_bin[_aps_bin>_smps_bin[-1]].copy()

		_aps		  = _aps[_aps_bin].copy()
		_aps_bin_ary  = n.full(_aps.shape,_aps_bin)

		## merge
		##  un-equal length data
		## 	loop in each time
		_smps_bin_shift = n.full(_smps.shape,_smps_bin)*_shift ## different to origin algorithm
		_max_bin_num = (_aps_bin_ary>_smps_bin_shift[:,-1].reshape(-1,1)).sum(axis=1).max()
		
		_bins_lst, _data_lst, _smps_inte_lst = [], [], []
		_bin_aps = _aps_bin_ary[0]
		for _bin_smps, _dt_aps, _dt_smps in zip(_smps_bin_shift,_aps.values,_smps.values):

			## shift data and shift bin
			_condi = _bin_aps>=_bin_smps[-1]
			_append_ary = n.full(_max_bin_num-_condi.sum(),n.nan)

			_bins_lst.append(n.hstack((_bin_smps,_bin_aps[_condi],_append_ary)))
			_data_lst.append(n.hstack((_dt_smps,_dt_aps[_condi],_append_ary)))

			## linear interpolate and extrapolate for range of smps
			_inte_fc = interp1d(_bin_smps,_dt_smps,kind='linear',fill_value='extrapolate')
			_smps_inte_lst.append(_inte_fc(_smps_bin))

		_data_inte = concat([DataFrame(_smps_inte_lst,columns=_smps_bin,index=_aps.index),_aps],axis=1)

		## shift factor to rho
		_rho_list = (_shift**2).flatten()

		## process output df
		## average, align with index
		def _out_df(*_df_arg,**_df_kwarg):
			_df = DataFrame(*_df_arg,**_df_kwarg).set_index(_aps.index).reindex(self.out_index)
			_df.index.name = 'time'
			return _df
	
		return _out_df(_bins_lst), _out_df(_data_lst), _out_df(_rho_list), _out_df(_data_inte)


	## aps_fit_highbound : the diameter I choose randomly
	def merge_data(self,aps_fit_highbound=1382,smps_overlap_lowbound=523):
		print(f'\nMerge data :')
		print(f' APS fittint higher diameter : {aps_fit_highbound:4d} nm')
		print(f' SMPS overlap lower diameter : {smps_overlap_lowbound:4d} nm')
		print(f' Average time                : {self.data_freq:>4s}\n')

		self.fout = None

		## read raw data
		raw_data = self.__read_data()

		## pre-process data
		aps, smps = self.__pre_process(*raw_data,aps_fit_highbound,smps_overlap_lowbound)

		## shift infomation, calculate by powerlaw fitting
		shift = self.__overlap_fitting(aps,smps,aps_fit_highbound,smps_overlap_lowbound)

		## process data by shift infomation, and average data
		aps, smps, shift = self.__shift_data_process(aps,smps,shift)

		## merge aps and smps
		bins, data, density, comp_data = self.__merge_data(aps,smps,shift)
		
		self.fout = {'bins'    : bins,
					 'data'    : data,
					 'density' : density,
					 'all'     : comp_data,
					}

	## save data to given path
	def save_data(self,path_save=None):
		path_save = Path(path_save) if path_save is not None else self.path
		
		print(f'\nSave data :')
		print(f' Save path : {realpath(path_save)}')

		if self.fout is not None:
			## save to pickle
			with open(path_save/self.pkl_nam,'wb') as f:
				pkl.dump(self.fout,f,protocol=pkl.HIGHEST_PROTOCOL)
			
			## save to csv
			for nam, df in self.fout.items():
				df.to_csv(path_save/self.csv_nam(nam))
		else:
			print('Please Use Function(merge_data) First !!!')



	## maybe not necessary
	def get_data(self,**kwarg):
		pass
		
		default_par = {start_time : self.start_time,
					   end_time   : self.end_time,
					   path_data  : self.path,
					   ave_time   : '1h',
					}
		default_par.update(kwarg)
		



