# initial function 
# version : 

# target : 
# 1. read data and process

from datetime import datetime as dtm
from os import listdir, mkdir
from os.path import join as pth, exists, dirname, realpath
from pathlib import PurePath as Path
import pickle as pkl
from numpy import array, nan
from pandas import date_range, concat, read_table
import json as jsn

## bugs box
"""




# """


__all__ = [
		'save_data',
	]


# parameter
cur_file_path = dirname(realpath(__file__))
with open(Path(cur_file_path)/'meta.json','r') as f:
	meta_dt = jsn.load(f)

# class
## parant class (read file)
## list the file in the path and 
## read pickle file if it exisits, else read raw data and dump the pickle file
class reader:
	
	## initial setting
	## input : file path, 
	## 		   start time,
	## 		   final time,
	## 		   reset switch
	## 
	## because the pickle file will be generated after read raw data first time,
	## if want to re-read the rawdata, please set 'reser=True'

	def __init__(self,path_raw,path_save,stara_time,end_time,reset=False):
		print(f'\nSMPS and APS')
		print('='*65)
		print(f"Reading file and process data")

		## class parameter
		self.index = lambda _freq: date_range(stara_time,end_time,freq=_freq)
		self.path  = Path(path_raw)
		self.path_out  = Path(path_save)
		self.reset = reset
		self.meta_read = meta_dt['read']
		self.out_nam = 'smps_aps_raw.pkl'
		self.__time  = (stara_time,end_time)
		
		print(f" from {stara_time.strftime('%Y-%m-%d %X')} to {end_time.strftime('%Y-%m-%d %X')}")
		print('='*65)
		print(f"{dtm.now().strftime('%m/%d %X')}")

	def __smps_reader(self,_file):
		## customize each instrument
		## read one filess
		with open(_file,'r',encoding='utf-8',errors='ignore') as f:
			_df = read_table(f,skiprows=18,parse_dates={'Time':['Date','Start Time']}).set_index('Time')
			_key = list(_df.keys()[7:-26])

			_newkey = {}
			for _k in _key: _newkey[_k] = float(_k).__round__(4)
			_newkey['Total Conc.(#/cm)'] = 'total'

		return _df[_key+['Total Conc.(#/cm)']].rename(_newkey,axis=1)

	def __aps_reader(self,_file):
		## customize each instrument
		## read one filess
		with open(_file,'r',encoding='utf-8',errors='ignore') as f:

			_df = read_table(f,skiprows=6,parse_dates={'Time':['Date','Start Time']}).set_index('Time')
			_key = list(_df.keys()[2:54])

			_newkey = {}
			for _k in _key: _newkey[_k] = (float(_k.strip('<>'))*1e3).__round__(4) ## to nm
			_newkey['Total Conc.'] = 'total'

			_df = _df[_key+['Total Conc.']].rename(_newkey,axis=1)
			_df['total'] = (_df.total.copy().map(lambda _: _.strip('(#/cm3)'))).astype(float)

		return _df
		 

	def __raw_process(self,_df,_freq):
		## customize each instrument
		out = _df.resample(_freq).mean().reindex(self.index(_freq))
		return out

	## read raw data
	def __reader(self,_reader,_nam):
		
		_path_raw = self.path/_nam
		_pkl_nam  = f'{_nam}.pkl'

		## read pickle if pickle file exisits and 'reset=False' or process raw data
		if (self.out_nam in listdir(self.path_out))&(~self.reset):
			print(f"\n\t{dtm.now().strftime('%m/%d %X')} : \033[96m{self.out_nam}\033[0m already exisits")
			return None
		else: 
			print(f"\n\t{dtm.now().strftime('%m/%d %X')} : Reading \033[96mRAW DATA\033[0m of {_nam} and process it")

		##=================================================================================================================
		## metadata parameter
		ext_nam, dt_freq = self.meta_read.values()

		## read raw data
		_df_con = None
		
		for file in listdir(_path_raw):
			if ext_nam not in file.lower(): continue
			print(f"\r\t\treading {file}",end='')

			_df = _reader(_path_raw/file)

			if _df is not None:
				_df_con = concat([_df_con,_df]) if _df_con is not None else _df

		_df_prcs = self.__raw_process(_df_con,dt_freq)
		print()

		return {'data':_df_prcs[_df_prcs.keys()[:-1]],'total':_df_prcs['total'].to_frame()}

	## get process data
	def save_data(self):

		_smps = self.__reader(self.__smps_reader,'smps')

		if (_smps is not None):
			_aps  = self.__reader(self.__aps_reader,'aps')

			with open(self.path_out/'smps_aps_raw.pkl','wb') as f:
				pkl.dump({'smps':_smps,'aps':_aps},f,protocol=pkl.HIGHEST_PROTOCOL)