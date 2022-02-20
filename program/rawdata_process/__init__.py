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

	def __init__(self,path_raw,path_data,start_time,end_time,reset=False,input_process_data=False,**kwarg):
		print(f'\nSMPS and APS')
		print('='*65)
		print(f"Reading file and process data")

		## class parameter

		self.path  	   = Path(path_raw)
		self.path_out  = Path(path_data)

		self.meta_read = meta_dt['read']

		self.index = date_range(start_time,end_time,freq=self.meta_read['dt_freq'])
		start_time, end_time = self.index[[0,-1]]

		self.out_nam = 'smps_aps_raw.pkl'
		self.__time  = (start_time,end_time)

		self.reset = reset
		self.input_process_data = input_process_data
		
		print(f" from {start_time.strftime('%Y-%m-%d %X')} to {end_time.strftime('%Y-%m-%d %X')}")
		print('='*65)
		print(f"{dtm.now().strftime('%m/%d %X')}")

	def __smps_reader(self,_file):
		## customize each instrument
		## read one filess
		with open(_file,'r',encoding='utf-8',errors='ignore') as f:
			_df  = read_table(f,skiprows=18,parse_dates={'Time':['Date','Start Time']}).set_index('Time')
			_key = list(_df.keys()[7:-26])

			_newkey = {}
			for _k in _key: _newkey[_k] = float(_k).__round__(4)
			_newkey['Total Conc.(#/cm)'] = 'total'

		return _df[_key+['Total Conc.(#/cm)']].rename(_newkey,axis=1)

	def __aps_reader(self,_file):
		## customize each instrument
		## read one filess
		with open(_file,'r',encoding='utf-8',errors='ignore') as f:

			_df  = read_table(f,skiprows=6,parse_dates={'Time':['Date','Start Time']}).set_index('Time')
			_key = list(_df.keys()[3:54]) ## 542 ~ 1981

			_newkey = {}
			for _k in _key: _newkey[_k] = (float(_k.strip('<>'))*1e3).__round__(4) ## to nm
			_newkey['Total Conc.'] = 'total'

			_df = _df[_key+['Total Conc.']].rename(_newkey,axis=1)
			_df['total'] = (_df.total.copy().map(lambda _: _.strip('(#/cm3)'))).astype(float)

		return _df
		 

	def __raw_process(self,_df,_freq):
		## customize each instrument
		out = _df.resample(_freq).mean().reindex(self.index)
		return out

	## read raw data
	def __reader(self,_nam):
		print(f"\n\t{dtm.now().strftime('%m/%d %X')} : Reading \033[96mRAW DATA\033[0m of {_nam} and process it")
		
		_path_raw = self.path/_nam

		_reader = {'smps' : self.__smps_reader,
				   'aps'  : self.__aps_reader
				}
		_reader = _reader[_nam]

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

		return _df_prcs


	
	def __prcs_dt_reader(self,_file):
		## same templet for smps and aps
		
		
		
		return _df

	## read data
	## processed data, has been QC, should has same templet
	def __prcs_reader(self,):
		from pandas import read_excel

		print(f"\n\t{dtm.now().strftime('%m/%d %X')} : Reading \033[96mPROCESSED DATA\033[0m of {_nam}")
		
		## parameter
		_file = lambda _: listdir(self.path/_)[0]

		## read smps
		with open(_file('smps'),'rb') as f:
			_smps = read_excel(_file,parse_dates=['Time']).set_index('Time')

		## read aps
		## remove first key(<0.523)
		with open(_file('smps'),'rb') as f:
			_aps = read_excel(_file,parse_dates=['Time']).set_index('Time')
			_aps = _aps[_aps.keys()[1::]].copy()

		return _smps, _aps
		

	## get process data
	def save_data(self):


		## read pickle if pickle file exisits and 'reset=False' or process raw data
		if (self.out_nam in listdir(self.path_out))&(~self.reset):
			print(f"\n\t{dtm.now().strftime('%m/%d %X')} : \033[96m{self.out_nam}\033[0m already exisits")
			return None

		## read data
		if self.input_process_data:
			_smps, _aps = self.__prcs_reader()
		else:
			_smps, _aps = self.__reader('smps'), self.__reader('aps')


		## output data
		with open(self.path_out/'smps_aps_raw.pkl','wb') as f:
			pkl.dump({'smps':_smps,'aps':_aps},f,protocol=pkl.HIGHEST_PROTOCOL)
