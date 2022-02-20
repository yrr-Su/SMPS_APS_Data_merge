from rawdata_process import reader
from data_merge.smps2aps import merger
from datetime import datetime as dtm
from datetime import timedelta as dtmdt
from pathlib import PurePath as Path




## time decorater
def __timer(func):
	def __wrap(*arg,**kwarg):
		print(f'\nPROGRAM : {__file__}\n')
		__st = dtm.now()

		## main function
		__out = func(*arg,**kwarg)

		__fn = dtm.now()
		__run = (__fn-__st).seconds

		print(f'\nProgram done\nrunning time = {__run//60:3d} min {__run%60:6.3f} s')
		return __out

	return __wrap


@__timer
def run(**kwarg):

	raw_process = reader(**kwarg)
	raw_process.save_data()

	merge_process = merger(**kwarg)
	merge_process.merge_data()
	merge_process.save_data()








#=============================================================================
if __name__=='__main__':

	# function 
	## run(path_raw, path_save, stara_time, end_time, reset=False, process_data_input=False)
	## path_raw   		  : str or path object
	##						Any valid string path is acceptable, input the path of raw data
	## path_data  		  : str or path object
	##						Any valid string path is acceptable, input the path of output data
	## start_time 		  : str or datetime-like
	##						Left bound for data index time
	## end_time  		  : str or datetime-like
	##						Right bound for data index time
	## reset   			  : bool, default is False
	##						Reload the raw data and re-produce the merge data
	## input_process_data : bool, default is False
	##						Use the processed data(made by YE, JUN-FA) which has been QC
	## 						rather than raw data

	# test - raw data
	# env paramerter
	'''
	ST_TIME = dtm(2021,11,8,18)
	ED_TIME = dtm(2021,11,9,16)

	PATH_SAVE = Path('..')/'data'/'test'
	PATH_RAW  = Path('..')/'data'/'test'/'raw_data'

	run(path_raw   = PATH_RAW,
		path_data  = PATH_SAVE,
		start_time = ST_TIME,
		end_time   = ED_TIME,
	)
	# '''

	# test2 - processed data
	# env paramerter
	# '''
	ST_TIME = '20200904 23:00:00'
	ED_TIME = '20210831 23:00:00'

	PATH_SAVE = Path('..')/'data'/'test2'
	PATH_RAW  = Path('..')/'data'/'test2'/'raw_data'

	run(path_raw   = PATH_RAW,
		path_data  = PATH_SAVE,
		start_time = ST_TIME,
		end_time   = ED_TIME,
		input_process_data = True
	)
	# '''