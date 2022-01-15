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
def run(path_data,start,end,**kwarg):

	default_par = { 'path_rawdata' : path_data/'raw_data',
					}

	raw_process = reader(default_par['path_rawdata'],path_data,start,end,reset=0)
	raw_process.save_data()

	merge_process = merger(path_data,start,end)
	merge_process.merge_data()








#=============================================================================
if __name__=='__main__':

	# env paramerter
	ST_TIME = dtm(2021,11,8,18)
	ED_TIME = dtm(2021,11,9,16)

	PATH_DATA = Path('..')/'data'/'test'





	run(PATH_DATA,ST_TIME,ED_TIME)



