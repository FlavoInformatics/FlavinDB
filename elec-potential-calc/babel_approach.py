import os

def command(file_name, start, end):

	# converts to an sdf format
	os.system("babel -ipdb {file_in} -osdf first_step.sdf".format(file_in = file_name))

	for idx in range(start-end):
		
		s_idx = str(start + idx)

		# makes seperate sdf files for relevant molecules ie flavins
		os.system("babel -isdf first_step.sdf -f {f_idx} -l {f_idx} -osdf step2{s}.sdf".format(f_idx = s_idx, s = str(idx))

		# makes seperate mol2 files for each flavin
		os.system("babel -isdf step2{}.sdf -omol2 final{}.mol2".format(str(idx)))
