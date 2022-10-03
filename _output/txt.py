# IMPORT PACKAGES
import numpy as np
from os import mkdir,listdir,remove


# IMPORT MODULES

# ------------------------------------------------------------------ txt_save()
def txt_save(opti,options,run):
	path = 'Output/txt/'+str(run)+'/'
	try:
		for f in listdir(path):
			remove(path+f)
	except:
		pass
	
	try:
		mkdir(path)
	except:
		pass

	nevals = len(opti.variables.all)
	evals = np.arange(1,nevals+1)
	best = opti.best.evolution
	nvars = len(opti.ref['initial'])
	ndig = len(str(nevals))
	ndig2 = len(str(opti.best.evolution[-1]))
	deli = 5
	spc = ' '.ljust(deli)
	evalt = 'EVAL'
	bestt = 'BEST'
	
	# EVOLUTION OF VARIABLES BY ALL EVALUATIONS
	lines = '{:s}'.format('-'*((18+deli)*nvars+ndig))
	header = [lines+'\n']
	title = []
	fmt = ['%-{:d}d'.format(ndig)] 
	for i in range(0,nvars):
		title.append('{:>18}'.format('VAR {:d}'.format(i+1)))
		fmt.append('%.12e')
	add = deli+ndig-len(str(evalt))
	title = spc.join(title)
	title = ' '.ljust(add)+ title
	header.append(evalt)
	header.append(title)
	header.append('\n'+lines)
	header = ''.join(header)

	fdata = np.column_stack((evals,opti.variables.all))
	with open(path+'Var_Evol_All.txt','w') as f:
		np.savetxt(f,fdata,fmt=fmt,header=header,delimiter=spc,comments='')

	# EVOLUTION OF VARIABLES BY BEST EVALUATIONS
	lines = '{:s}'.format('-'*((18+deli)*nvars+ndig+ndig2+deli))
	header = [lines+'\n']
	title = []
	fmt = ['%-{:d}d'.format(ndig2),'%-{:d}d'.format(ndig)] 
	for i in range(0,nvars):
		title.append('{:>18}'.format('VAR {:d}'.format(i+1)))
		fmt.append('%.12e')
	add = deli+ndig-len(str(evalt))
	title = spc.join(title)
	title = ' '.ljust(add) + title
	add2 = deli+ndig2-len(str(bestt))
	header.append(bestt)
	header.append(' '.ljust(add2) + evalt)
	header.append(title)
	header.append('\n'+lines)
	header = ''.join(header)

	fdata = np.column_stack((best,evals,opti.variables.best))
	with open(path+'Var_Evol_Best.txt','w') as f:
		np.savetxt(f,fdata,fmt=fmt,header=header,delimiter=spc,comments='')

	# EVOLUTION OF OBJECTIVE FUNCTION BY ALL EVALUATIONS OF EACH TEST
	ntests = len(options.info.name)

	total = opti.objfunc.all.total
	femu = opti.objfunc.all.femu
	vfm = opti.objfunc.all.vfm

	for n in range(0,ntests):
		fdata = np.column_stack((evals,total[:,n],femu[:,n],vfm[:,n]))
		fname = options.info.name[n]

		lines = '{:s}'.format('-'*((18+deli)*3+ndig))
		header = [lines+'\n']
		title = []
		titlenames = ['TOTAL','FEMU','VFM']
		fmt = ['%-{:d}d'.format(ndig)] 
		for i in range(0,3):
			title.append('{:>18}'.format(titlenames[i]))
			if ntests == 1:
				tmp1 = femu[0]
				tmp2 = vfm[0]
			else:
				tmp1 = femu[0,n]
				tmp2 = vfm[0,n]
			if i == 1:
				if np.isnan(tmp1):
					fmt.append('%18f')
				else:
					fmt.append('%.12e')
			elif i == 2:
				if np.isnan(tmp2):
					fmt.append('%18f')
				else:
					fmt.append('%.12e')
			else:
				fmt.append('%.12e')

		add = deli+ndig-len(str(evalt))
		title = spc.join(title)
		title = ' '.ljust(add)+ title
		header.append(evalt)
		header.append(title)
		header.append('\n'+lines)
		header = ''.join(header)
		
		with open(path+'Obj_Evol_All_{:s}.txt'.format(fname),'w') as f:
			np.savetxt(f,fdata,fmt=fmt,header=header,delimiter=' '.ljust(deli),
					   comments='')

	# EVOLUTION OF OBJECTIVE FUNCTION BY BEST EVALUATIONS OF EACH TEST
	ntests = len(options.info.name)

	total = opti.objfunc.best.total
	femu = opti.objfunc.best.femu
	vfm = opti.objfunc.best.vfm

	for n in range(0,ntests):
		fdata = np.column_stack((best,evals,total[:,n],femu[:,n],vfm[:,n]))
		
		fname = options.info.name[n]
		
		lines = '{:s}'.format('-'*((18+deli)*nvars+ndig+ndig2+deli))
		header = [lines+'\n']
		title = []
		titlenames = ['TOTAL','FEMU','VFM']
		fmt = ['%-{:d}d'.format(ndig2),'%-{:d}d'.format(ndig)]
		for i in range(0,3):
			title.append('{:>18}'.format(titlenames[i]))
			if ntests == 1:
				tmp1 = femu[0]
				tmp2 = vfm[0]
			else:
				tmp1 = femu[0,n]
				tmp2 = vfm[0,n]
			if i == 1:
				if np.isnan(tmp1):
					fmt.append('%16f')
				else:
					fmt.append('%.12e')
			elif i == 2:
				if np.isnan(tmp2):
					fmt.append('%18f')
				else:
					fmt.append('%.12e')
			else:
				fmt.append('%.12e')

		add = deli+ndig-len(str(evalt))
		title = spc.join(title)
		title = ' '.ljust(add) + title
		add2 = deli+ndig2-len(str(bestt))
		header.append(bestt)
		header.append(' '.ljust(add2) + evalt)
		header.append(title)
		header.append('\n'+lines)
		header = ''.join(header)
		
		with open(path+'Obj_Evol_Best_{:s}.txt'.format(fname),'w') as f:
			np.savetxt(f,fdata,fmt=fmt,header=header,delimiter=' '.ljust(deli),
					   comments='')

	# EVOLUTION OF OBJECTIVE FUNCTION BY ALL EVALUATIONS OF TOTAL
	total = opti.objfunc.all.total
	fdata = np.column_stack((evals,np.sum(total,axis=1)))

	lines = '{:s}'.format('-'*((18+deli)+ndig))
	header = [lines+'\n']
	title = []
	titlenames = ['TOTAL']
	fmt = ['%-{:d}d'.format(ndig)] 
	title.append('{:>18}'.format('TOTAL'))
	tmp1 = fdata[-1,0]
	fmt.append('%.12e')

	add = deli+ndig-len(str(evalt))
	title = spc.join(title)
	title = ' '.ljust(add)+ title
	header.append(evalt)
	header.append(title)
	header.append('\n'+lines)
	header = ''.join(header)
		
	with open(path+'Obj_Evol_All.txt','w') as f:
		np.savetxt(f,fdata,fmt=fmt,header=header,delimiter=' '.ljust(deli),
				   comments='')

	# EVOLUTION OF OBJECTIVE FUNCTION BY BEST EVALUATIONS OF TOTAL
	total = opti.objfunc.best.total
	fdata = np.column_stack((best,evals,np.sum(total,axis=1)))

	lines = '{:s}'.format('-'*((18+deli)+ndig+ndig2+deli))
	header = [lines+'\n']
	title = []
	titlenames = ['TOTAL']
	fmt = ['%-{:d}d'.format(ndig2),'%-{:d}d'.format(ndig)]
	title.append('{:>18}'.format('TOTAL'))
	fmt.append('%.12e')

	add = deli+ndig-len(str(evalt))
	title = spc.join(title)
	title = ' '.ljust(add) + title
	add2 = deli+ndig2-len(str(bestt))
	header.append(bestt)
	header.append(' '.ljust(add2) + evalt)
	header.append(title)
	header.append('\n'+lines)
	header = ''.join(header)
		
	with open(path+'Obj_Evol_Best.txt','w') as f:
		np.savetxt(f,fdata,fmt=fmt,header=header,delimiter=' '.ljust(deli),
				   comments='')