#!/usr/bin/env/python

from sys import argv


def create_defines_file():
	with open('inputs/params.opt','r') as f:
		temp = [x.split('+') for x in f.readlines()]
		defs=[]
		for x in temp:
			if x[0] == '' and '#' not in x:
				defs.append(x[-1].split('\n')[0])
	with open('src/defines.h','w') as g:
		if defs != []:
			for x in defs:
				g.write('#define ' + x + '\n')
		else:
			g.write('\n')
		if 'SELFGRAV2D' in defs and 'INFINITE' in defs:
			g.write('#define NEQNS 4\n')
		else:
			g.write('#define NEQNS 3\n')
	
				
	return defs			

def generate_extra_files(defs):

		
	print 'EDISK code set up with the flags'
	print defs
	if 'IMPLICIT' in defs and 'SPLIT' in defs:
		print 'Error cannot have both IMPLICIT and SPLIT flags'
		return ''
	
	if 'INFINITE' in defs:
		if 'SPLIT' in defs:
			print 'Error cannot have both INFINITE and SPLIT flags'
		if 'WAVEKILL' in defs or 'KILLIN' in defs or 'KILLOUT' in defs:
			print 'Error cannot have both INFINITE and WAVEKILL flags'
		if 'IMPLICIT' not in defs:
			print 'Error must have IMPLICIT enabled with INFINITE'
	
	
	
	if 'IMPLICIT' in defs or 'SPLIT' in defs:
		algfile = 'algo_driver.c implicit.c'
	else:
		algfile = 'algo.c rkf.c rk45.c'

	if 'SPLIT' in defs:
		algfile += ' rktvd.c'
	
	
	if 'INDIRECT' in defs:
		algfile += ' star.c'
	
	if 'COMPANION' in defs:
		algfile += ' companion.c'	
		
	if 'SELFGRAV' in defs:
		algfile += ' poisson.c'
	
	if 'SELFGRAV2D' in defs:
		algfile += ' poisson2d.c'
	
	print 'Adding ', algfile, 'to the Makefile'
	return algfile


def create_makefile(algfile):	
	with open("Makefile.in","r") as f:
		with open("Makefile","w") as g:
			for line in f.readlines():
				if line.split('=')[0] == 'SOURCES':
					files = line.split('=')[-1].split()
					files.append(algfile)
					files = ' '.join(files)
					files += '\n'
					line = '='.join([line.split('=')[0],files])
				g.write(line)
	return 0
	
	
defs = create_defines_file()	
algfile = generate_extra_files(defs)
create_makefile(algfile)


		