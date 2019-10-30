# make batch / submission scripts for HPC jobs


def get_chunks(files, end=-1, start=-1, max=50):
	if start < 0 and end < 0:
		start = 1
		if files > max:
			chunks = files / max
			if files % max > 0:
				chunks = int(chunks) + 1
			if chunks == 0:
				chunks = 1
		else:
			chunks = 1
			end = files
	else:
		if end - start >= max:
			chunks = files / max
			if files % max > 0:
				chunks = int(chunks) + 1
		else:
			chunks = 1
	chunks = int(chunks)

	return chunks

def make_HPC_header(name='auto-ENRICH', system='BC3', nodes=1, ppn=1, walltime="100:00:00", mem='3'):

	strings = []

	if system == 'BC3':
		strings.append('# submission script for BC3')
		strings.append("#PBS -l nodes={0:<1d}:ppn={1:<1d}".format(nodes, ppn))
		strings.append("#PBS -l walltime={0:<9s}".format(walltime))
		strings.append("#PBS -l mem={0:<1d}GB".format(mem))
		strings.append("#PBD -N {0:<10s}".format(name))
		strings.append("cd $PBS_O_WORKDIR")
	elif system == 'BC4':
		# sbatch version
		strings.append('# submission script for BC4')
	elif system == 'localbox':
		strings.append('# submission script for local linux box')

	return strings

def make_HPC_orca_batch_submission(prefs, in_array, start, end, ck=0):
	strings = []
	strings.append('#!/bin/bash')
	if prefs['comp']['system'] == 'BC3':
		strings.append("#PBS -l nodes={0:<1d}:ppn={1:<1d}".format(nodes, ppn))
		strings.append("#PBS -l walltime={0:<9s}".format(walltime))
		strings.append("#PBS -l mem={0:<1d}GB".format(mem))
		if parallel:
			strings.append("#PBS -N {0:>1s}_{1:<1d}".format(molname, ck))
			strings.append("#PBS -t {0:>1d}-{1:<1d}".format(start, end))
			strings.append("cd $PBS_O_WORKDIR")
			strings.append("NMRNAME=$(gawk -v y=${{PBS_ARRAYID}} 'NR == y' {0:<5s})".format(in_array))
			strings.append("OUTNAME=$( echo $NMRNAME | sed 's/.in/.log/')")
			strings.append("orca ${NMRNAME} > ${OUTNAME}")
		else:
			strings.append("cd $PBS_O_WORKDIR")
			strings.append("for i in $(seq {0:>1d} 1 {1:<1d});".format(start, end))
			strings.append("do")
			strings.append("  NMRNAME=$(gawk -v y=${i} 'NR == y' {0:<5s})".format(in_array))
			strings.append("  base=${NMRNAME::${#NMRNAME}-4}")
			strings.append("  FILE=${base}.log")
			strings.append("OUTNAME=$( echo $NMRNAME | sed 's/.in/.log/')")
			strings.append("  if test -f '$FILE'; then")
			strings.append("    continue")
			strings.append("  else")
			strings.append("    orca ${NMRNAME} > ${OUTNAME}")
			strings.append("  fi")
			strings.append("done")

	elif prefs['comp']['system'] =='BC4':
		print('not done yet . . .')
	elif prefs['comp']['system'] == 'localbox':
		strings.append("NUMBERS=$(seq {0:>1d} {1:<1d})".format(start, end))
		strings.append("for NUM in ${NUMBERS}; do")
		strings.append("  NMRNAME=$(head -n${{NUM}}" + " {0:<5s} | tail -1)".format(in_array))
		strings.append("  OUTNAME=$( echo $NMRNAME | sed 's/.in/.log/')")
		strings.append("  orca ${NMRNAME} > ${OUTNAME}")
		strings.append("done")

	return strings
