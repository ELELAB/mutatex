#!/usr/env/python

#    mutatex: automate FoldX in-silico mutagenesis experiments
#    Copyright (C) 2015, Matteo Tiberti <matteo.tiberti@gmail.com>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
import os
import argparse
import subprocess as sp
import multiprocessing as mp
import logging as log
import shutil
import re
import numpy as np
from Bio import PDB

class ResList:
	def __init__(self, reslist=None, fname=None):

		if reslist and fname:
			log.warning("Warning; mutation list AND file specified. Mutation list will be ignored.")
		if fname:
			self.reslist = self.parse_list_file(fname=fname)
			return
		if reslist:
			self.reslist = tuple(reslist)
			return
		self.reslist = tuple()

	def __repr__(self):
		return ", ".join(self.reslist)

	__str__ = __repr__


	def parse_list_file(self, fname):
		try: 
			fh = open(fname, 'r')
		except:
			log.error("Couldn't open mutation list file %s." % fname)
			raise

		restypes = []
		for line in fh:
			if line and not line.startswith("#"):
				restypes.append(line.strip()[0])

		fh.close()
		return tuple(restypes)


class FoldXVersion:
	def __init__(self, binary=None):
		if binary:
			self.binary = os.path.abspath(binary)
		else:
			self.binary = None

	version = ""
	runfile_string = ""

	out_ext = "fxout"

	repaired_pdb_prefix = "RepairPDB_"
	average_fxout_prefix = "Average_BuildModel_"
	dif_fxout_prefix = "Dif_BuildModel_"

	mut_list_file = "individual_list.txt"
	
	len_dif_file_header = 9

	def parse_fxout(self, fname):
		pattern = '(\s+([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)){12}'

		energies = []
		log.info("Parsing file %s ..." %fname)
		try: 
			fh = open(fname, 'r')
		except:
			log.warning("Couldn't open file %s." % fname)
			raise

		for line in fh:
			if re.search(pattern, line.strip()):
				energies.append(float(line.strip().split()[2]))

		if len(energies) < 1:
			fh.close()
			log.warning("No energy values found in file %s!" % fname)
			raise

		fh.close()
		return energies

class FoldXVersion3b6(FoldXVersion):
	version="3b6"
	runfile_string = "-runfile"

class FoldXRun:

	ready = False

	finished = False

	runfile_name = "runfile.txt"

	runfile_content = ""

	def __init__(self, 
				name, 
				foldx_binary, 
				foldx_version, 
				base_directory, 
				pdbs, 
				runfile_name="runfile.txt",
				runfile_content="",
				runfile_processing={},
				prepare_finalization={},
				output_processing={},
				rotabase=None, 
				link_files=False,
				write_log=False):

		self.name = name
		self.pdbs = pdbs
		self.base_directory = os.path.abspath(base_directory)
		self.working_directory = self.base_directory+"/"+self.name
		self.runfile_name = runfile_name
		self.runfile_content = runfile_content
		self.runfile_processing = runfile_processing
		self.prepare_finalization = prepare_finalization
		self.foldx_binary = foldx_binary
		self.foldx_version = foldx_version
		self.output_processing = output_processing
		self.rotabase = rotabase
		self.link_files = link_files
		self.write_log = write_log
		self.ready = False


	def prepare(self):

		# Check if base directory exists
		if not os.path.exists(self.base_directory):
			log.warning("base directory %s does not exist; run %s will be skipped." % (self.base_directory, self.name))
			self.ready = False
			return False

		# Check status. Possible statuses: "already_done", "conflicting", "not_done"
		status = self.check_status()

		if status == "already_done":
			self.finished = True
			return True

		elif status == "conflicting":
			self.finished = False
			return False

		# Create the working directory 
		try: 
			safe_makedirs(self.working_directory, doexit=False)
		except:
			log.warning("Could not create working directory %s; run %s will be skipped." % (	self.working_directory, self.name))
			self.ready = False
			return False

		# Copy the PDB file(s) in the working directory
		for pdb in self.pdbs:
			try:
				safe_cp(os.path.abspath(pdb), self.working_directory+"/"+os.path.basename(os.path.abspath(pdb)), dolink=self.link_files, doexit=False)
			except:
				log.warning("Couldn't copy essential files for run %s; it will be skipped." % self.name)
				self.ready = False
				return False

		if self.rotabase:
			try:
				safe_cp(os.path.abspath(self.rotabase), self.working_directory+"/"+os.path.basename(os.path.abspath(self.rotabase)), dolink=self.link_files, doexit=False)
			except:
				log.warning("Couldn't copy essential rotabase.txt for run %s; it will be skipped." % self.name)
				self.ready = False
				return False

		# process runfile (if required)
		try: 
			self.process_runfile(**self.runfile_processing)
		except:
			log.warning("Couldn't process the runfile as required for run %s; it will be skipped." % self.name)
			self.ready = False
			return False

		# Write runfile
		try:
			with open(self.working_directory+"/"+self.runfile_name, 'w') as fh:
				fh.write(self.runfile_content)
		except:
			log.warning("Couldn't write runfile for run %s; it will be skipped." % self.name)
			self.ready = False
			return False

		try:
			self.finalize_prepare(**self.prepare_finalization)
		except:
			log.warning("Couldn't finalize preparation step for run %s; it will be skipped." % self.name)
			self.ready = False
			return False

		# Finally ...
		self.ready = True
		return True

	def run(self):
		if self.finished:
			log.warning("Run %s has already been performed; will not be run" % self.name)
			return True

		if not self.ready:
			log.warning("Run %s is not ready; will not be run." % self.name)
			return False

		runline = [self.foldx_binary, self.foldx_version.runfile_string, self.runfile_name]
		log.info("now running: %s" % self.name)

		if self.write_log:
			this_stdout_str = self.working_directory+"/"+'FoldXrun.log'
		else:
			this_stdout_str = os.devnull

		with open(this_stdout_str, 'w') as this_stdout:
			returncode = sp.call(runline, cwd=self.working_directory, stderr=sp.STDOUT, stdout=this_stdout)

		if returncode == 0:
			log.info("run %s has completed successfully" % self.name)
			self.finished = True
			return True
		else: 
			log.warning("FoldX exited with error for run %s!" % self.name)
			return False

		self.process_output(**self.output_processing)

	def status_check(self, **kwargs):
		if os.path.exists(self.working_directory):
			log.warning("working directory %s already exists." % self.working_directory)
			return "not_done"
			for f in os.listdir(self.working_directory):
				if f.endswith(".fxout"):
					log.warning(".fxout files already exist in directory %s; run %s will be skipped" % (self.working_directory, self.name))
					self.finished = True
					return "conflicting"

		return "not_done"

	def process_runfile(self, **kwargs):
		pass

	def process_output(self, **kwargs):
		pass

	def finalize_prepare(self, **kwargs):
		pass

class FoldXRepairRun(FoldXRun):
	def check_status(self):
		if os.path.exists(self.working_directory):
			log.warning("working directory %s already exists." % self.working_directory)
			if self.foldx_version.repaired_pdb_prefix + self.pdbs[0] in os.listdir(self.working_directory):
				log.warning("PDB output file already present; run %s will be skipped" % self.name)
				return "already_done"
			else:
				return "not_done"
		else:
			return "not_done"
		
	def process_runfile(self, **kwargs):
		pdbs = [os.path.basename(i) for i in self.pdbs]
		pdbstring = ",".join(pdbs)
		self.runfile_content = self.runfile_content.replace('$PDBS$',pdbstring)

class FoldXMutateRun(FoldXRun):
	def check_status(self):
		# Check if working directory exists
		if os.path.exists(self.working_directory):
			log.warning("working directory %s already exists." % self.working_directory)
			
			mutlist = []
			if True: #try:
				with open(self.working_directory+"/"+self.foldx_version.mut_list_file, 'r') as fh:
					for line in fh:
						if line:
							mutlist.append(line.strip()[-2])
			else: #except:
				log.warning("Couldn't open mutation file, so it's not clear what's in here. Run will be repeated.")
				return "not_done"

			if tuple(mutlist) != self.prepare_finalization['mutlist']:
				log.warning("mutation list file does not agree with the current mutation list! Run will be skipped.")
				return "conflicting"

			dif_files = [f for f in os.listdir(self.working_directory) if f.startswith(self.foldx_version.dif_fxout_prefix) and f.endswith(self.foldx_version.out_ext)]

			if len(dif_files) > 1:
				log.warning("More than one Dif_*.fxout files present in dir! Something fishy is going on. Run will be skipped.")
				return "conflicting"
			elif len(dif_files) == 0:
				return "not_done"
			if True:
				with open(self.working_directory+"/"+dif_files[0],'r') as fh:
					dif_file_len = len(fh.readlines())
					if (dif_file_len - self.foldx_version.len_dif_file_header) != self.runfile_processing["nruns"]*len(mutlist):
						log.warning("The  Dif_*.fxout file wasn't large as expected. The run didn't complete; it will be rerun.")
						return "not_done"
					else:
						return "already_done"
			else:
				log.warning("Couldn't read Dif_*.fxout out file; run will be skipped")
				return "conflicting"

			mutlist = []
		else:
			return "not_done"

	
	def process_runfile(self, **kwargs):
		pdbs = [os.path.basename(i) for i in self.pdbs]
		pdbstring = ",".join(pdbs)		
		self.runfile_content = self.runfile_content.replace('$PDBS$',pdbstring)
		self.runfile_content = self.runfile_content.replace('$NRUNS$', str(kwargs["nruns"]))

	def finalize_prepare(self, **kwargs):
		mutations_filename = self.foldx_version.mut_list_file

		try:
			fh = open(self.working_directory+"/"+mutations_filename, 'w')
		except:
			log.error("couldn't write file %s for run %s" % (mutations_filename, self.name))
			raise

		for res in kwargs['mutlist']:
			fh.write("%s%s;\n" % (self.name, res))
		fh.close()


def safe_makedirs(dirname, doexit=True):
	if os.path.exists(dirname):
		if not os.path.isdir(dirname):
			if doexit:
				log.error("%s exists but is not a directory. Exiting..." % dirname) 
				exit(1)
			else:
				log.warning("%s exists but is not a directory" % dirname) 
				raise
		else:
			log.warning("directory %s already exists" % dirname)
			return
	else:
		try: 
			os.makedirs(dirname)
		except:
			if doexit:
				log.error("Could not create directory %s. Exiting..." % dirname )
				exit(1)
			else:
				log.warning("Could not create directory %s." % dirname)
				raise


def safe_cp(source, destination, dolink=True, doexit=True):

	if dolink:
		verb = "link"
	else:
		verb = "copy"

	if not os.path.exists(source):
		log.error("Couldn't %s file %s; no such file or directory" % (verb, source)) 
		if doexit:
			exit(1)
		else:
			raise

	if not os.path.isfile(source):
		log.error("Couldn't %s file %s; it is not a file" % (verb, source))
		if doexit:
			exit(1)
		else:
			raise

	if not dolink:
		if os.path.exists(destination):
			log.warning("Destination file %s already exist; it will be overwritten." % destination)
		try:
			shutil.copyfile(source, destination)
		except:
			log.error("Couldn't copy file %s to %s" % (source, destination))
			if doexit:
				exit(1)
			else:
				raise
	else:
		if os.path.exists(destination):
			log.warning("Destination file %s already exist; it will not be overwritten by a link." % destination)
		else:
			try:
				os.symlink(source, destination)
			except:
				log.error("Couldn't link file %s to %s" % (source, destination))
				if doexit:
					exit(1)
				else:
					raise


def load_models(pdb, check_models=False):
	pdb_parser = PDB.PDBParser()

	log.info("loading pdb file %s" % pdb)

	try:
		structure = pdb_parser.get_structure('structure',pdb)
	except:
		log.error("Couldn't open file %s." % pdb)
		raise 

	if len(structure.get_list()) == 0:
		log.error("File %s doesn't contain any useful model." % pdb)
		raise

	if check_models:
		log.info("checking models in pdb file")
		for model in structure:
			for chain in model:
				if chain.id == ' ':
					log.warning('at least one residue in model %d in pdb file has no chain identifier. Will be defaulted to A.' % model.id)
					chain.id = 'A'

	return structure

def load_runfile(runfile):
	try:
		with open(runfile, 'r') as fh:
			data = fh.read()
	except:
		log.error("Couldn't open file %s." % runfile)
		raise 

	return data 

def foldx_worker(work_queue, done_queue):
	for run in iter(work_queue.get, 'STOP'):
		log.info("starting FoldX run %s" % run.name)
		run_out = run.run()
		done_queue.put((run.name, run_out))
	return True

def parallel_foldx_run(foldx_runs, np):
	work_queue = mp.Queue()
	done_queue = mp.Queue()

	for run in foldx_runs:
		work_queue.put(run)

	processes = []
	for w in xrange(np):
		work_queue.put('STOP')
		p = mp.Process(target = foldx_worker, args=[work_queue, done_queue])
		p.start()
		processes.append(p)

	for p in processes:
		p.join()

	out = []
	done_queue.put('STOP')

	for status in iter(done_queue.get, 'STOP'):
		out.append(status)

	return out

def split_pdb(filename, structure, checked):
	pdb_list = []
	writer = PDB.PDBIO()
	for model in structure:
		tmpstruc = PDB.Structure.Structure('structure')
		tmpstruc.add(model)
		writer.set_structure(tmpstruc)
		if checked:
			checked_str = "_checked"
		else:
			checked_str = ""
		writer.save("%s_model%d%s.pdb" % (os.path.splitext(filename)[0], model.id, checked_str))
		pdb_list.append("%s_model%d%s.pdb" % (os.path.splitext(filename)[0], model.id, checked_str))

	return pdb_list

def get_foldx_sequence(pdb):
	pdb_parser = PDB.PDBParser()

	log.info("loading pdb file %s" % pdb)

	try:
		structure = pdb_parser.get_structure('structure',pdb)
	except:
		log.error("Couldn't open file %s." % pdb)
		raise 

	residue_list = []
	for model in structure:
		for chain in model:
			for residue in chain:
				try:
					res_code = PDB.Polypeptide.three_to_one(residue.get_resname())
				except: 
					log.warning("Residue %s in file %s couldn't be recognized; it will be skipped" %(residue, pdb))
					continue
				res_number = residue.get_id()[1]
				residue_list.append("%s%s%d" % (res_code, chain.get_id(), residue.get_id()[1]))

	return tuple(residue_list)

def save_energy_file(fname, data, fmt="%.5f"):
	try:
		np.savetxt(fname, data, fmt=fmt)
	except:
		log.warning("Couldn't write file %s" % fname)	

# Gather info about FoldX version


def main():

	foldx_versions = [FoldXVersion3b6]

	supported_foldx_versions = {v.version:v for v in foldx_versions}

	foldx_binary_var = os.getenv('FOLDX_BINARY')
	foldx_rotabase_var = os.getenv('FOLDX_ROTABASE')
	
	if foldx_binary_var is None:
		foldx_binary_var = ""
	if foldx_rotabase_var is None:
		foldx_rotabase_var = ""

	parser = argparse.ArgumentParser(description='Setup and run in silico saturation mutagenesis with FoldX.')

	parser.add_argument('pdb', metavar='PDBFILE', type=str, nargs='+',  )
	parser.add_argument('--skip-check', '-r', dest='skip_check', default=False, action='store_true', help="Skip PDB checking phase")
	parser.add_argument('--skip-repair', dest='skip_repair', default=False, action='store_true', help="Skip PDB repair phase")
	parser.add_argument('--skip-mutate', dest='skip_mutate', default=False, action='store_true', help="Skip PDB mutation phase")
	parser.add_argument('--skip-report', dest='skip_report', default=False, action='store_true', help="Skip reporting phase")
	parser.add_argument('--reslist', '-s', dest='reslist', default=None, action='store', type=str, help="Mutation list file")
	parser.add_argument('--multiple-models', '-u', dest='multiple_models', default=False, action='store_true', help="Indipendently use mulitple models found in the PDB")
	parser.add_argument('--nruns', '-n', dest='nruns', default=5, action='store', type=int, help="number of FoldX mutation runs per mutation (default: 5)")
	parser.add_argument('--np', dest='np', default=1, type=int, help="Number of FoldX processes to be run at the same time")
	parser.add_argument('--mutlist','-m', dest="mutlist", default=None, type=str, help="File containing the residue types that each ")
	parser.add_argument('-x', '--foldx-binary', dest="foldx_binary", action='store', type=str, help="Location of the FoldX binary (default: content of the FOLDX_BINARY system variable", default=foldx_binary_var)
	parser.add_argument('--rotabase', dest="rotabase", action="store", type=str, help="Location of the FoldX rotabase.txt file (default: content of the FOLDX_ROTABASE system variable", default=foldx_rotabase_var)
	parser.add_argument('--foldx-version', dest="foldx_version", action='store', choices=supported_foldx_versions.keys(), default=supported_foldx_versions.keys()[0], help="FoldX version to be used (possible options: %s" % ", ".join(supported_foldx_versions.keys()))
	parser.add_argument('-v', '--verbose', dest="verbose", default=False, action='store_true', help="Verbose mode")
	parser.add_argument('--foldx-log', dest="write_log", default=False, action='store_true', help="Write FoldX standard output on file for each run")
	parser.add_argument('-l','--use-links', dest="use_links", default=False, action='store_true', help="Use links instead of copying files as much as possibile")
	parser.add_argument('--repair-runfile-template','--repair', dest="repair_runfile_template", type=str, default="repair_runfile_template.txt", help="Template runfile for repair runs (default: ./repair_runfile_template.txt)")
	parser.add_argument('--mutate-runfile-template','--mutate', dest="mutate_runfile_template", type=str, default="mutate_runfile_template.txt", help="Template runfile for mutation runs (default: ./mutate_runfile_template.txt)")


# XXX: remove default verbose mode
	log.basicConfig(level=log.INFO)

	args = parser.parse_args()

# Test whether the FoldX binary is available and executable
	if not args.foldx_binary:
		log.error("The FoldX binary must be provided, either by setting the FOLDX_BINARY system variable to the path of the executable or by using the --foldx-binary option. Exiting...")
		exit(1)

	if not os.path.isfile(args.foldx_binary):
		log.error("The specified FoldX binary (%s) is not a file. Exiting..." % args.foldx_binary)
		exit(1)

	if not os.access(args.foldx_binary, os.X_OK):
		log.error("The specified FoldX binary (%s) is not executable. Exiting..." % args.foldx_binary)
		exit(1)

# Test whether the rotabase.txt file exists and is not a link
	if not os.path.isfile(args.rotabase):
		log.error("The rotabase.txt (%s) is not a file or could not be found. Exiting..." % args.rotabase)
		exit(1)

	if os.path.islink(args.rotabase):
		log.warning("The provided rotabase.txt file is a link. Be advised: in some occasions, FoldX has been reported to misbehave when rotabase.txt is provided as a link.")

# Set up FoldX version object
	if args.foldx_version not in supported_foldx_versions.keys():
		log.error("FoldX version %s not supported by this release. Exiting...")
		exit(1)

		current_version = supported_foldx_versions[args.foldx_version](binary=args.foldx_binary)

# defaults 
		repair_dirname = "repair"
		mutations_dirname = "mutations"
		results_dirname = "results"
		averages_dirname = "final_averages"
		DEVNULL = open(os.devnull, 'w')
		default_mutlist = ('G','A','V','L','I','M','F','W','P','S','T','C','Y','N','Q','D','E','K','R','H')

		pdb_models = []
		try:
			for pdb in args.pdb:
				pdb_models.append(load_models(pdb, check_models = not args.skip_check))
		except:
			exit(1)

	for p in range(len(pdb_models)):
		if len(pdb_models[p].get_list()) > 1 and not args.multiple_models:
			log.error("input pdb %s has more than one model, and --multiple-models was not selected. Please use --multiple-models or input a single-model PDB file." % args.pdb[p])
			exit(1)

# PHASE ONE: prepare PDB file(s)
	log.info("starting phase 1: structure check and generation of required PDB files")

	main_dir = os.getcwd()

	log.info("working directory will be: %s" % main_dir)


	pdbs_list = []

	for p in range(len(pdb_models)):
		pdbs_list.extend( split_pdb(args.pdb[p], pdb_models[p], not args.skip_check) )

	log.info("starting phase 2: REPAIR")

	working_directory = main_dir+"/"+repair_dirname

	safe_makedirs(working_directory)

	log.info("working dir: %s" %working_directory)

	repair_runs = []

	for pdb in pdbs_list:
		repair_runs.append(FoldXRepairRun(name = "repair_"+os.path.splitext(pdb)[0],
										foldx_binary = current_version.binary,
										foldx_version = current_version,
										base_directory = working_directory,
										pdbs = [pdb],
										runfile_content=load_runfile(args.repair_runfile_template),
										rotabase = args.rotabase,
										link_files = args.use_links,
										write_log = args.write_log,
										))


# Prepare repair runs

	if not args.skip_repair:
		for r in repair_runs:
			log.info("Preparing for repair run %s" % (r.name))
			r.prepare()

	# Run repair runs
		log.info("Running repair runs")
		repair_outcome = parallel_foldx_run(repair_runs, np=args.np)

		if False in zip(*repair_outcome)[1]:
			log.error("The following repair runs failed to complete. Exiting...")
		for i in filter(lambda x: x[1] == False, repair_outcome):
			log.error("\t%s" % i[0])
	else:
		log.info("Repair phase skipped, as requested.")

# PHASE TWO: mutate + energy

# Prepare mutation list
	if args.mutlist:
		log.info("File %s will be used as mutation list" % args.mutlist)
		try:
			mutlist = ResList(fname=args.mutlist)
		except:
			exit(1)
	else:
		log.info("default mutation list will be used.")
		mutlist = ResList(reslist=default_mutlist)

	log.info("mutation list is: %s" % mutlist)

# Select PDB(s) to be used
	if args.skip_repair:
		repaired_pdbs_list = args.pdb
	else:
		repaired_pdbs_list = [repair_runs[i].working_directory+"/"+current_version.repaired_pdb_prefix+pdbs_list[i] for i in range(len(pdbs_list))]
	log.info("list of PDBs to be used: %s" % ", ".join(repaired_pdbs_list))

# create working directory
	working_directory = main_dir+"/"+mutations_dirname
	log.info("Working directory is: %s" % working_directory)

	safe_makedirs(working_directory)

	residues_list = []
	dir_list = []
	for pdb in repaired_pdbs_list:
		this_dirname = working_directory+"/"+os.path.splitext(os.path.basename(pdb))[0] 
		try:
			safe_makedirs(this_dirname)
		except:
			logging.error("PDB %s will be ignored." % pdb)
			repaired_pdbs_list.remove(pdb)
			continue
		residues_list.append(get_foldx_sequence(pdb))
		log.info("Sequence for pdb %s is: %s)" % (pdb, ",".join(residues_list[-1])))
		    

	residue_sets = [set(l) for l in residues_list]
	for i in range(len(residue_sets)):
		for j in range(len(residue_sets)):
			if residue_sets[i] != residue_sets[j]:
				log.warning("PDB %s has a different sequence respect to %s!" % (repaired_pdbs_list[i], repaired_pdbs_list[j]))


	unique_residues = sorted(list(set.intersection(*residue_sets)), key=lambda x: ord(x[1])*10**12+int(x[2:]))

	if len(unique_residues) == 0:
		logging.error("ERROR: The models provided don't have any residue in common! Exiting...")
		exit(1)

	mutation_runs = []

	for r in unique_residues:
		this_pdbs = []
		for i in range(len(residue_sets)):
			if r in residue_sets[i]:
				this_pdbs.append(repaired_pdbs_list[i])
	#log.info("these PDBs will be used for run %s:\n\t%s" % (r, "\n\t".join(this_pdbs)))
		for pdb in this_pdbs:
			this_workdir = working_directory+"/"+os.path.splitext(os.path.basename(pdb))[0]
			mutation_runs.append(FoldXMutateRun(name = r,
										foldx_binary = current_version.binary,
										foldx_version = current_version,
										base_directory = this_workdir,
										pdbs = [pdb],
										runfile_content = load_runfile(args.mutate_runfile_template),
										rotabase = args.rotabase,
										link_files = args.use_links,
										write_log = args.write_log,
										runfile_processing = {"nruns":args.nruns},
										prepare_finalization = {"mutlist":mutlist.reslist}
										))

	if not args.skip_mutate:

		for r in mutation_runs:
			log.info("Preparing for mutation run %s" % r.name)
			r.prepare()

		log.info("Running mutate runs")

		mutate_outcome = parallel_foldx_run(mutation_runs, np=args.np)
		if False in zip(*mutate_outcome)[1]:
			log.warning("The following mutation runs failed to complete. The corresponding positions in sequence will be skipped.")
		for i in filter(lambda x: x[1] == False, mutate_outcome):
			log.warning("\t%s" % i[0])
			unique_residues.remove(i[0])
	else:
		log.info("Mutation phase skipped, as requested")

	working_directory = main_dir+"/"+results_dirname
	log.info("Working directory is: %s" % working_directory)

	if not args.skip_report:
		safe_makedirs(working_directory)

		for pdb in repaired_pdbs_list:
			this_pdb_dir = working_directory+"/"+os.path.splitext(os.path.basename(pdb))[0]
			safe_makedirs(this_pdb_dir)
		safe_makedirs(working_directory+"/"+averages_dirname)
	
		for res in unique_residues:
			energies = []
			this_runs = filter(lambda x: x.name == res, mutation_runs)
			for r in this_runs:
				for pdb in r.pdbs:
					results_basedir = os.path.splitext(os.path.basename(pdb))[0]
					try:
						energies.append( current_version.parse_fxout(r.working_directory+"/"
														+current_version.average_fxout_prefix
														+results_basedir
														+".fxout"))
					except:
						log.warning("Couldn't parse energy file for PDB %s; mutation site %s will be skipped." % (pdb, r.name))
						continue
					save_energy_file(working_directory+"/"+results_basedir+"/"+r.name, energies[-1])
			save_energy_file(working_directory+"/"+averages_dirname+"/"+r.name, np.average(energies, axis=0))
	else:
		log.info("Reporting phase was skipped, as requested.")

	log.info("All done!")

if __name__ == '__main__':
	main()
