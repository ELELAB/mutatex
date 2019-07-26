#!/usr/bin/env python

# core.py: classes and functions for the main script
# Copyright (C) 2015, Matteo Tiberti <matteo.tiberti@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
import os
import argparse
import subprocess as sp
import logging as log
import re
import numpy as np
from Bio import PDB
from six import iteritems
from mutatex.utils import * 

class MutationList(object):
    def __init__(self, res_groups, mutations, name="", selfmutate=False):
        self.name = name
        if selfmutate:
            self.res_groups = tuple([res_groups])
            mutations = []
            for rg in self.res_groups:
                mutations.append(tuple([ rrg[0] for rrg in rg]))
            self.mutations = tuple(mutations)
            return

        if type(mutations) is ResList:
            self.res_groups = tuple([res_groups] * len(mutations.reslist))
            self.mutations = tuple([ tuple([mutations.reslist[i]] * len(self.res_groups[i])) for i in range(len(self.res_groups))])
        else:
            self.res_groups = tuple(res_groups)
            self.mutations = tuple(mutations)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented
        if set(self.res_groups) != set(other.res_groups):
            return False
        if set(self.mutations) != set(other.mutations):
            return False
        return True

    def __ne__(self, other):
        return not self == other


class ResList(object):
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

    def parse_list_file(self, fname):
        return parse_mutlist_file(fname)

class EnergyReport:
    def __init__(self, pdbs=None):
        self.residues = {}
        self.energies = {}
        if not pdbs is None:
            for pdb in pdbs:
                self.energies[pdb] = []
                self.residues[pdb] = []

    def add_residue(self, res, energy, pdb, do_avg=None, do_std=None, do_max=None, do_min=None):
        assert(energy.shape[0] == 1)
        energy = energy[0]

        if not pdb in list(self.energies):
            self.energies[pdb] = []
            self.residues[pdb] = []

        self.energies[pdb].append(energy)
        self.residues[pdb].append(res)

    def save(self, directory, fname="selfmutation_energies.dat", do_avg=True, do_std=True, do_min=True, do_max=True):

        header_cols = []

        if do_avg:
            header_cols.append("avg")
        if do_std:
            header_cols.append("std")
        if do_min:
            header_cols.append("min")
        if do_max:
            header_cols.append("max")

        dtype = [('res', 'U8')] + [(i, 'f') for i in header_cols]
        header = "\t".join(header_cols)
        fmt=["%8s"] + ["%10f" for f in header_cols]

        for pdb,energies in iteritems(self.energies):

            out = [self.residues[pdb]]
            print(out)

            if do_avg:
                out.append(np.average(energies, axis=1))
            if do_std:
                out.append(np.std(energies, axis=1))
            if do_min:
                out.append(np.min(energies, axis=1))
            if do_max:
                out.append(np.max(energies, axis=1))

            out = np.array(list(zip(*out)), dtype=dtype)

            try:
                np.savetxt(os.path.join(directory, pdb, fname), out, fmt=fmt, header=header)
            except:
                log.error("Couldn't write file %s" % fname)

class FoldXVersion:
    def __init__(self, binary=None, rotabase=None):
        if not binary is None:
            self.binary = os.path.abspath(binary)
        else:
            self.binary = None

        if not rotabase is None:
            self.rotabase = os.path.abspath(rotabase)
        else:
            self.rotabase = None


    version = None
    runfile_string = None

    out_ext = "fxout"
    mut_list_file = "individual_list.txt"

class FoldXSuiteVersion4(FoldXVersion):
    version="suite4"
    runfile_string = "-f"
    can_generate_rotabase = False

    len_dif_file_header = 9
    pdblist_fxout_prefix = "PdbList_"
    summary_fxout_suffix= "_AC"
    repaired_pdb_prefix = "RepairPDB_"
    average_fxout_prefix = "Average_BuildModel_"
    dif_fxout_prefix = "Dif_"
    summary_fxout_prefix = "Summary_"
    mutation_output_pdb_WT_prefix = "WT_"
    mutation_output_pdb_prefix = ""
    mutlist_eol = ";\n"

    def save_mutlist(self, fname, mutlist):

        try:
            fh = open(fname, 'w')
        except:
            log.error("Couldn't open file %s for writing!" % fname)
            return None

        for i, rg in enumerate(mutlist.res_groups):
            fh.write(",".join(["%s%s" % (rg[j], mutlist.mutations[i][j]) for j,r in enumerate(rg) ]))
            fh.write(self.mutlist_eol)

        fh.close()

    def parse_mutlist(self, fname):

        try:
            fh = open(fname, 'r')
        except:
            log.error("Couldn't open file %s for reading!" % fname)
            raise IOError

        res_groups = []
        mutations = []

        for line in fh:
            if line:
                tmp = line.strip(self.mutlist_eol).split(",")
                for t in tmp:
                    res_groups.append(tuple([t[:-1] for t in tmp]))
                    mutations.append(tuple([t[-1] for t in tmp]))

        mutations = tuple(mutations)
        res_groups = tuple(res_groups)
        ml = MutationList(res_groups=res_groups, mutations=mutations)
        return ml

    def repair_pdb_output_fname(self, basename):
        return "%s_Repair.pdb" % "".join(os.path.splitext(basename)[:-1])

    def mutate_average_fxout_output_fname(self, basename, *args, **kwargs):
        return "Average_%s.fxout" % basename

    def mutate_dif_fxout_output_fname(self, basename, *args, **kwargs):
        return "Dif_%s.fxout" % basename

    def mutate_pdblist_fxout_output_fname(self, basename, *args, **kwargs):
        return "%s%s.fxout" % (self.pdblist_fxout_prefix, basename)

    def ac_summary_fxout_output_fname(self, basename, *args, **kwargs):
        return "%s%s%s.fxout" % (self.summary_fxout_prefix, basename, self.summary_fxout_suffix)

    def get_mutation_fxout_fnames(self, directory, pdbs):
        basenames = ["".join(os.path.splitext(os.path.basename(pdb))[:-1]) for pdb in pdbs]
        return [os.path.join(directory,self.mutate_dif_fxout_output_fname(basename)) for basename in basenames]

    def parse_mutations_fxout(self, this_run):

        pattern = '(\s+([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)){22}'

        directory = this_run.working_directory

        pdb = this_run.pdbs[0]

        fname = self.mutate_dif_fxout_output_fname("".join(os.path.splitext(os.path.basename(pdb))[:-1]))

        energies = []

        log.info("Parsing file %s ..." %fname)

        try:
            with open(os.path.join(directory,fname), 'r') as fh:
                for line in fh:
                    if re.search(pattern, line.strip()):
                        energies.append(float(line.strip().split()[1]))
        except:
            log.error("Couldn't open or parse file %s" % fname)

        if len(energies) < 1:
            log.error("No energy values found in energy files!")

        energies = np.array(energies)
        energies = energies.reshape(len(this_run.mutlist.mutations),
                                        energies.shape[0]//len(this_run.mutlist.mutations))

        return energies

    def get_mutation_pdb_fnames(self, directory, pdbs, run, WT=True, include_original=False):
        fnames = []
        files = os.listdir(directory)
        if WT:
            prefixes = (self.mutation_output_pdb_WT_prefix, self.mutation_output_pdb_prefix)
        else:
            prefixes = (self.mutation_output_pdb_prefix)

        for pdb in pdbs:
            results_basedir = os.path.splitext(os.path.basename(pdb))[0]
            fnames.append([[] for i in range(len(prefixes)+int(include_original))])

            for p,prefix in enumerate(prefixes):
                for i in range(1,len(run.mutlist.mutations)+1):
                    for j in range(run.runfile_processing['nruns']):
                        fnames[-1][p].append("%s/%s%s_%d_%d.pdb" % (directory, prefix, results_basedir, i, j))
                        fnames[-1][-1].append(pdb)

        if len(fnames) == 1:
            if len(fnames[0]) == 1:
                return fnames[0][0]
            return fnames[0]
        return fnames

    def get_interaction_fxout_fnames(self, directory, pdbs, run, original_pdb=False):

        fnames = [[],[]]

        for pdb in pdbs:
            pdb_basename = "".join(os.path.splitext(os.path.basename(pdb))[:-1])
            if self.mutation_output_pdb_WT_prefix in pdb:
                fnames[0].append(os.path.join(directory, self.ac_summary_fxout_output_fname(pdb_basename)))
            else:
                fnames[1].append(os.path.join(directory, self.ac_summary_fxout_output_fname(pdb_basename)))

        return fnames

    def parse_interaction_energy_summary_fxout(self, directory, pdbs, run):

        fnames = self.get_interaction_fxout_fnames(directory, pdbs, run, original_pdb=True)

        energies = {}
        delta_energies = {}
        groups = []
        types = ['wt','mutated']

        for i,prefix in enumerate(types):
            energies[prefix] = {}
            for fname in fnames[i]:
                log.info("Parsing file %s ..." %fname)
                try:
                    fh = open(fname, 'r')
                except:
                    log.warning("Couldn't open file %s." % fname)
                    raise

                fh_lines = fh.readlines()
                for line in fh_lines[self.len_dif_file_header:]:
                    tmp = line.strip().split()
                    idx = frozenset((tmp[1], tmp[2]))
                    if idx in list(energies[prefix]):
                        energies[prefix][idx].append(float(tmp[5]))
                    else:
                        energies[prefix][idx] = [float(tmp[5])]
                fh.close()

            for k,v in iteritems(energies[prefix]):
                v = np.array(v)
                energies[prefix][k] = v.reshape((len(run.mutlist.mutations),
                                      v.shape[0]//len(run.mutlist.mutations)))

        interaction_groups = set(list(energies[types[0]]) + list(energies[types[1]]))

        for ig in interaction_groups:
            delta_energies[ig] = energies[types[1]][ig] - energies[types[0]][ig]

        return delta_energies

    def check_dif_file_size(self, cwd, fname, nmuts, nruns):

        with open("%s/%s" % (cwd,fname)) as fh:
            fsize = len(fh.readlines())

        if fsize-self.len_dif_file_header == nruns*nmuts:
            return True
        return False

    def check_pdb_file_size(self, cwd, fname, nmuts, nruns):

        with open("%s/%s" % (cwd,fname)) as fh:
            fsize = len(fh.readlines())

        if fsize == 2*nruns*nmuts:
            return True
        return False

class FoldXSuiteVersion5(FoldXSuiteVersion4):
    version = "suite5"
    can_generate_rotabase = True

class FoldXRun(object):

    runfile_name = "runfile.txt"

    logfile_name = "FoldXrun.log"

    runfile_content = ""

    finished = False
    ready = False

    def __init__(self,
                name,
                foldx_version,
                base_directory,
                pdbs,
                runfile_name="runfile.txt",
                runfile_content="",
                runfile_processing={},
                prepare_finalization={},
                output_processing={},
                link_files=False,
                write_log=False,
                clean='partial',
                *args,
                **kwargs):

        self.name = name
        self.pdbs = pdbs
        self.base_directory = os.path.abspath(base_directory)
        self.working_directory = os.path.join(self.base_directory, self.name)
        self.runfile_name = runfile_name
        self.runfile_content = runfile_content
        self.runfile_processing = runfile_processing
        self.prepare_finalization = prepare_finalization
        self.foldx_version = foldx_version
        if output_processing is None:
            self.output_processing = {}
        else:
            self.output_processing = output_processing
        self.link_files = link_files
        self.write_log = write_log
        self.do_clean = clean

    def prepare(self):

        # Check if base directory exists
        if not os.path.exists(self.base_directory):
            log.warning("base directory %s does not exist; run %s will be skipped." % (self.base_directory, self.name))
            self.ready = False
            return False

        # Check status. Possible statuses: "already_done", "conflicting", "not_done", "broken", "interface_missing_data"
        status = self.check_status()

        if status == "interface_missing_data":
            log.warning("The corresponding mutation run %s hasn't been completed successfully; skipping" % self.name)
            self.ready = False
            return False

        if status == "already_done":
            self.finished = True
            return True

        elif status == "conflicting" or status == "broken":
            log.warning("Something went wrong with run %s; it will be skipped." % self.name)
            self.finished = False
            return False

        # Create the working directory
        try:
            safe_makedirs(self.working_directory)
        except:
            log.warning("Could not create working directory %s; run %s will be skipped." % (    self.working_directory, self.name))
            self.ready = False
            return False

        # Copy the PDB file(s) in the working directory
        for pdb in self.pdbs:
            try:
                #print 'source', os.path.abspath(pdb)
                #print 'dest', self.working_directory+"/"+os.path.basename(os.path.abspath(pdb))
                safe_cp(os.path.abspath(pdb), os.path.join(self.working_directory, os.path.basename(os.path.abspath(pdb))), dolink=self.link_files)
            except:
                log.warning("Couldn't copy essential files for run %s; it will be skipped." % self.name)
                self.ready = False
                return False

        if self.foldx_version.rotabase:
            try:
                safe_cp(os.path.abspath(self.foldx_version.rotabase), os.path.join(self.working_directory, os.path.basename(os.path.abspath(self.foldx_version.rotabase))), dolink=self.link_files)
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

        runline = [self.foldx_version.binary, self.foldx_version.runfile_string, self.runfile_name]
        log.info("now running: %s" % self.name)
        log.info("command is %s" % " ".join(runline))

        if self.write_log:
            this_stdout_str = self.working_directory+"/"+self.logfile_name
        else:
            this_stdout_str = os.devnull

        with open(this_stdout_str, 'w') as this_stdout:
            returncode = sp.call(runline, cwd=self.working_directory, stderr=sp.STDOUT, stdout=this_stdout)

        if returncode == 0:
            log.info("run %s has completed successfully" % self.name)
            self.finished = True

            self.process_output(**self.output_processing)
            return True
        else:
            log.warning("FoldX exited with error for run %s!" % self.name)
            return False


    def check_status(self, **kwargs):
        if not os.path.exists(self.working_directory):
            return "not_done"
        else:
            log.warning("working directory %s already exists." % self.working_directory)
            for f in os.listdir(self.working_directory):
                if f.endswith(".fxout"):
                    log.warning(".fxout files already exist in directory %s; run %s will be skipped" % (self.working_directory, self.name))
                    self.finished = True
                    return "conflicting"

        return "not_done"

    def process_runfile(self, **kwargs):
        pass

    def process_output(self, **kwargs):
        log.info("Processing output ...")
        if self.do_clean == 'deep':
            self.clean()
        elif self.do_clean == 'partial':
            self.partial_clean()
        elif self.do_clean == 'none':
            log.info("No cleaning will be performed at this stage.")


    def finalize_prepare(self, **kwargs):
        pass

    def partial_clean(self):
        log.info("Doing partial cleaning of the working directory %s" % self.working_directory)
        for f in os.listdir(self.working_directory):
            fname = os.path.join(self.working_directory, f)
            if f.startswith("WT_") and f.endswith(".pdb"):
                log.info("removing %s" % fname)
                os.remove(fname)

    def clean(self):
        log.info("Cleaning working directory %s" % self.working_directory)
        for f in os.listdir(self.working_directory):
            fname = os.path.join(self.working_directory, f)
            if os.path.isfile(fname) and not fname.endswith(".fxout") and not fname.endswith(".log") and not fname == os.path.basename(fname):
                log.info("removing %s" % fname)
                os.remove(fname)

    def reset_working_directory(self):
        ftypes = [".pdb", ".fxout", ".log"]
        individual_files = ["individual_list.txt", "rotabase.txt", "runfile.txt"]
        for f in os.listdir(self.working_directory):
            if os.path.splitext(f)[-1] in ftypes:
                try:
                    os.remove(os.path.join(self.working_directory, f))
                except:
                    pass
        for f in individual_files:
            full_path = os.path.join(self.working_directory, f)
            if os.path.isfile(full_path):
                log.info("removing %s" % os.path.join(self.working_directory, f))
                try:
                    os.remove(full_path)
                except:
                    pass

class FoldXRepairRun(FoldXRun):
    def check_status(self):
        if os.path.exists(self.working_directory):
            log.warning("working directory %s already exists." % self.working_directory)
            if self.foldx_version.repair_pdb_output_fname(self.pdbs[0]) in os.listdir(self.working_directory):
                log.warning("PDB output file already present; run %s will be skipped" % self.name)
                return "already_done"
        return "not_done"

    def process_runfile(self, **kwargs):
        pdbs = [os.path.basename(i) for i in self.pdbs]
        pdbstring = ",".join(pdbs)
        self.runfile_content = self.runfile_content.replace('$PDBS$',pdbstring)

    def clean(self):
        pass

class FoldXMutateRun(FoldXRun):

    logfile_name = "MutateFoldXRun.log"

    def __init__(self, mutlist=None, *args, **kwargs):
        super(FoldXMutateRun, self).__init__(*args, **kwargs)
        self.mutlist = mutlist

    def prepare(self):

        dir_reset = False

        # Check if base directory exists
        if not os.path.exists(self.base_directory):
            log.warning("base directory %s does not exist; run %s will be skipped." % (self.base_directory, self.name))
            self.ready = False
            return False

        # Check status. Possible statuses: "already_done", "conflicting", "not_done", "broken"
        status = self.check_status()

        if status == "already_done":
            self.finished = True
            return True

        elif status == "broken" or status == "conflicting" or status == "partially_done":
            self.finished = False
            log.warning("Working directory of run %s was left in an undefined state; it will be reset" % self.name)
            self.reset_working_directory()
            dir_reset = True

        # Create the working directory
        if not dir_reset:
            try:
                safe_makedirs(self.working_directory)
            except:
                log.warning("Could not create working directory %s; run %s will be skipped." % (    self.working_directory, self.name))
                self.ready = False
                return False

        # Copy the PDB file(s) in the working directory
        for pdb in self.pdbs:
            try:
                #print 'source', os.path.abspath(pdb)
                #print 'dest', self.working_directory+"/"+os.path.basename(os.path.abspath(pdb))
                safe_cp(os.path.abspath(pdb), self.working_directory+"/"+os.path.basename(os.path.abspath(pdb)), dolink=self.link_files)
            except:
                log.warning("Couldn't copy essential files for run %s; it will be skipped." % self.name)
                self.ready = False
                return False

        if self.foldx_version.rotabase:
            try:
                safe_cp(os.path.abspath(self.foldx_version.rotabase), self.working_directory+"/"+os.path.basename(os.path.abspath(self.foldx_version.rotabase)), dolink=self.link_files)
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

    def check_status(self):
        # Check if working directory exists
        if os.path.exists(self.working_directory):
            log.warning("working directory %s already exists." % self.working_directory)

            try:
                old_mutlist = self.foldx_version.parse_mutlist(os.path.join(self.working_directory, self.foldx_version.mut_list_file))
                if old_mutlist != self.mutlist:
                    return "conflicting"
            except:
                log.warning("Couldn't open mutation file, so it's not clear what's in here. Continuing anyway")


            dif_files = self.foldx_version.get_mutation_fxout_fnames(self.working_directory, self.pdbs)
            if not set(map(os.path.basename, dif_files)).issubset(os.listdir(self.working_directory)):
                return "partially_done"

            dif_file_ok = True
            try:
                for dif_file in set(map(os.path.basename, dif_files)):
                    if not self.foldx_version.check_dif_file_size(self.working_directory, dif_file, len(self.mutlist.mutations), self.runfile_processing['nruns']):
                        log.warning("File %s wasn't large as expected. The run didn't complete; it will be rerun." % dif_file)
                        dif_file_ok = False
            except:
                log.warning("Couldn't read Dif_*.fxout out file; run will be skipped")
                return "conflicting"

            pdb_basename = "".join(os.path.splitext(os.path.basename(self.pdbs[0]))[:-1])
            pdb_list_fname = self.foldx_version.mutate_pdblist_fxout_output_fname(pdb_basename)
            pdb_file_ok = True
            if not self.foldx_version.check_pdb_file_size(self.working_directory, pdb_list_fname, len(self.mutlist.mutations), self.runfile_processing['nruns']):
                log.warning("File %s wasn't large as expected. The run didn't complete; it will be rerun." % pdb_list_fname)
                pdb_file_ok = False

            if dif_file_ok and pdb_file_ok:
                return "already_done"
            else:
                return "broken"
        return "not_done"

    def process_runfile(self, **kwargs):
        pdbs = [os.path.basename(i) for i in self.pdbs]
        pdbstring = ",".join(pdbs)
        self.runfile_content = self.runfile_content.replace('$PDBS$',pdbstring)
        self.runfile_content = self.runfile_content.replace('$NRUNS$', str(kwargs["nruns"]))

    def finalize_prepare(self, **kwargs):
        self.foldx_version.save_mutlist(self.working_directory+"/"+self.foldx_version.mut_list_file, self.mutlist)

class FoldXInterfaceRun(FoldXRun):

    logfile_name = "InterfaceFoldXRun.log"

    def __init__(self, mr, pdbs=None, *args, **kwargs):

        self.mr = mr
        self.name = mr.name
        self.original_pdb = os.path.basename(mr.pdbs[0])

        if pdbs:
            self.pdbs = pdbs
        else:
            pdb_basename = "".join(os.path.splitext(self.original_pdb)[:-1])
            self.pdb_list = mr.foldx_version.mutate_pdblist_fxout_output_fname(pdb_basename)
            self.pdbs    = [os.path.join(mr.working_directory, x) for x in self.parse_pdb_list(os.path.join(mr.working_directory, self.pdb_list))]

        self.base_directory = mr.base_directory
        self.working_directory = mr.working_directory
        self.runfile_name = mr.runfile_name
        self.runfile_content = ""
        self.prepare_finalization = mr.prepare_finalization
        self.foldx_version = mr.foldx_version
        self.output_processing = {}
        self.link_files = mr.link_files
        self.write_log = mr.write_log
        self.ready = False
        self.runfile_processing = {'pdb_list' : self.pdb_list}
        self.do_clean = mr.do_clean
        self.mutlist = mr.mutlist

    def parse_pdb_list(self, pdb_list=None):
        if not pdb_list:
            pdb_list = self.pdb_list

        try:
            with open(pdb_list) as fh:
                return fh.read().splitlines()
        except:
            log.error("Couldn't parse PdbList file %s!" % pdb_list)

    def check_status(self):

        if not self.mr.finished:
            return "interface_missing_data"

        if os.path.exists(self.working_directory):
            log.warning("working directory %s already exists." % self.working_directory)
            dif_files = self.foldx_version.get_interaction_fxout_fnames(self.working_directory, self.pdbs, self)
            if not set(map(os.path.basename,dif_files[0]+dif_files[1])).issubset(set(os.listdir(self.working_directory))):
                return "not_done"
        else:
            return "not_done"

        return "already_done"

    def process_runfile(self, **kwargs):
        self.runfile_content = self.runfile_content.replace('$PDBLIST$', kwargs['pdb_list'])

