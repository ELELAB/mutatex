#!/usr/bin/env python

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
            # self.res_grups = [(A130, B130), (A130,B130)]
            self.res_groups = tuple([res_groups] * len(mutations.reslist))
            # self.mutations = [(L,L),(C,C)]
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

class EnergyReport:
    def __init__(self, pdbs=None):
        self.residues = {}
        self.energies = {}
        if pdbs:
            for pdb in pdbs:
                self.energies[pdb] = []
                self.residues[pdb] = []

    def add_residue(self, res, energy, pdb, do_avg=None, do_std=None, do_max=None, do_min=None):
        assert(energy.shape[0] == 1)
        energy = energy[0]

        if not pdb in self.energies.keys():
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

        dtype = [('res', 'S8')] + [(i, 'f') for i in header_cols]
        header = "\t".join(header_cols)
        #fmt= "%8s\t" + "\t".join(["%.10f" for f in header_cols])
        fmt=["%8s"] + ["%10f" for f in header_cols]

        for pdb,energies in self.energies.iteritems():

            out = [self.residues[pdb]]

            #print "ASAZ", energies,

            if do_avg:
                out.append(np.average(energies, axis=1))
            if do_std:
                out.append(np.std(energies, axis=1))
            if do_min:
                out.append(np.min(energies, axis=1))
            if do_max:
                out.append(np.max(energies, axis=1))

            out = np.array(zip(*out), dtype=dtype)

            try:
                np.savetxt(os.path.join(directory, pdb, fname), out, fmt=fmt, header=header)
            except:
                log.error("Couldn't write file %s" % fname)

class FoldXVersion:
    def __init__(self, binary=None):
        if binary:
            self.binary = os.path.abspath(binary)
        else:
            self.binary = None

    version = None
    runfile_string = None

    out_ext = "fxout"
    mut_list_file = "individual_list.txt"

    len_dif_file_header = 9

class FoldXSuiteVersion4(FoldXVersion):
    version="suite4"
    runfile_string = "-f"

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

    def parse_mutations_fxout(self, directory, pdb, this_run):

        pattern = '(\s+([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)){22}'

        pdb = pdb[0]

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

        #this_run.prepare_finalization['mutlist']

        energies = np.array(energies)
        energies = energies.reshape(len(this_run.mutlist.mutations),
                                        energies.shape[0]/len(this_run.mutlist.mutations))

        print "NRG", energies

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

        #prefixes = (self.mutation_output_pdb_WT_prefix, self.mutation_output_pdb_prefix):

        fnames = [[],[]]

        for pdb in pdbs:
            pdb_basename = "".join(os.path.splitext(os.path.basename(pdb))[:-1])
            if self.mutation_output_pdb_WT_prefix in pdb:
                fnames[0].append(os.path.join(directory, self.ac_summary_fxout_output_fname(pdb_basename)))
            else:
                fnames[1].append(os.path.join(directory, self.ac_summary_fxout_output_fname(pdb_basename)))

        return fnames

    def parse_interaction_energy_summary_fxout(self, directory, pdbs, run):
        #pattern = '\s+([A-Z])\s+([A-Z])(\s+([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)){5}'
        #pattern='\s+([A-Z])\s+([A-Z])\s+[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s+([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s+([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?){3}'

        #print "IUIU3", directory, pdb

        fnames = self.get_interaction_fxout_fnames(directory, pdbs, run, original_pdb=True)
        #print "IUIU1",     fnames[0][-1], len(fnames), len(fnames[0])
        #print "IUIU2",     fnames[1][-1], len(fnames), len(fnames[1])

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
                    #match = re.search(pattern, line.strip())
                    #if match:
                    tmp = line.strip().split()
                    idx = frozenset((tmp[1], tmp[2]))
                    if idx in energies[prefix].keys():
                        energies[prefix][idx].append(float(tmp[5]))
                    else:
                        energies[prefix][idx] = [float(tmp[5])]


                fh.close()
                print "From file %s: " %fname , energies[prefix][idx][-1], prefix, idx

            print energies, '\n-------\n'
            for k,v in energies[prefix].iteritems():
                #print energies[prefix][k], 'QQ'
                v = np.array(v)
                #print '--', v
                print directory, len(run.mutlist.mutations), v.shape, k
                energies[prefix][k] = v.reshape((len(run.mutlist.mutations),
                                      v.shape[0]/len(run.mutlist.mutations)))

        #print '====', energies



        #print "AABB", energies[types[0]].keys(), energies[types[1]].keys()
        #print "AABB", energies[types[0]].keys(), energies[types[1]].keys()
        #print "CC", energies.keys()

        interaction_groups = set(energies[types[0]].keys() + energies[types[1]].keys())

        #print interaction_groups, "IG"
        #print energies, "ENE"

        for ig in interaction_groups:
            print "GR1", energies[types[0]][ig]
            print "GR2", energies[types[1]][ig]
            print "GR1-2", energies[types[0]][ig] - energies[types[1]][ig]

            delta_energies[ig] = energies[types[1]][ig] - energies[types[0]][ig]

            #print delta_energies, "DELTAS"

        return delta_energies

    def check_dif_file_size(self, cwd, fname, nmuts, nruns):
        #if self.foldx_version(check_dif_file_size(dif_file)):
        with open("%s/%s" % (cwd,fname)) as fh:
            fsize = len(fh.readlines())

        print fsize, self.len_dif_file_header, nruns*nmuts, "YY"
        if fsize-self.len_dif_file_header == nruns*nmuts:
            return True
        return False



class FoldXRun(object):

    ready = False

    finished = False

    runfile_name = "runfile.txt"

    logfile_name = "FoldXrun.log"

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
                write_log=False,
                clean=False,
                *args,
                **kwargs):

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
        if output_processing is None:
            self.output_processing = {}
        else:
            self.output_processing = output_processing
        self.rotabase = rotabase
        self.link_files = link_files
        self.write_log = write_log
        self.ready = False
        self.do_clean = clean


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
            log.warning("Could not create working directory %s; run %s will be skipped." % (    self.working_directory, self.name))
            self.ready = False
            return False

        # Copy the PDB file(s) in the working directory
        for pdb in self.pdbs:
            try:
                #print 'source', os.path.abspath(pdb)
                #print 'dest', self.working_directory+"/"+os.path.basename(os.path.abspath(pdb))
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
            print "listdir", os.listdir(self.working_directory)
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
        self.partial_clean()
        if self.do_clean:
            self.clean()

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
            if os.path.isfile(fname) and not fname.endswith(".fxout") and not fname == os.path.basename(fname):
                log.info("removing %s" % fname)
                os.remove(fname)

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


    def check_status(self):
        # Check if working directory exists
        if os.path.exists(self.working_directory):
            log.warning("working directory %s already exists." % self.working_directory)

            try:
                old_mutlist = self.foldx_version.parse_mutlist(self.working_directory+"/"+self.foldx_version.mut_list_file)
            except:
                log.warning("Couldn't open mutation file, so it's not clear what's in here. Run will be repeated.")
                return "not_done"

            if old_mutlist != self.mutlist:
                return "conflicting"

            dif_files = self.foldx_version.get_mutation_fxout_fnames(self.working_directory, self.pdbs)
            if not set(map(os.path.basename, dif_files)).issubset(os.listdir(self.working_directory)):
                return "not_done"
            #dif_files = [f for f in os.listdir(self.working_directory) if f.startswith(self.foldx_version.dif_fxout_prefix) and f.endswith(self.foldx_version.out_ext)]

            #if len(dif_files) > 1:
            #    log.warning("More than one Dif_*.fxout files present in dir! Something fishy is going on. Run will be skipped.")
            #    return "conflicting"
            #elif len(dif_files) == 0:
            try:
                #print "diffile", dif_file
                for dif_file in set(map(os.path.basename, dif_files)):
                    print len(self.mutlist.mutations), self.runfile_processing['nruns'], "CHECK"
                    if self.foldx_version.check_dif_file_size(self.working_directory, dif_file, len(self.mutlist.mutations), self.runfile_processing['nruns']):
                        return "already_done"
                    else:
                        log.warning("File %s wasn't large as expected. The run didn't complete; it will be rerun." % dif_file)
                        return "not_done"
            except:
                log.warning("Couldn't read Dif_*.fxout out file; run will be skipped")
                return "conflicting"
        else:
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
        # super(FoldXInterfaceRun, self).__init__(*args, **kwargs)

        self.name = mr.name
        self.original_pdb = os.path.basename(mr.pdbs[0])

        if pdbs:
            self.pdbs = pdbs
        else:
            pdb_basename = "".join(os.path.splitext(self.original_pdb)[:-1])
            self.pdb_list = mr.foldx_version.mutate_pdblist_fxout_output_fname(pdb_basename)
            self.pdbs    = map(lambda x: os.path.join(mr.working_directory, x), self.parse_pdb_list(os.path.join(mr.working_directory, self.pdb_list)))

        self.base_directory = mr.base_directory
        self.working_directory = mr.working_directory
        self.runfile_name = mr.runfile_name
        self.runfile_content = ""
        self.prepare_finalization = mr.prepare_finalization
        self.foldx_binary = mr.foldx_binary
        self.foldx_version = mr.foldx_version
        self.output_processing = {}
        self.rotabase = mr.rotabase
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
                #print "SS",fh.read().splitlines()
                return fh.read().splitlines()
        except:
            log.error("Couldn't parse PdbList file %s!" % pdb_list)

    def check_status(self):
        # Check if working directory exists
        print "WD", self.working_directory
        if os.path.exists(self.working_directory):
            log.warning("working directory %s already exists." % self.working_directory)
            dif_files = self.foldx_version.get_interaction_fxout_fnames(self.working_directory, self.pdbs, self)
            #print "diffiles", dif_files
            #print dif_files
            #print "SETO", set(map(os.path.basename,dif_files[0]+dif_files[1]))
            #print os.listdir(self.working_directory)
            #print "TUTO", set(os.listdir(self.working_directory))
            #print dif_files[0]+dif_files[1]
            #print "INTERSECTO", set(map(os.path.basename,dif_files[0]+dif_files[1])).difference(set(os.listdir(self.working_directory)))
            if not set(map(os.path.basename,dif_files[0]+dif_files[1])).issubset(set(os.listdir(self.working_directory))):
                #print "NODONE1"
                return "not_done"
        else:
            #print "NODONE2"
            return "not_done"
        #print "DONEE"
        return "already_done"

    def process_runfile(self, **kwargs):
        print kwargs, 'PP'
        self.runfile_content = self.runfile_content.replace('$PDBLIST$', kwargs['pdb_list'])

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

    if os.path.abspath(source) == os.path.abspath(destination):
        return

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

def foldx_worker(run):
    log.info("starting FoldX run %s" % run.name)
    return (run.name, run.run())

def parallel_foldx_run(foldx_runs, np):

    pool = mp.Pool(np)

    print foldx_runs

    result = pool.imap_unordered(foldx_worker, foldx_runs)

    pool.close()
    pool.join()

    print result

    return list(result)


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

def get_foldx_sequence(pdb, multimers=True):
    pdb_parser = PDB.PDBParser()

    log.info("loading pdb file %s" % pdb)

    try:
        structure = pdb_parser.get_structure('structure',pdb)
    except:
        log.error("Couldn't open file %s." % pdb)
        exit(1)

    residue_list = []
    sequences = {}
    for model in structure:
        for chain in model:
            chain_name = chain.get_id()
            sequences[chain_name] = ''
            for residue in chain:
                try:
                    res_code = PDB.Polypeptide.three_to_one(residue.get_resname())
                except:
                    log.warning("Residue %s in file %s couldn't be recognized; it will be skipped" %(residue, pdb))
                    continue
                if not multimers:
                    residue_list.append(("%s%s%d") % (res_code, chain.get_id(), residue.get_id()[1]))
                else:
                    sequences[chain_name] += res_code

    print "SEQUENCES", sequences
    if multimers:
        print "MERING"
        collated_chains = []
        seq_ids, seqs = zip(*list(sequences.iteritems()))
        seq_ids = np.array(seq_ids)
        unique_seqs, unique_idxs = np.unique(seqs, return_inverse=True)

        for i in np.unique(unique_idxs):
            collated_chains.append(seq_ids[unique_idxs == i])
        print "COLLATED", collated_chains

        for cg in collated_chains:
            for model in structure:
                for residue in model[cg[0]]:
                    resid = residue.get_id()[1]
                    try:
                        res_code = PDB.Polypeptide.three_to_one(residue.get_resname())
                    except:
                        log.warning("Residue %s in file %s couldn't be recognized; it will be skipped" %(residue, pdb))
                        continue
                    residue_list.append(tuple([ "%s%s%d" % (res_code, c, resid) for c in cg ]))

    print "AA", residue_list
    return tuple(residue_list)

def save_energy_file(fname, data, fmt="%.5f", do_avg=True, do_std=False, do_min=False, do_max=False, axis=1):

    out = []
    header_cols = []

    if do_avg:
        out.append(np.average(data, axis=axis))
        header_cols.append("avg")
    if do_std:
        out.append(np.std(data, axis=axis))
        header_cols.append("std")
    if do_min:
        out.append(np.min(data, axis=axis))
        header_cols.append("min")
    if do_max:
        out.append(np.max(data, axis=axis))
        header_cols.append("max")

    header = "\t".join(header_cols)

    out = np.array(out).T

    try:
        np.savetxt(fname, out, fmt=fmt, header=header)
    except:
        log.error("Couldn't write file %s" % fname)

def save_interaction_energy_file(fname, data, fmt="%.5f", do_avg=True, do_std=False, do_min=False, do_max=False, axis=1):

    out = []
    header_cols = []

    if do_avg:
        out.append(np.average(data, axis=axis))
        header_cols.append("avg")
    if do_std:
        out.append(np.std(data, axis=axis))
        header_cols.append("std")
    if do_min:
        out.append(np.min(data, axis=axis))
        header_cols.append("min")
    if do_max:
        out.append(np.max(data, axis=axis))
        header_cols.append("max")

    header = "\t".join(header_cols)

    out = np.array(out).T

    if True: #try:
        np.savetxt(fname, out, fmt=fmt, header=header)
    else: #except:
        log.error("Couldn't write file %s" % fname)


# Gather info about FoldX version

mutatex_version = "0.6"

name_separator = "_"

def main():

    foldx_versions = [FoldXSuiteVersion4]

    supported_foldx_versions = {v.version:v for v in foldx_versions}

    foldx_binary_var = os.getenv('FOLDX_BINARY')
    foldx_rotabase_var = os.getenv('FOLDX_ROTABASE')

    if foldx_binary_var is None:
        foldx_binary_var = ""
    if foldx_rotabase_var is None:
        foldx_rotabase_var = ""

    parser = argparse.ArgumentParser(description='Setup and run in silico saturation mutagenesis with FoldX. Release %s' % mutatex_version)

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
    parser.add_argument('--self-mutate', dest="selfmutate", default=False, action="store_true", help="Ignore mutation list and perform mutation to the same residue")
    parser.add_argument('-x', '--foldx-binary', dest="foldx_binary", action='store', type=str, help="Location of the FoldX binary (default: content of the FOLDX_BINARY system variable", default=foldx_binary_var)
    parser.add_argument('--rotabase', dest="rotabase", action="store", type=str, help="Location of the FoldX rotabase.txt file (default: content of the FOLDX_ROTABASE system variable", default=foldx_rotabase_var)
    parser.add_argument('--foldx-version', dest="foldx_version", action='store', choices=supported_foldx_versions.keys(), default=supported_foldx_versions.keys()[0], help="FoldX version to be used (possible options: %s" % ", ".join(supported_foldx_versions.keys()))
    parser.add_argument('-v', '--verbose', dest="verbose", default=False, action='store_true', help="Verbose mode")
    parser.add_argument('--foldx-log', dest="write_log", default=False, action='store_true', help="Write FoldX standard output on file for each run")
    parser.add_argument('-l','--use-links', dest="use_links", default=False, action='store_true', help="Use links instead of copying files as much as possibile")
    parser.add_argument('--repair-runfile-template','--repair', dest="repair_runfile_template", type=str, default="repair_runfile_template.txt", help="Template runfile for repair runs (default: ./repair_runfile_template.txt)")
    parser.add_argument('--mutate-runfile-template','--mutate', dest="mutate_runfile_template", type=str, default="mutate_runfile_template.txt", help="Template runfile for mutation runs (default: ./mutate_runfile_template.txt)")
    parser.add_argument('--interface-runfile-template','--interface', dest="interface_runfile_template", type=str, default="interface_runfile_template.txt", help="Template runfile for mutation runs (default: ./interface_runfile_template.txt)")
    parser.add_argument('--binding-interface', dest="interface", action='store_true', default=False, help="Do calculate binding DDG with mutations")
    parser.add_argument('--clean', dest="clean", action='store_true', default=False, help="Clean output directories after calculation")
    parser.add_argument('-c','--multimers', dest='multimers', action='store_false', default=True, help="Whether to consider multimers")

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

    log.info("FoldX version %s will be used" % current_version.version)

# defaults
    repair_dirname = "repair"
    mutations_dirname = "mutations"
    results_dirname = "results"
    mutations_results_dirname = "mutation_ddgs"
    interface_results_dirname = "interface_ddgs"
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
                                        clean = args.clean
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
    if args.selfmutate:
        log.info("Mutation to self will be performed; mutlist will be ignored")
        mutation_reslist = None
    elif args.mutlist:
        log.info("File %s will be used as mutation list" % args.mutlist)
        try: # XXX to be changed
            mutation_reslist = ResList(fname=args.mutlist)
        except:
            log.error("Couldn't parse mutation list")
            exit(1)
    else:
        log.info("default mutation list will be used.")
        mutation_reslist = ResList(reslist=default_mutlist)

    log.info("mutation list is: %s" % mutation_reslist)

# Select PDB(s) to be used
    #if args.skip_repair:
        #repaired_pdbs_list = args.pdb
    #else:
    repaired_pdbs_list = [repair_runs[i].working_directory+"/"+current_version.repair_pdb_output_fname(pdbs_list[i]) for i in range(len(pdbs_list))]
    log.info("list of PDBs to be used: %s" % ", ".join(repaired_pdbs_list))

# create working directory
    if args.selfmutate:
        mutations_dirname = "selfmutations"

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
            log.error("PDB %s will be ignored." % pdb)
            repaired_pdbs_list.remove(pdb)
            continue
        residues_list.append(get_foldx_sequence(pdb, multimers=args.multimers))
        print residues_list

    print residues_list,"ZZS"
    unique_residues = tuple(set(residues_list))
    if len(unique_residues) != 1:
        log.error("The supplied PDB files must have identical sequences. Exiting...")
        exit(1)

    unique_residues = list(unique_residues[0])

    str_unique_residues = ""
    for res in unique_residues:
        str_unique_residues += "(%s) " % ",".join(res)
    log.info("The following residue groups will be considered: %s" % str_unique_residues)
    print unique_residues

    mutation_runs = []

    if args.interface:
        mutate_clean = False
    else:
        mutate_clean = args.clean

    for r in unique_residues:
        print "errR", r
        print "mres", mutation_reslist

        mutlist = MutationList(r, mutation_reslist, selfmutate=args.selfmutate)

        print mutlist.res_groups
        print mutlist.mutations
        name = name_separator.join(r)
        for pdb in repaired_pdbs_list:
            this_workdir = working_directory+"/"+os.path.splitext(os.path.basename(pdb))[0]
            mutation_runs.append(FoldXMutateRun(name = name,
                                        foldx_binary = current_version.binary,
                                        foldx_version = current_version,
                                        base_directory = this_workdir,
                                        pdbs = [pdb],
                                        runfile_content = load_runfile(args.mutate_runfile_template),
                                        rotabase = args.rotabase,
                                        link_files = args.use_links,
                                        write_log = args.write_log,
                                        runfile_processing = {"nruns":args.nruns},
                                        mutlist=mutlist,
                                        clean = mutate_clean
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

    #working_directory = os.path.join(main_dir, results_dirname)
    log.info("Working directory is: %s" % working_directory)

    if args.interface:
        log.debug("Calculating interfaces")
        interface_runs = []
        for mr in mutation_runs:
            #this_pdbs = zip(*current_version.get_mutation_pdb_fnames(mr.working_directory, mr.pdbs, mr, include_original=True))
            #print "pdbs:", this_pdbs
            #for i,pdbs in enumerate(this_pdbs):
                #original_pdb = pdbs[-1]
                #pdbs = pdbs[:-1]
            this_interaction_run = FoldXInterfaceRun(mr)
            this_interaction_run.do_clean = args.clean
            this_interaction_run.runfile_content = load_runfile(args.interface_runfile_template)
            this_interaction_run.runfile_name = "interaction_energy_runfile.txt"
            this_interaction_run.prepare()
            interface_runs.append(this_interaction_run)
        log.info("Running interface runs")
        interface_outcome = parallel_foldx_run(interface_runs, np=args.np)
        if False in zip(*interface_outcome)[1]:
            log.warning("Some interface runs failed to complete.")

    if not args.skip_report:

        if args.selfmutate:
            report = EnergyReport()

        energies = []

        working_directory = os.path.join(main_dir, results_dirname, mutations_results_dirname)
        safe_makedirs(working_directory)

        for pdb in repaired_pdbs_list:
            this_pdb_dir = working_directory+"/"+os.path.splitext(os.path.basename(pdb))[0]
            safe_makedirs(this_pdb_dir)
        safe_makedirs(working_directory+"/"+averages_dirname)

        for res in unique_residues:
            name = name_separator.join(res)
            energies = []

            #print "UQ", unique_residues
            #for i in mutation_runs: print i.name
            this_runs = filter(lambda x: x.name == name, mutation_runs)
            #print this_runs, "THIS_RUNS"
            dobreak = False
            for r in this_runs:
                #print "NONO", r.name
                for pdb in r.pdbs:
                    try:
                        energies.append(current_version.parse_mutations_fxout(r.working_directory, [os.path.basename(pdb)], r))
                    except:
                        log.warning("Couldn't parse energy file for PDB %s; mutation site %s will be skipped." % (pdb, r.name))
                        dobreak=True
			break
                    if args.selfmutate:
                        report.add_residue(pdb="".join(os.path.splitext(os.path.basename(pdb))[:-1]), res=r.name, energy=energies[-1])
                    else:
                        save_energy_file(working_directory+"/"+"".join(os.path.splitext(os.path.basename(pdb))[:-1])+"/"+r.name, energies[-1], do_avg=True, do_std=True, do_max=True, do_min=True)
            if dobreak:
                continue
            if not args.selfmutate:
                print "SHAPPO", np.array(energies).shape
                save_energy_file(working_directory+"/"+averages_dirname+"/"+r.name, np.average(energies, axis=2), axis=0, do_avg=True, do_std=True, do_max=True, do_min=True)
            else:
                report.save(working_directory)

        if args.interface:
            working_directory = os.path.join(main_dir, results_dirname, interface_results_dirname)

            for pdb in repaired_pdbs_list:
                this_pdb_dir = os.path.join(working_directory, os.path.splitext(os.path.basename(pdb))[0])
                safe_makedirs(this_pdb_dir)
            safe_makedirs(os.path.join(working_directory, averages_dirname))

            unique_residues_str = [ "_".join(i) for i in unique_residues ]
            for res in unique_residues_str:
                energies = {}
                this_runs = []
                tmp_this_runs = filter(lambda x: x.name == res, interface_runs)
                original_pdbs = list(set([r.original_pdb for r in tmp_this_runs]))
                for op in original_pdbs:
                    this_runs.append(next((r for r in tmp_this_runs if r.original_pdb == op)))

                #print "RUNNY", [r.name for r in this_runs]
                for run in this_runs:
                                        #print "RUNNN", run, run.name, run.pdbs
                    try:
                        this_energies = run.foldx_version.parse_interaction_energy_summary_fxout(    run.working_directory,
                                                                                                    run.pdbs,
                                                                                                    run    )
                    except:
                        log.warning("Couldn't parse energy file for PDB %s; mutation site %s will be skipped." % (pdb, r.name))
                        continue

                    for k in this_energies.keys():
                        labels = tuple(sorted(list(k)))
                        this_wd = os.path.join(    working_directory,
                                                "".join(os.path.splitext(os.path.basename(run.original_pdb))[:-1]),
                                                "%s-%s" % labels)
                        safe_makedirs(this_wd)
                        #print "KEYIN", k
                        if k not in energies.keys():
                            energies[k] = [this_energies[k]]
                            #print "NOVEL"
                        else:
                            #print "OLDEL"
                            energies[k].append(this_energies[k])

                                                #print "SALVO", "%s/%s/%s_%s-%s" % (working_directory,"".join(os.path.splitext(os.path.basename(run.original_pdb))[:-1]),run.name,labels[0],labels[1])
                        save_interaction_energy_file(os.path.join(this_wd, run.name),
                                                        this_energies[k],
                                                        axis=1,
                                                        do_avg=True,
                                                        do_std=True,
                                                        do_max=True,
                                                        do_min=True)

# working_directory+"/"+averages_dirname+"/"+r.name, np.average(energies, axis=0), axis=1, do_avg=True, do_std=True, do_max=True, do_min=True)

                for k,e in energies.iteritems():
                    labels = tuple(sorted(list(k)))
                    this_wd = os.path.join(working_directory, averages_dirname,
                                            "%s-%s" % labels)
                    safe_makedirs(this_wd)

                    save_interaction_energy_file(os.path.join(this_wd, run.name),
                            np.average(energies[k],axis=2),
                            axis=0,
                            do_avg=True,
                            do_std=True,
                            do_max=True,
                            do_min=True)


            # for run in interface_runs:
            #     #print "=== run %s" % run
            #     # def get_interaction_fxout_names(self, directory, pdbs, run):
            #     energies = run.foldx_version.parse_interaction_energy_summary_fxout(    run.working_directory,
            #                                                                             [os.path.basename(run.original_pdb)],
            #                                                                             run
            #                                                                         )
            #     for k,v in energies.iteritems():
            #         k = sorted(list(k))
            #         working_directory+"/"+"".join(os.path.splitext(os.path.basename(pdb))[:-1])+"/"+r.name
            #         save_interaction_energy_file("%s/%s/%s_%s-%s" % (working_directory,"".join(os.path.splitext(os.path.basename(pdb))[:-1]),run.name,k[0],k[1]),
            #                                             v,
            #                                             axis=1,
            #                                             do_avg=True,
            #                                             do_std=True,
            #                                             do_max=True,
            #                                             do_min=True)

    else:
        log.info("Reporting phase was skipped, as requested.")

    log.info("All done!")

if __name__ == '__main__':
    main()
