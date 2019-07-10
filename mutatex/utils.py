#!/usr/bin/env python

#    utils.py: miscellaneous utilities for MutateX plotting scripts
#    Copyright (C) 2015, Matteo Tiberti <matteo.tiberti@gmail.com> 
#                        Thilde Bagger Terkelsen <ThildeBT@gmail.com>
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

import numpy as np
import logging as log
from six import iteritems
from Bio import PDB
import sys
import os
import argparse
import multiprocessing as mp
import logging as log
import shutil
import re
import numpy as np
import tarfile as tar
import platform

def init_arguments(arguments, parser):

    assert len(set(arguments)) == len(arguments)

    for arg in arguments:
        if   arg == 'pdb':
            parser.add_argument("-p","--pdb", dest="in_pdb", help="Input PDB file")
        elif arg == 'data':
            parser.add_argument("-d","--data-directory", dest="ddg_dir", type=str, help="Input DDG data directory")
        elif arg == 'mutation_list':
            parser.add_argument("-l","--mutation_list", dest="mutation_list",  help="MutateX mutation list file")
        elif arg == 'multimers':
            parser.add_argument("-M","--multimers", dest="multimers", default=True, action='store_false', help="Do not use multimers (default: yes)")
        elif arg == 'labels':
            parser.add_argument("-b","--label-list", dest="labels", help="Residue label list")
        elif arg == 'fontsize':
            parser.add_argument("-f","--fontsize",dest='fontsize',action='store', type=int, default=8, help="Axis label font size")
        elif arg == 'verbose':
            parser.add_argument("-v","--verbose", dest="verbose", action="store_true", default=False, help="Toggle verbose mode")
        elif arg == 'title':
            parser.add_argument("-i","--title", dest='title', type=str, default=None, help="Title for the output image file")
        elif arg == 'color':
            parser.add_argument("-c","--color", dest='mycolor', type=str, default="black", help="Color used for plotting")
        elif arg == 'splice':
            parser.add_argument("-s","--splice",dest='sv',action='store', type=int, default=20, help="Divide data in multiple plots, use -s residues per plot")
        else:
            raise NameError

    return parser

def parse_ddg_file(fname, reslist, full=False):
    try:
        fh = open(fname, 'r')
    except:
        log.error("Couldn't open energy file %s" % fname)
        raise IOError

    ddgs = []

    if full:
        stds = []
        mins = []
        maxs = []

    for line in fh:
        if line and not line.startswith("#"):
            tmp = [ float(i) for i in line.split() ]
            ddgs.append(tmp[0])
            if full:
                stds.append(tmp[1])
                mins.append(tmp[2])
                maxs.append(tmp[3])

    if reslist is not None:
        if len(ddgs) != len(reslist):
            log.error("file %s has %d values, with %d required." % (fname, len(ddgs), len(reslist)))
            raise TypeError

    fh.close()

    if full:
        return ddgs, stds, mins, maxs

    return ddgs

def parse_mutlist_file(fname):

    try:
        fh = open(fname, 'r')
    except IOError:
        log.error("Couldn't open mutation list file %s" % fname)
        raise IOError

    restypes = []
    
    for line in fh:
        if line and not line.startswith("#"):
            str_line = line.strip()
            if len(str_line) == 0:
                continue
            if len(str_line) > 1:
                log.warning("more than one character per line found in mutation list file; only the first letter will be considered")
            restypes.append(line.strip()[0])

    fh.close()

    if len(set(restypes)) != len(restypes):
        log.error("mutation list file contains duplicates")
        raise TypeError

    if len(restypes) == 0:
        log.error("No residue types found in mutation list")
    
    return restypes

def get_residue_list(infile, multimers=True, get_structure=False):

    parser = PDB.PDBParser()

    try:
        structure = parser.get_structure("structure", infile)
    except IOError:
        log.error("couldn't read or parse your PDB file")
        raise IOError

    models = list(structure.get_models())

    if len(models) > 1:
        log.warning("%d models are present in the input PDB file; only the first will be used." % len(models))
    if len(models) < 1:
        log.error("the input PDB file does not contain any model. Exiting ...")
        raise IOError

    model = models[0]

    residue_list = []
    sequences = {}

    for chain in model:
        chain_name = chain.get_id()
        sequences[chain_name] = ''
        for residue in chain:
            try:
                res_code = PDB.Polypeptide.three_to_one(residue.get_resname())
            except:
                log.warning("Residue %s couldn't be recognized; it will be skipped" % residue )
                continue
            if not multimers:
                residue_list.append(("%s%s%d") % (res_code, chain.get_id(), residue.get_id()[1]))
            else:
                sequences[chain_name] += res_code

    if multimers:
        collated_chains = []
        seq_ids, seqs = list(zip(*list(iteritems(sequences))))
        seq_ids = np.array(seq_ids)
        unique_seqs, unique_idxs = np.unique(seqs, return_inverse=True)

        for i in np.unique(unique_idxs):
            collated_chains.append(seq_ids[unique_idxs == i])

        for cg in collated_chains:
            for model in structure:
                for residue in model[cg[0]]:
                    resid = residue.get_id()[1]
                    try:
                        res_code = PDB.Polypeptide.three_to_one(residue.get_resname())
                    except:
                        log.warning("Residue %s couldn't be recognized; it will be skipped" % residue)
                        continue
                    residue_list.append(tuple([ "%s%s%d" % (res_code, c, resid) for c in cg ]))

    if get_structure:
        return residue_list, structure
    return residue_list

#########################################################

def get_foldx_sequence(pdb, multimers=True):

    parser = PDB.PDBParser()
    try:
        structure = parser.get_structure("structure", infile)
    except:
        log.error("couldn't read or parse your PDB file")
        raise IOError

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
                    residue_list.append(tuple(["%s%s%d" % (res_code, chain.get_id(), residue.get_id()[1])]))
                else:
                    sequences[chain_name] += res_code

    if multimers:
        collated_chains = []
        seq_ids, seqs = list(zip(*list(iteritems(sequences))))
        seq_ids = np.array(seq_ids)
        unique_seqs, unique_idxs = np.unique(seqs, return_inverse=True)

        for i in np.unique(unique_idxs):
            collated_chains.append(seq_ids[unique_idxs == i])

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

    return tuple(residue_list)


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

    parser = PDB.PDBParser()
    try:
        structure = parser.get_structure("structure", infile)
    except:
        log.error("couldn't read or parse your PDB file")
        raise IOError

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
        log.error("Couldn't open runfile %s." % runfile)
        raise IOError

    return data

def foldx_worker(run):
    log.info("starting FoldX run %s" % run.name)
    return (run.name, run.run())

def parallel_foldx_run(foldx_runs, np):

    pool = mp.Pool(np)

    result = pool.imap_unordered(foldx_worker, foldx_runs)

    pool.close()
    pool.join()

    return list(result)


def split_pdb(filename, checked):

    pdb_list = []
    
    writer = PDB.PDBIO()
    parser = PDB.PDBParser()

    try:
        structure = parser.get_structure("structure", infile)
    except:
        log.error("couldn't read or parse your PDB file")
        raise IOError

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
        log.error("Couldn't write energy file %s" % fname)
        raise IOError

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

    try:
        np.savetxt(fname, out, fmt=fmt, header=header)
    except:
        log.error("Couldn't interaction energy write file %s" % fname)
        raise IOError

def compress_mutations_dir(cwd, mutations_dirname, mutations_archive_fname='mutations.tar.gz'):

    archive_path = os.path.join(cwd, mutations_archive_fname)
    mutations_dir_path = mutations_dirname

    log.info("Compressing mutations directory as per user request")
    if not os.path.isdir(cwd):
        log.warning("Directory mutations doesn't exist; it won't be compressed.")

    try:
        fh = tar.open(archive_path, 'w:gz')
    except:
        log.warning("Couldn't open compressed file %s for writing." % mutations_archive_fname)
        return

    try:
        fh.add(mutations_dir_path)
    except:
        log.warning("Couldn't build compressed archive. This step will be skipped.")
        fh.close()
        os.remove(archive_path)
        return

    fh.close()
    log.info("Removing mutations directory ...")
    shutil.rmtree(mutations_dir_path)
    return



"""
def seq_to_res(seq):
    if options.multimers:
        return "_".join(seq)
    else:
        return seq
"""
