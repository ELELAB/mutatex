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

def parse_ddg_file(fname, reslist, full=False):
    fh = open(fname, 'r')
    
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
            log.error("file %s has %d values, with %d required. Exiting..." % (fname, len(ddgs), len(reslist)))
            exit(1)

    if full:
        return ddgs, stds, mins, maxs

    return ddgs

def parse_mutlist_file(fname):
    fh = open(fname, 'r')

    restypes = []
    
    for line in fh:
        if line and not line.startswith("#"):
            restypes.append(line.strip()[0])
    fh.close()
    
    return restypes

def get_residue_list(structure, multimers=True):

    models = list(structure.get_models())

    if len(models) > 1:
        log.warning("%d models are present in the input PDB file; only the first will be used." % len(models))
    if len(models) < 1:
        log.error("The input PDB file does not contain any model. Exiting ...")
        exit(1)

    model = models[0]

    residue_list = []
    sequences = {}

    for chain in model:
        chain_name = chain.get_id()
        sequences[chain_name] = ''
        for residue in chain:
            if True:
                res_code = PDB.Polypeptide.three_to_one(residue.get_resname())
            else:
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
    return residue_list

def get_foldx_sequence(model, multimers=True):

    sequences = {}
    residue_list = []
    
    for chain in model:
        chain_name = chain.get_id()
        sequences[chain_name] = ''
        for residue in chain:
            try:
                res_code = PDB.Polypeptide.three_to_one(residue.get_resname())
            except:
                log.warning("Residue %s couldn't be recognized; it will be skipped" % residue)
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

    return tuple(residue_list)

"""
def seq_to_res(seq):
    if options.multimers:
        return "_".join(seq)
    else:
        return seq
"""
