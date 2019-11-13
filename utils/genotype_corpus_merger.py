#//////////////////////////////////////////////////////////////////////////#
#                                                                          #
#   Copyright (C) 2019 by David B. Blumenthal                              #
#                                                                          #
#   This file is part of EpiGEN.                                           #
#                                                                          #
#   EpiGEN is free software: you can redistribute it and/or modify         #
#   it under the terms of the GNU General Public License as published by   #
#   the Free Software Foundation, either version 3 of the License, or      #
#   (at your option) any later version.                                    #
#                                                                          #
#   EpiGEN is distributed in the hope that it will be useful,              #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of         #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the           #
#   GNU General Public License for more details.                           #
#                                                                          #
#   You should have received a copy of the GNU General Public License      #
#   along with EpiGEN. If not, see <http://www.gnu.org/licenses/>.         #
#                                                                          #
#//////////////////////////////////////////////////////////////////////////#

"""Contains definition of GenotypeCorpusMerger class."""

import os.path
import numpy as np
import json
import bz2
import matplotlib.pyplot as plt

class GenotypeCorpusMerger(object):
    """Merges pre-computed genotype corpora.
    
    This class is employed by the script generate_genotype_corpus.py.
    Not intended for use outside of this script.
    Expects to be imported from EpiGEN's root directory.
    
    Attributes:
        corpus_ids (list of int): A list of integers that represents the IDs of the corpora that should be merged.
        pops (list of str): A list of strings representing the HAPMAP3 population codes of the corpora that should be merged.
        corpus_id (int): An integer that represents the ID of the generated corpus.
        pop (str): A string representing the HAPMAP3 population code of the merged corpus.
        genotype (numpy.array): A numpy.array with entries from range(3) that contains the merged genotypes. 
            The rows represent SNPs, the columns represent individuals.
        snps (list of (list of str)): A list with one entry for each row of self.genotype. The entries provide information 
            about the corresponding SNP.
        mafs (numpy.array): A numpy.array of floats representing the MAFs of all rows of self.genotype.
        cum_mafs (list of (float,int)): A list of pairs of the form (MAF,count) representing the cumulative MAF distribution.
        axis (int): An integer representing the axis of the merge (0 for merge along SNPs, 1 for merge along individuals). 
        num_snps (int): Number of SNPs in merged corpus.
        num_inds (int): Number of individuals in merged corpus.
        compress (bool): If True, the merged corpus is compressed. 
     """


    def __init__(self, corpus_ids, pops, corpus_id, axis, compress):
        """Initializes GenotypeCorpusGenerator.
        
        Args:
            corpus_ids (list of int): A list of integers that represents the IDs of the corpora that should be merged.
            pop (str): A list of strings representing the HAPMAP3 population codes of the corpora that should be merged.
            corpus_id (int): An integer that represents the ID of the generated corpus. 
            axis (int): An integer representing the axis of the merge (0 for merge along SNPs, 1 for merge along individuals)
            compress (bool): If True, the merged corpus is compressed.
        """
        
        # Print information.
        print("Initializing the merger ... ")
        
        self.corpus_ids = corpus_ids
        self.pops = pops
        self.pop = "MIX"
        if len(set(self.pops)) == 0:
            self.pop = self.pops[0]
        self.corpus_id = corpus_id
        self.genotype = None
        self.snps = []
        self.mafs = None
        self.cum_mafs = []
        self.axis = axis
        self.num_snps = 0
        self.num_inds = 0
        self.compress = compress
        
    
    def merge_corpora(self):
        print("Merging genotype corpora ...")
        
        # Load the first genotype corpus.
        corpus_id = self.corpus_ids[0]
        pop = self.pops[0]
        if os.path.exists("corpora/" + str(corpus_id) + "_" + pop + "_genotype.json"):
            with open("corpora/" + str(corpus_id) + "_" + pop + "_genotype.json", "rt") as jsonfile:
                self.genotype = np.asarray(json.load(jsonfile), dtype=np.uint8)
        elif os.path.exists("corpora/" + str(corpus_id) + "_" + pop + "_genotype.json.bz2"):
            with bz2.open("corpora/" + str(corpus_id) + "_" + pop + "_genotype.json.bz2", "rt", encoding="ascii") as zipfile:
                self.genotype = np.asarray(json.load(zipfile), dtype=np.uint8)
        else:
            msg = "Neither the file corpora/" + str(corpus_id) + "_" + pop + "_genotype.json.bz2 "
            msg += "nor the file corpora/" + str(corpus_id) + "_" + pop + "_genotype.json exists. "
            msg += "Change corpus or population ID or re-run generate_genotype_corpus.py." 
            raise OSError(msg)
        
        # Load SNPs.
        if os.path.exists("corpora/" + str(corpus_id) + "_" + pop + "_snps.json"):
            with open("corpora/" + str(corpus_id) + "_" + pop + "_snps.json", "rt") as jsonfile:
                self.snps = json.load(jsonfile)
        elif os.path.exists("corpora/" + str(corpus_id) + "_" + pop + "_snps.json.bz2"):
            with bz2.open("corpora/" + str(corpus_id) + "_" + pop + "_snps.json.bz2", "rt", encoding="ascii") as zipfile:
                self.snps = json.load(zipfile)
        else:
            msg = "Neither the file corpora/" + str(corpus_id) + "_" + pop + "_snps.json.bz2 "
            msg += "nor the file corpora/" + str(corpus_id) + "_" + pop + "_snps.json exists. "
            msg += "Change corpus or population ID or re-run generate_genotype_corpus.py." 
            raise OSError(msg)
        
        # Merge the corpora.
        for pos in range(1, len(self.corpus_ids)):
            corpus_id = self.corpus_ids[pos]
            pop = self.pops[0]
            if os.path.exists("corpora/" + str(corpus_id) + "_" + pop + "_genotype.json"):
                with open("corpora/" + str(corpus_id) + "_" + pop + "_genotype.json", "rt") as jsonfile:
                    try:
                        self.genotype = np.append(self.genotype, np.asarray(json.load(jsonfile), dtype=np.uint8), axis=self.axis)
                    except:
                        raise ValueError("Wrong array dimensions. Cannot merge along axis {}.".format(self.axis))
            elif os.path.exists("corpora/" + str(corpus_id) + "_" + pop + "_genotype.json.bz2"):
                with bz2.open("corpora/" + str(corpus_id) + "_" + pop + "_genotype.json.bz2", "rt", encoding="ascii") as zipfile:
                    try:
                        self.genotype = np.append(self.genotype, np.asarray(json.load(zipfile), dtype=np.uint8), axis=self.axis)
                    except:
                        raise ValueError("Wrong array dimensions. Cannot merge along axis {}.".format(self.axis))
            else:
                msg = "Neither the file corpora/" + str(corpus_id) + "_" + pop + "_genotype.json.bz2 "
                msg += "nor the file corpora/" + str(corpus_id) + "_" + pop + "_genotype.json exists. "
                msg += "Change corpus or population ID or re-run generate_genotype_corpus.py." 
                raise OSError(msg)
            if self.axis == 0:
                if os.path.exists("corpora/" + str(corpus_id) + "_" + pop + "_snps.json"):
                    with open("corpora/" + str(corpus_id) + "_" + pop + "_snps.json", "rt") as jsonfile:
                        self.snps += json.load(jsonfile)
                elif os.path.exists("corpora/" + str(corpus_id) + "_" + pop + "_snps.json.bz2"):
                    with bz2.open("corpora/" + str(corpus_id) + "_" + pop + "_snps.json.bz2", "rt", encoding="ascii") as zipfile:
                        self.snps += json.load(zipfile)
                else:
                    msg = "Neither the file corpora/" + str(corpus_id) + "_" + pop + "_snps.json.bz2 "
                    msg += "nor the file corpora/" + str(corpus_id) + "_" + pop + "_snps.json exists. "
                    msg += "Change corpus or population ID or re-run generate_genotype_corpus.py." 
                    raise OSError(msg)
                
        # Set number of SNPs and individuals.
        self.num_snps = float(np.shape(self.genotype)[0])
        self.num_inds = float(np.shape(self.genotype)[1])
        
        
    
    def compute_mafs(self):
        """Computes MAFs for all SNPs contained in self.genotype."""
        
        # Print information.
        print("Computing MAFs ... ")
        
        # Compute the minor allele frequencies.
        self.mafs = np.apply_along_axis(np.sum, 1, np.sign(self.genotype)) / self.num_inds
        
        # Compute cumulative MAF distribution.
        self.cum_mafs = []
        sorted_mafs = sorted(self.mafs.tolist())
        current_maf = sorted_mafs[0]
        counter = 1
        for pos in range(1,len(sorted_mafs)):
            if sorted_mafs[pos] != current_maf:
                self.cum_mafs.append([current_maf, counter])
                current_maf = sorted_mafs[pos]
            counter += 1
        self.cum_mafs.append([current_maf, counter])
            
        
        
    def dump_corpus(self):
        """Dumps the generated corpus to zipped JSON files."""
        
        # Print information.
        print("Serializing the genotype corpus ... ")
        
        # Dump genotype.
        if self.compress:
            with bz2.open("corpora/" + str(self.corpus_id) + "_" + self.pop + "_genotype.json.bz2", "wt", encoding="ascii") as zipfile:
                json.dump(self.genotype.tolist(), zipfile)
                
            # Dump SNPs.
            with bz2.open("corpora/" + str(self.corpus_id) + "_" + self.pop + "_snps.json.bz2", "wt", encoding="ascii") as zipfile:
                json.dump(self.snps, zipfile)
                
            # Dump MAFs.
            with bz2.open("corpora/" + str(self.corpus_id) + "_" + self.pop + "_mafs.json.bz2", "wt", encoding="ascii") as zipfile:
                json.dump(self.mafs.tolist(), zipfile)
                
            # Dump cumulative MAF distribution.
            with bz2.open("corpora/" + str(self.corpus_id) + "_" + self.pop + "_cum_mafs.json.bz2", "wt", encoding="ascii") as zipfile:
                json.dump(self.cum_mafs, zipfile)
        else:
            with open("corpora/" + str(self.corpus_id) + "_" + self.pop + "_genotype.json.bz2", "wt", encoding="ascii") as jsonfile:
                json.dump(self.genotype.tolist(), jsonfile)
                
            # Dump SNPs.
            with open("corpora/" + str(self.corpus_id) + "_" + self.pop + "_snps.json.bz2", "wt", encoding="ascii") as jsonfile:
                json.dump(self.snps, jsonfile)
                
            # Dump MAFs.
            with open("corpora/" + str(self.corpus_id) + "_" + self.pop + "_mafs.json.bz2", "wt", encoding="ascii") as jsonfile:
                json.dump(self.mafs.tolist(), jsonfile)
                
            # Dump cumulative MAF distribution.
            with open("corpora/" + str(self.corpus_id) + "_" + self.pop + "_cum_mafs.json.bz2", "wt", encoding="ascii") as jsonfile:
                json.dump(self.cum_mafs, jsonfile)
        
        # Plot cumulative MAF distribution.
        fig, ax = plt.subplots()
        plt.xlim(0,1)
        ax.plot(*zip(*self.cum_mafs))
        ax.set(xlabel="MAF", ylabel="# SNPs with MAF(SNP) <= MAF", title="Cumulative MAF Distribution")
        ax.grid()
        fig.savefig("corpora/" + str(self.corpus_id) + "_" + self.pop + "_cum_mafs.pdf")
