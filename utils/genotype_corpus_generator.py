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

"""Contains definition of GenotypeCorpusGenerator class."""

import platform
import subprocess
import os.path
import csv
import numpy as np
import json
import bz2
import matplotlib.pyplot as plt

class GenotypeCorpusGenerator(object):
    """Generates genotype corpus by calling HAPGEN2 and merging the obtained chromosome-wise genotypes.
    
    This class is employed by the script generate_genotype_corpus.py.
    Not intended for use outside of this script.
    Expects to be imported from EpiGEN's root directory.
    
    Attributes:
        hapmap_dir (str): The path to the HAPMAP3 directory.
        hapgen_dir (str): The path to the directory that contains the HAPGEN2 binary.
        chroms (list of int): A list of integers between 1 and 22 representing the chromosomes 
            for which HAPGEN2 generates the genotypes.
        num_inds (int): An integer that represents the number of individuals for which 
            genotypes are constructed.
        num_snps (int): An integer that represents the number of SNPs in the generated corpus.
        corpus_id (int): An integer that represents the ID of the generated corpus.
        pop (str): A string representing the HAPMAP3 population.
        dummy_disease_snps (dict of (int,int)): A dict that, for each chromosome, specifies the position 
            of a SNP that can be passed to HAPGEN2 as argument of the -dl option.
        genotype (numpy.array): A numpy.array with entries from range(3) that contains the merged genotypes. 
            The rows represent SNPs, the columns represent individuals.
        snps (list of (list of str)): A list with one entry for each row of self.genotype. The entries provide information 
            about the corresponding SNP.
        mafs (numpy.array): A numpy.array of floats representing the MAFs of all rows of self.genotype.
        cum_mafs (list of (float,int)): A list of pairs of the form (MAF,count) representing the cumulative MAF distribution.
        compress (bool): If True, the generated corpus is compressed.
     """


    def __init__(self, chroms, num_inds, corpus_id, pop, compress):
        """Initializes GenotypeCorpusGenerator.
        
        Args:
            chroms (list of int): A list of integers between 1 and 22 representing the chroms for which
                HAPGEN2 should be called. Duplicates are ignored.
            num_inds (int): An integer representing the number of individuals for which genotypes should
                be constructed.
            corpus_id (int): An integer that represents the ID of the generated corpus.
            compress (bool): If True, the generated corpus is compressed.
        """
        
        # Print information.
        print("Initializing the generator ... ")
        
        # Set attributes.
        self.hapmap_dir = "ext/HAPMAP3/"
        self.hapgen_dir = "ext/HAPGEN2/"
        if platform.system() == "Linux":
            self.hapgen_dir += "Linux/"
        elif platform.system() == "Darwin":
            v, _, _ = platform.mac_ver()
            v = float('.'.join(v.split('.')[:2]))
            if v > 10.14:
                raise Exception("Unsupported macOS version {}. Only macOS < 10.15 is supported. Use pre-computed corpora instead.".format(v))
            else:
                self.hapgen_dir += "Darwin/"
        else:
            raise Exception("Unsupported system {}. Only Linux and macOS < 10.15 are supported. Use pre-computed corpora instead.".format(platform.system()))
        self.chromosomes = set(chroms)
        self.num_inds = num_inds
        self.num_snps = 0
        self.corpus_id = corpus_id
        self.pop = pop
        self.genotype = None
        self.snps = []
        self.mafs = None
        self.cum_mafs = []
        self.compress = compress
        
    def call_hapgen2(self):
        """Calls HAPGEN2 to generate genotypes for all chromosomes."""
        
        # Initialize list of child processes (one for each chromosome).
        child_processes = []
        
        # Decompress HAPMAP3 data unless it is already decompressed.
        print("Decompressing HAPMAP3 data ...")
        for chrom in self.chromosomes:
            gen_map = self.hapmap_dir + "genetic_map_chr" + str(chrom) + "_combined_b36.txt.bz2"
            if not os.path.exists(gen_map):
                gen_map = ""
            hap = self.hapmap_dir + self.pop + ".chr" +str(chrom) + ".hap.bz2"
            if not os.path.exists(hap):
                hap = ""
            legend = self.hapmap_dir + "hapmap3.r2.b36.chr" + str(chrom) + ".legend.bz2"
            if not os.path.exists(legend):
                legend = ""
            if (gen_map != "") or (hap != "") or (legend != ""):
                p = subprocess.Popen("bunzip2 " + gen_map + " " + hap + " " + legend, shell=True)
                child_processes.append(p)
            
        # Wait for all data to be decompressed.
        for p in child_processes:
            p.wait()
        
        # Print information.
        print("Calling HAPGEN2 ... ")
        p = []
        
        # Call HAPGEN2 for each chromosome.
        for chrom in self.chromosomes:
            
            # Get filenames of HAPMAP3 data.
            gen_map = self.hapmap_dir + "genetic_map_chr" + str(chrom) + "_combined_b36.txt"
            hap = self.hapmap_dir + self.pop + ".chr" +str(chrom) + ".hap"
            legend = self.hapmap_dir + "hapmap3.r2.b36.chr" + str(chrom) + ".legend"
            
            # Build command string to call HAPGEN2. Since we are only interested in the genotype,
            # we set up HAPGEN2 to generate a null-model without cases.
            command = "cd " + self.hapgen_dir + "; ./hapgen2"
            command += " -m ../../../" + gen_map
            command += " -h ../../../" + hap
            command += " -n " + str(self.num_inds) + " 0"
            command += " -l ../../../" + legend
            
            # With the -dl option of HAPGEN2, the user can specify the disease SNPs.
            # HAPGEN2 requires this option, even if only controls are generated,
            # although in this case -dl has no effect. Moreover, the selected SNP
            # has to have a non-zero MAF. To ensure that HAPGEN2 is happy, we here select
            # a dummy disease SNP that meets these requirements.
            dummy_snp = 0
            with open(hap, "r") as csvfile:
                reader = csv.reader(csvfile, delimiter=" ")
                found_dummy_snp = False
                while not found_dummy_snp:
                    dummy_snp += 1
                    found_zero = False
                    found_one = False
                    for gen in next(reader):
                        if gen == "0":
                            found_zero = True
                        elif gen == "1":
                            found_one = True
                        if found_zero and found_one:
                            found_dummy_snp = True
                            break
            with open(legend, "r") as csvfile:
                reader = csv.reader(csvfile, delimiter=" ")
                for i in range(dummy_snp):
                    next(reader)
                command += " -dl " + next(reader)[1] + " 1 1 1"
                
            # The created genotype will be saved in the file temp/<corpus_id>_chr<chrom>.controls.gen.
            command += " -o ../../../temp/" + str(self.corpus_id) + "_" + self.pop + "_chr" + str(chrom)
            
            # Call HAPGEN2 in len(self.chromosomes) many child processes that run in parallel.
            with open("temp/" + str(self.corpus_id) + "_" + self.pop + "_chr" + str(chrom) + "_hapgen2.log", mode="wb") as logfile:
                p = subprocess.Popen(command, shell=True, stdout=logfile, stderr=logfile)
                child_processes.append(p)
                
        # Wait for HAPGEN2 to terminate.
        for p in child_processes:
            p.wait()
            
        print("Re-compressing HAPMAP3 data ...")
        for chrom in self.chromosomes:
            gen_map = self.hapmap_dir + "genetic_map_chr" + str(chrom) + "_combined_b36.txt"
            if os.path.exists(gen_map + ".bz2"):
                gen_map = ""
            hap = self.hapmap_dir + self.pop + ".chr" +str(chrom) + ".hap"
            if os.path.exists(hap + ".bz2"):
                hap = ""
            legend = self.hapmap_dir + "hapmap3.r2.b36.chr" + str(chrom) + ".legend"
            if os.path.exists(legend + ".bz2"):
                legend = ""
            if (gen_map != "") or (hap != "") or (legend != ""):
                p = subprocess.Popen("bzip2 " + gen_map + " " + hap + " " + legend, shell=True)
                child_processes.append(p)
            
        # Wait for all data to be compressed.
        for p in child_processes:
            p.wait()
            
        
            
    def merge_hapgen2_output(self):
        """Merges the output of HAPGEN2."""
        
        # Clear the genotype array and the list of SNP names.
        self.genotype = []
        self.snps = []
        
        # Merge the genotypes of all chromosomes.
        for chrom in self.chromosomes:
            
            # Read the output generated by HAPGEN2.
            with open("temp/" + str(self.corpus_id) + "_" + self.pop + "_chr" + str(chrom) + ".controls.gen", "r") as csvfile:
                print("Merging HAPGEN2 output for chromosome " + str(chrom) + " ...")
                reader = csv.reader(csvfile, delimiter=" ")
                for row in reader:
                    
                    # The SNP information is contained in the second to fifth column of the output 
                    # generated by HAPGEN2.
                    snp = row[1:5]
                    snp.insert(1, "chr" + str(chrom))
                    self.snps.append(snp)
                    
                    # Retrieve the genotype for the current SNP.
                    # The HAPGEN2 output contains genotype information from position 5 onwards. 
                    # For each individual, the HAPGEN2 output contains a binary triplet with exactly one 1.
                    # (1,0,0) stands for no minor allele, (0,1,0) for one minor allele, and (0,0,1) for 
                    # two minor alleles.
                    row = row[5:]
                    if (len(row) != 3 * self.num_inds):
                        msg = "Input file " + str(self.corpus_id) + "_chr" + str(chrom) + "controls.gen has wrong number of columns.\n"
                        msg += "Expected: " + str(3 * self.num_inds + 5) + ".\n"
                        msg += "Actual: " + str(len(row) + 5) + "."
                        raise ValueError(msg)
                    snp_geno = []
                    for i in range(self.num_inds):
                        for geno_at_snp in range(3):
                            if int(row[3 * i + geno_at_snp]) == 1:
                                snp_geno.append(geno_at_snp)
                                break
                            
                    # Append the genotype of the current SNP to the genotype corpus.
                    self.genotype.append(snp_geno)
        
        # Set the number of SNPs.
        self.num_snps = len(self.snps)
        
        # Transform to genotype to numpy.array.
        self.genotype = np.asarray(self.genotype, dtype=np.uint8)
        
    
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
        
            
        