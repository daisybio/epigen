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

"""This script constitutes the main user interface of EpiGEN -- run it, to simulate epistasis data.

This script has to be run on top of a pre-computed genotype corpus. For each population code ``<POP>``
and each chromosome ``<CHROM>``, EpiGEN contains a pre-computed corpus for 10000 individuals. These corpora
can be selected by running the script with the options ``--corpus-id <CRHOM> --pop <POP>``.
You can also use your own corpora -- simply run the script ``generate_genotype_corpus.py`` before running this script.

**Usage**::

    python3 simulate_data.py [required arguments] [optional arguments]

**Required Arguments:**
    ``--corpus-id CORPUS_ID``
        Description:
            ID of selected genotype corpus.
        Accepted Arguments:
            Non-negative integers.
        Effect:
            Together with ``--pop``, this option selects the genotype corpus with the prefix ``./corpora/<CORPUS_ID>_<POP>``. 
            If this corpus does not exist, the script raises an error. If necessary, run the script ``./generate_genotype_corpus.py``
            to generate the desired corpus.
    ``--pop POP``    
        Description:
            HAPMAP3 population code of selected genotype corpus.
        Accepted Arguments: 
            ASW, CEU, CEU+TSI, CHD, GIH, JPT+CHB, LWK, MEX, MKK, TSI, and MIX (for merged corpora).
        Effect:
            Together with ``--corpus-id``, this option selects the genotype corpus with the prefix ``./corpora/<CORPUS_ID>_<POP>``. 
            If this corpus does not exist, the script raises an error. If necessary, run the script ``./generate_genotype_corpus.py``
            to generate the desired corpus.
    ``--model MODEL``
        Description:
            Path to epistasis model given as INI or XML file.
            INI files are used to provide extensionally defined models, XML files specify parametrized models.
        Accepted Arguments:
            Strings ending in ``.ini`` or ``.xml`` that represent paths to existing model files.
        Effect:
            Specifies the epistasis model used by the simulator.
        Format of INI files for extensionally defined models:
            For each genotype of length ``<size>``, parameters of a Normal distribution must be provided.
            INI files for quantitative phenotypes have to be of the following format:::
            
                [Model Type]
                size = <size>
                phenotype = quantitative 
                [Model Definition]
                0,...,0 = <mu>,<stdev>
                .
                .
                .
                2,...,2 = <mu>,<stddev>
                 
            For each genotype of length ``<size>``, parameters of a categorical distribution must be provided
            INI files for categorical phenotypes have to be of the following format:::
            
                [Model Type]
                size = <size>
                phenotype = <c> 
                [Model Definition]
                0,...,0 = <p_1>,...,<p_c>
                .
                .
                .
                2,...,2 = <p_1>,...,<p_c>
                
        Format of XML files for parametrized models:
            XML files have to match the document type definition ``models/ParametrizedModel.dtd``.
        Examples:
            Cf. the files ``models/ext_model.ini`` and ``models/param_model.xml``.        
    ``--snps SNPS``:
        Description:
            Number of SNPs in simulated data.
        Accepted Arguments:
            Positive integers. If larger than the number of SNPs in the selected corpus, it is lowered to this number.
            Should be set to a number that is significantly smaller than the number of SNPs in the selected corpus,
            because otherwise, EpiGEN's subsampling techniques have no effect. 
        Effect:
            Determines how many SNPs from the selected corpus are included in the simulated data.
    ``--inds INDS``
        Description:
            Number of individuals in simulated data.
        Accepted Arguments:
            Positive integers. If larger than the number of individuals in the selected corpus, it is lowered to this number.
            Should be set to a number that is significantly smaller than the number of individuals in the selected corpus,
            because otherwise, EpiGEN's subsampling techniques have no effect. 
        Effect:
            Determines how many individuals from the selected corpus are included in the simulated data.
            
**Required Group of Mutually Exclusive Arguments:**
    ``--sim-ids SIM_ID [SIM_ID ...]``
        Description:
            IDs of simulated data.
        Accepted Arguments:
            One or several non-negative integers.
        Effect: 
            If N IDs <SIM_ID_1> ... <SIM_ID_N> are provided, N simulated epistasis instances with prefices 
            ``./sim/<SIM_ID_1>_<CORPUS_ID>_<POP>``, ..., ``./sim/<SIM_ID_N>_<CORPUS_ID>_<POP>`` are generated.
    ``--num-sims NUM_SIMS``
        Description:
            The number of epistasis instances that should be generated.
        Accepted Arguments:
            Positive integer.
        Effect:
            If specified, <NUM_SIMS> simulated epistasis instances with prefices
            ``./sim/0_<CORPUS_ID>_<POP>``, ..., ``./sim/<NUM_SIMS-1>_<CORPUS_ID>_<POP>`` are generated.

**Optional Arguments:**
    ``--global-maf-range LB UB``
        Description:
            Range of acceptable MAFs for noise SNPs.
        Accepted Arguments:
            Floats ``<LB>`` and ``<UB>`` with ``0 <= <LB> < <UB> <= 1``.
        Default:
            [0,1]
        Effect:
            All SNPs except for the disease SNPs are randomly sampled from those SNPs in the corpus whose MAFs 
            fall into the specified range. If the range is too narrow, it is dynamically extended at runtime.
    ``--biased-distr PARAM [PARAM ...]``
        Description:
            Biased target distribution for simulated phenotypes.
        Accepted Arguments for Quantitative Phenotypes:
            A white-space separated list of floats of length 2 whose elements represent the mean (first element) and 
            standard deviation (second element) of a Normal distribution.
        Accepted Arguments for Categorical Phenotypes:
            A white-space separated list of floats of length ``<c>`` whose entries represent the probabilities of the ``<c>`` 
            categories.
        Effect:
            If provided, the individuals are subsampled after generating the phenotypes such that the obtained phenotype distribution 
            matches the biased distribution. This option can hence be used to model observation bias.
    ``--seed SEED``
        Description:
            Seed for numpy.random. 
        Accepted Arguments:
            Non-negative integers.
        Effect:
            If provided, the simulator always generates the same data given the same input.
    ``--compress``
        Description:
            Compress the generated output files.
        Accepted Arguments:
            None.
        Effect:
            Determines the suffix ``<SUFFIX>`` of the generated files. If provided, ``<SUFFIX>`` is set to ``json.bz2``.
            Otherwise, it is set to ``json``.
    ``-h, --help``
        Effect:
            Show help message and exit.

**Optional Mutually Exclusive Arguments:**
    ``--disease-snps SNP [SNP ...]``
        Description:
            Position of disease SNPs in selected genotype corpus. 
        Accepted Arguments: 
            White space separated list of non-negative integers whose length matches the size of the model specified 
            via the option ``--model``. All integers must be smaller than the number of SNPs in the selected corpus.
        Effect:
            If provided, the selected SNPs form the disease SNP set employed by the simulator.
    ``--disease-maf-range LB UB``
        Description:
            Range of acceptable MAFs for disease SNPs.
        Accepted Arguments:
            Floats ``<LB>`` and ``<UB>`` with ``0 <= <LB> < <UB> <= 1``.
        Default:
            [0,1]
        Effect:
            Unless ``--disease-snps`` is provided, the disease SNPs are randomly sampled from those SNPs in the corpus whose MAFs 
            fall into the specified range. If the range is too narrow, it is dynamically extended at runtime.

**Output:**
    *Output File:*
        Path:
            ``./sim/<ID>_<CORPUS_ID>_<POP>.<SUFFIX>``
        Format for Quantitative Phenotypes:
            (Compressed) JSON file of the form ``{"num_snps": <NUM_SNPS>, "num_inds": <NUM_INDS>, "model_type": "quantitative", "genotype": <GENOTYPE_DATA>, "phenotype": <PHENOTYPE_DATA>, "snps": <SNP_DATA>, "disease_snps": <DISEASE_SNP_DATA>, "mafs": <MAF_DATA>}``.
        Format for Categorical Phenotypes:
            (Compressed) JSON file of the form ``{"num_snps": <NUM_SNPS>, "num_inds": <NUM_INDS>, "model_type": "categorical", "num_categories": <NUM_CATEGORIES>, "genotype": <GENOTYPE_DATA>, "phenotype": <PHENOTYPE_DATA>, "snps": <SNP_DATA>, "disease_snps": <DISEASE_SNP_DATA>, "mafs": <MAF_DATA>}``. 
    ``<NUM_SNPS>``
        Key:
            ``"num_spns"``
        Content and Format:
            Integer representing the number of SNPs.
    ``<NUM_INDS>``
        Key:
            ``"num_inds"``
        Content and Format:
            Integer representing the number of individuals.
    ``<NUM_CATEGORIES>``
        Key:
            ``"num_categories"``
        Content and Format:
            Integer representing the number of categories for categorical phenotypes.
    ``<GENOTYPE_DATA>``
        Key: 
            ``"genotype"``
        Content and Format: 
            JSON field of the form ``[[G_0_0, ..., G_0_<INDS-1>] ... [G_<SNPS-1>_0, ..., G_<SNPS-1>_<INDS-1>]]``,
            where ``G_S_I`` encodes the number of minor alleles of the individual with index ``I`` at the SNP with index ``S``.
    ``<PHENOTYPE_DATA>``
        Key: 
            ``"phenotype"``
        Content and Format: 
            JSON field of the form ``[P_0, ..., P_<INDS-1>]``, where ``P_I`` encodes the phenotype of 
            the individual with index ``I``.
    ``<SNP_DATA>``
        Key:
            ``"snps"``
        Content and Format:
            JSON field of the form ``[INFO_0, ..., INFO_<SNPS-1>]``, where ``INFO_S`` contains the following
            information about the SNP with index ``S``: RS identifier, chromosome number, position on chromosome, major allele, minor allele.
    ``<DISEASE_SNP_DATA>``
        Key: 
            ``"disease_snps"``
        Content and Format:
            JSON field of the form ``[S_0, ..., S_<MODELSIZE-1>]``, where ``S_POS`` encodes the index of the SNP at position ``POS``
            in the epistasis model.
    ``<MAF_DATA>``
        Key:
            ``"mafs"``
        Content and Format:
            JSON field of the form ``[MAF_0, ..., MAF_<SNPS-1>]``, where ``MAF_S`` encodes the MAF of
            the SNP with index ``S``.
    
"""

from utils.data_simulator import DataSimulator as DataSim
import utils.argparse_checks as checks
import argparse

def run_script():
    """Runs the script."""
    
    descr = "\n############################################################################\n"
    descr += "######################### EpiGEN - simulate_data.py ########################\n"
    descr += "\nRun this script to simulate epistasis data.\n"
    descr += "\nusage: python3 %(prog)s [required arguments] [optional arguments]\n"
    epilo = "The simulated data can be found in the ./sim directory:\n"
    epilo += "Generated data:\t./sim/<ID>_<CORPUS_ID>_<POP>.<SUFFIX>\n"
    epilo += "\n############################################################################\n"
    parser = argparse.ArgumentParser(description=descr,formatter_class=argparse.RawTextHelpFormatter, epilog=epilo, usage=argparse.SUPPRESS)
    required_args = parser.add_argument_group("required arguments")
    required_args.add_argument("--corpus-id", type=int, required=True, help="ID of selected genotype corpus.", action=checks.check_non_negative("--corpus-id"))
    required_args.add_argument("--pop", required=True, choices=["ASW","CEU","CEU+TSI","CHD","GIH","JPT+CHB","LWK","MEX","MKK","TSI","MIX"], metavar="POP", help="HAPMAP3 population code of selected genotype corpus.")
    required_args.add_argument("--model", required=True, help="Path to model file (see examples in ./model directory).")
    required_args.add_argument("--snps", type=int, required=True, help="Number of SNPs in simulated data.", action=checks.check_positive("--snps"))
    required_args.add_argument("--inds", type=int, required=True, help="Number of individuals in simulated data.", action=checks.check_positive("--inds"))
    required_exclusive_args = parser.add_argument_group("required group of mutually exclusive arguments")
    sim_options = required_exclusive_args.add_mutually_exclusive_group(required=True)
    sim_options.add_argument("--sim-ids", type=int, nargs="+", help="IDs of simulated data.", action=checks.check_non_negative("--sim-ids"))
    sim_options.add_argument("--num-sims", type=int, help="Number of simulations.", action=checks.check_positive("--num-sims"))
    optional_args = parser.add_argument_group("optional arguments")
    optional_args.add_argument("--noise-maf-range", type=float, nargs=2, default=[0,1], metavar=("LB", "UB"), help="Range of acceptable MAFs for noise SNPs. Default = [0,1].", action=checks.check_interval("--global-maf-range"))
    optional_args.add_argument("--biased-distr", default=[], type=float, nargs="+", action=checks.check_length("--biased-distr"), help="Biased target distribution for simulated phenotypes. Default = [].", metavar="PARAM")
    optional_args.add_argument("--seed", type=int, default=None, help="Seed for numpy.random. Default = None.", action=checks.check_non_negative("--seed"))
    optional_args.add_argument("--compress", help="Compress generated output files.", action="store_true")
    exclusive_args = parser.add_argument_group("optional group of mutually exclusive arguments")
    disease_snps = exclusive_args.add_mutually_exclusive_group()
    disease_snps.add_argument("--disease-snps", default=[], type=int, nargs="+", action=checks.check_length("--disease-snps"), metavar="SNP", help="Position of disease SNPs in selected genotype corpus. Default = [].")
    disease_snps.add_argument("--disease-maf-range", type=float, nargs=2, default=[0,1], metavar=("LB", "UB"), help="Range of acceptable MAFs for disease SNPs. Default = [0,1].", action=checks.check_interval("--disease-maf-range"))
    args = parser.parse_args()
    
    print("\n############################################################################")
    print("######################### EpiGEN - simulate_data.py ########################\n")
    sim_ids = []
    if args.sim_ids:
        sim_ids = args.sim_ids
    else:
        sim_ids = [i for i in range(args.num_sims)]
    sim = DataSim(args.corpus_id, args.pop, args.model, args.snps, args.inds, args.disease_snps, args.biased_distr, args.noise_maf_range, args.disease_maf_range, args.seed, args.compress)
    for index in range(len(sim_ids)):
        sim.set_sim_id(sim_ids[index])
        sim.sample_snps()
        sim.generate_phenotype()
        sim.dump_simulated_data()
    suffix = "json"
    if args.compress:
        suffix = "json.bz2"
    print("----------------------------------------------------------------------------")
    print("Finished simulation of epistasis data.")
    print("The generated data can be found in the ./sim directory:")
    for sim_id in sim_ids:
        print("Generated data:\t./sim/{}_{}_{}.{}".format(sim_id, args.corpus_id, args.pop, suffix))
    print("\n############################################################################")
    
if __name__ == "__main__":
    run_script()