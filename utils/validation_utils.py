# //////////////////////////////////////////////////////////////////////////#
#                                                                          #
#   Copyright (C) 2020 by David B. Blumenthal                              #
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
# //////////////////////////////////////////////////////////////////////////#

"""Contains utility functions for validating simulated data."""

import os.path
import json
import bz2
import numpy as np
import itertools as itt
import scipy.stats as sps


def load_data(filename):
    """Loads simulated data.

    Args:
        filename (str): Path to the file containing the data.

    Returns:
        dict: A dictionary containing the simulated data.
    """
    if not os.path.exists(filename):
        raise ValueError(f'The input file {filename} does not exist.')
    if filename.endswith('.json.bz2') :
        with bz2.open(filename, "rt", encoding="ascii") as zipfile:
            sim_data = json.load(zipfile)
    elif filename.endswith('.json'):
        with open(filename, "rt") as jsonfile:
            sim_data = json.load(jsonfile)
    else:
        ValueError(f'Invalid input file {filename}. Expecting .json or .json.bz2 file.')
    return sim_data


def generate_penetrance_table(sim_data):
    """Generates penetrance table.

    Args:
        sim_data (dict): The simulated data.

    Returns:
        dict: The penetrance table induced by the disease SNPs.
    """
    genotypes = np.asarray(sim_data['genotype'], dtype=np.uint8)
    phenotypes = sim_data['phenotype']
    disease_snps = sim_data['disease_snps']
    num_inds = sim_data['num_inds']
    penetrance_table = {genotype: [] for genotype in itt.product(range(3), repeat=len(disease_snps))}
    for ind in range(num_inds):
        genotype = tuple(genotypes[disease_snps, ind].tolist())
        penetrance_table[genotype].append(phenotypes[ind])
    return penetrance_table


def one_way_anova(penetrance_table):
    """Carries out one-way ANOVA F-test.

    Args:
         penetrance_table (dict): The penetrance table induced by the disease SNPs.

    Returns:
        float: The obtained p-value.
    """
    non_empty_cells = [phens for phens in penetrance_table.values() if len(phens) > 0]
    _, p_value = sps.f_oneway(*non_empty_cells)
    return p_value


def chi_square(penetrance_table, num_categories):
    """Carries out chi-square test.

        Args:
             penetrance_table (dict): The penetrance table induced by the disease SNPs.

        Returns:
            float: The obtained p-value.
        """
    non_empty_cells = [phens for phens in penetrance_table.values() if len(phens) > 0]
    contingency_table = np.zeros(shape=(len(non_empty_cells), num_categories), dtype=int)
    cell_id = 0
    for cell in non_empty_cells:
        for phenotype in cell:
            contingency_table[cell_id, phenotype] += 1
        cell_id += 1
    _, p_value, _, _ = sps.chi2_contingency(contingency_table, correction=False)
    return p_value


def write_to_log_file(logfilename, test, p_value, disease_mafs, penetrance_table):
    """Writes the results of the validation to a log-file.

    Args:
        logfilename (str): Name of the log-file.
        test (str): Name of the test.
        p_value (float): Obtained p-value.
        disease_mafs (list): List of MAFs of disease SNPs.
        penetrance_table (dict): Penetrance table induced by disease SNPs.
    """
    penetrance_table = {str(gen): phens for gen, phens in penetrance_table.items()}
    log = {'test': test, 'p_value': p_value, 'disease_mafs': disease_mafs, 'penetrance_table': penetrance_table}
    try:
        with open(logfilename, 'w') as jsonfile:
            json.dump(log, jsonfile)
            print(f'Wrote log-file to {logfilename}.')
    except IOError:
        print(f'Writing log-file to {logfilename} failed.')
