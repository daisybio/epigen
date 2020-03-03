#//////////////////////////////////////////////////////////////////////////#
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
#//////////////////////////////////////////////////////////////////////////#

"""Run this script to test to carry out runtime tests.

**Usage**:: 
    
    python3 test_runtime.py
    
**Output:**
    File: 
        ``./num_inds_vs_runtime.csv``
    Content: 
        Runtimes for simulating 1 dataset with varying number of individuals with number of SNPs fixed to 100000.
    File: 
        ``./num_snps_vs_runtime.csv``
    Content: 
        Runtimes for simulating 1 dataset with varying number of SNPs with number of individuals fixed to 10000.
    File: 
        ``./num_datasets_vs_runtime.csv``
    Content: 
        Runtimes for simulating varying number of datasets with number of SNPs and individuals fixed to 100.
"""

import subprocess
import time

def run_script():
    """Runs the script."""
    
    num_inds_vs_runtime = {10 ** i : 0 for i in range(2,5)}
    num_snps_vs_runtime = {10 ** i : 0 for i in range(2,6)}
    num_datasets_vs_runtime = {10 ** i : 0 for i in range(3)}
    command = "python3 simulate_data.py --corpus-id 1 --pop ASW --model models/param_model.xml --inds {} --snps {} --num-sims {}"
    for num_inds in num_inds_vs_runtime:
        start = time.perf_counter()
        subprocess.call(command.format(num_inds, 100000, 1), shell=True)
        end = time.perf_counter()
        num_inds_vs_runtime[num_inds] = end - start
        subprocess.call("rm sim/0_1_ASW.json", shell=True)
    for num_snps in num_snps_vs_runtime:
        start = time.perf_counter()
        subprocess.call(command.format(10000, num_snps, 1), shell=True)
        end = time.perf_counter()
        num_snps_vs_runtime[num_snps] = end - start
        subprocess.call("rm sim/0_1_ASW.json", shell=True)
    for num_datasets in num_datasets_vs_runtime:
        start = time.perf_counter()
        subprocess.call(command.format(100, 100, num_datasets), shell=True)
        end = time.perf_counter()
        num_datasets_vs_runtime[num_datasets] = end - start
        for sim_id in range(num_datasets):
            subprocess.call("rm sim/{}_1_ASW.json".format(sim_id), shell=True)
    with open("num_inds_vs_runtime.csv", mode="w") as output:
        output.write("num_inds,runtime\n")
        for num_inds in num_inds_vs_runtime:
            output.write("{},{}\n".format(num_inds, num_inds_vs_runtime[num_inds]))
    with open("num_snps_vs_runtime.csv", mode="w") as output:
        output.write("num_snps,runtime\n")
        for num_snps in num_snps_vs_runtime:
            output.write("{},{}\n".format(num_snps, num_snps_vs_runtime[num_snps]))
    with open("num_datasets_vs_runtime.csv", mode="w") as output:
        output.write("num_datasets,runtime\n")
        for num_datasets in num_datasets_vs_runtime:
            output.write("{},{}\n".format(num_datasets, num_datasets_vs_runtime[num_datasets]))

if __name__ == "__main__":
    run_script()