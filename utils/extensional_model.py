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

"""Contains definition of ExtensionalModel class."""

import numpy as np
import configparser as cp

class ExtensionalModel(object):
    """Represents an extensionally defined epistasis model.
    
    Extensional models can be used to model binary or non-binary categorical phenotypes, as well as quantitative phenotypes.
    This class is employed by the class DataSimulator.
    Not intended for use outside of this class.
    
    Attributes:
        size (int): An integer greater equal 2 representing the size of the model.
        model (numpy.array): A numpy.array of dimension self.size representing the epistasis model.
        phenotype (int/str): An integer greater equal 2 (for categorical phenotypes) 
            or the string "quantitative" (for quantitative phenotypes).
    """

    def __init__(self, model, seed):
        """Initializes the ExtensionalModel.
        
        Args:
            model (str): Path to an INI file that contains the full extensional definiton of an epistasis model.
                For categorical models, a discrete probability distribution must be provided for each possible genotype of dimension self.size. 
                For quantitative models, the mean and the standard deviation of normal distributions must be provided.
                Examples can be found in the directory epigen/models/.
                
            seed (int/None): The seed for np.random (possibly None).
        """
        
        # Set the seed.
        if seed != None:
            np.random.seed(seed)
        self.epsilon = 0.00001
        
        # Read the INI file.
        config = cp.ConfigParser()
        try:
            config.read(model)
        except:
            raise ValueError("Cannot read INI file " + model + ".")
        
        # Ensure that it contains the section "Model Type".
        if not "Model Type" in config:
            print(config.sections())
            raise ValueError("The INI file must contain the section [Model Type].")
        if not "size" in config["Model Type"]:
            raise ValueError("The INI file must contain the key \"size\" in the section [Model Type].")
        
        # Retrieve the size of the model.
        try:
            self.size = config.getint("Model Type", "size")
        except:
            raise ValueError("The value of the key \"size\" must be convertible to an int greater equal 2.")
        if self.size < 2:
            raise ValueError("The value of the key \"size\" must be convertible to an int greater equal 2.")
        
        # Retrieve the phenotype type.
        self.phenotype = config.get("Model Type", "phenotype")
        if self.phenotype != "quantitative":
            try:
                self.phenotype = int(self.phenotype)
            except:
                raise ValueError("The value of the key \"phenotype\" must be either \"quantitative\" or convertible to an int greater equal 2.")
            if self.phenotype < 2:
                raise ValueError("The value of the key \"phenotype\" must be either \"quantitative\" or convertible to an int greater equal 2.")
        
        # Get the number of model parameters.
        num_def= 0
        num_params = 0
        if self.phenotype == "quantitative":
            num_params = 2
        else:
            num_params = self.phenotype
        
        # Initialize the model as a numpy.array of appropriate shape.
        shape = [3] * self.size
        types = [("param" + str(i), float) for i in range(num_params)]
        self.model = np.empty(shape, dtype=types)
        defined = np.zeros(shape, dtype=bool)
        
        # Retrieve the model definition.
        for key in config["Model Definition"]:
            # Set error messages.
            wrong_key_msg = "The INI file contains a key \"" + key + "\" of the wrong format. "
            wrong_key_msg += "Expected: a comma-separated list of integers between 0 and 2 of length " + str(self.size) + "."
            wrong_value_msg = "The value of the key \"" + key + "\" must be convertible to a comma-separated list of floats of length " + str(num_params)
            if self.phenotype == "quantitative":
                wrong_value_msg += ". The second float must be positive."
            else:
                wrong_value_msg += ". All floats must be between 0 and 1 and sum up to 1."
            # Parse the key.
            try:
                index = tuple([int(i) for i in key.split(",")])
            except:
                raise ValueError(wrong_key_msg)
            if len(index) != self.size:
                raise ValueError(wrong_key_msg)
            for i in index:
                if i < 0 or i > 2:
                    raise ValueError(wrong_key_msg)
            if defined[index]:
                raise ValueError("The INI file contains multiple occurrences of the key \"" + key + "\".")
            # Parse the probability distribution for the current key.
            try:
                params = [float(p) for p in config.get("Model Definition", key).split(",")]
            except:
                raise ValueError(wrong_value_msg)
            if len(self.model[index]) != num_params:
                raise ValueError(wrong_value_msg)
            if self.phenotype == "quantitative":
                if params[1] <= 0:
                    raise ValueError(wrong_value_msg)
            else:
                sum_probs = 0
                for p in params:
                    if p < 0 or p > 1:
                        raise ValueError(wrong_value_msg)
                    sum_probs += p
                if (sum_probs > 1 + self.epsilon) or (sum_probs < 1 - self.epsilon):
                    raise ValueError(wrong_value_msg)      
            self.model[index] = tuple(params)
            defined[index] = True
            num_def += 1
            
        # Ensure that the definition is complete.
        if num_def < 3 ** self.size:
            raise ValueError("The INI file does not contain definitions for all genotypes s in {0,1,2}^" + str(self.size) + ".")
                
    def __call__(self, gen_at_snp_set):
        """Generates the phenotype for a genotype at a SNP set of appropriate size.
        
        Args:
            gen_at_snp_set (numpy.array): A one-dimensional numpy.array of length self.size containing integers between 0 and 2,
                which represents the genotype of an individual at a set of disease SNPs.
            
        Returns:
            float/int: The generated phenotype, i.e. a float if self.phenotype equals "quantitative",
                and an integer between 0 and self.phenotype - 1 otherwise.
        """
        
        # Get the parameters from self.model.
        params = self.model[tuple(gen_at_snp_set.tolist())]
        
        # Generate a quantitative phenotype.
        if self.phenotype == "quantitative":
            return np.random.normal(loc=params[0], scale=params[1])
        
        # Generate a categorical phenotype.
        return np.random.choice(a=self.phenotype, p=list(params))