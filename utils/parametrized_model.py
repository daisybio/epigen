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

"""Contains definition of ParametrizedModel class."""

import numpy as np
import xml.etree.ElementTree as ET

class ParametrizedModel(object):
    """Represents a parametrized epistasis model.
    
    Parametrized models can be used to model binary categorical phenotypes.
    This class is employed by the class DataSimulator.
    Not intended for use outside of this class.
    
    Attributes:
        size (int): An integer greater equal 2 representing the size of the model.
        baseline_alpha (float): A float greater 0 representing the baseline risk.
        marginal_models (dict of (int,callable)): The marginal models.
        interaction_models (dict of ((list of int),callable)): The interaction models.
        phenotype (int/str): The inetger 2 (for dichotomous phenotypes) or the 
            string "quantitative" (for quantitative phenotypes).
    """

    def __init__(self, model, seed):
        """Initializes the ExtensionalModel.
        
        Args:
            model (str): Path to an XML file that contains the specification of the parametrized model.
                Examples can be found in the directory epigen/models/.
                
            seed (int/None): The seed for np.random (possibly None).
        """
        
        # Set the seed.
        if seed != None:
            np.random.seed(seed)
            
        # Initialize the interaction model.
        self.baseline_alpha = 1.0
        self.marginal_models = {}
        self.interaction_models = {}
        
        # Parse the XML file.
        try:
            model = ET.parse(model).getroot()
        except:
            raise ValueError("Cannot read XML file " + model + ".")
        
        # Ensure that the XML document is of the correct type.
        if model.tag != "ParametrizedModel":
            raise ValueError("The root of the XML must be called \"ParametrizedModel\".")
        
        # Retrieve the size of the model.
        try:
            self.size = int(model.get("size"))
        except:
            raise ValueError("The XML element \"ParametrizedModel\" must contain an attribute \"size\" convertible to an int greater equal 2.")
        if self.size < 2:
            raise ValueError("The XML element \"ParametrizedModel\" must contain an attribute \"size\" convertible to an int greater equal 2.")
        
        # Retrieve the phenotype.
        try:
            self.phenotype = model.get("phenotype")
        except:
            raise ValueError("The XML element \"ParametrizedModel\" must contain an attribute \"phenotype\" which either equals \"2\" or \"quantitative\".")
        if self.phenotype != "2" and self.phenotype != "quantitative":
            raise ValueError("The XML element \"ParametrizedModel\" must contain an attribute \"phenotype\" which either equals \"2\" or \"quantitative\".")
        if self.phenotype == "2":
            self.phenotype = 2
            
        # Retrieve the standard deviation for generating quantitative phenotypes.
        self.stdev = None
        try:
            self.stdev = model.get("stdev")
        except:
            if self.phenotype == "quantitative":
                raise ValueError("For quantitative phenotypes, the XML element \"ParametrizedModel\" must contain an attribute \"stdev\" convertible to a positive float.")
        if self.stdev != None:
            try:
                self.stdev = float(self.stdev)
            except:
                raise ValueError("For quantitative phenotypes, the XML element \"ParametrizedModel\" must contain an attribute \"stdev\" convertible to a positive float.")
            if self.stdev <= 0:
                raise ValueError("For quantitative phenotypes, the XML element \"ParametrizedModel\" must contain an attribute \"stdev\" convertible to a positive float.")
        
        # Retrieve the baseline model, the marginal models, and the interaction models.
        baseline_set = False
        for child in model:
            # Retrieve the alpha of the model.
            try:
                alpha = float(child.get("alpha"))
            except:
                raise ValueError("The XML element \"" + child.tag + "\" must contain an attribute \"alpha\" convertible to a float greater equal 0.")
            if alpha < 0:
                raise ValueError("The XML element \"" + child.tag + "\" must contain an attribute \"alpha\" convertible to a float greater equal 0.")
            
            # Retrieve alpha of baseline model.
            if child.tag == "BaselineModel":
                if baseline_set:
                    raise ValueError("The XML file contains several XML elements \"ParametrizedModel\".")
                baseline_set = True
                self.baseline_alpha = alpha
            else:
                # Retrieve the type of the marginal or interaction model.
                try:
                    model_type = child.get("type")
                except:
                    raise ValueError("The XML element \"" + child.tag + "\" must contain an attribute \"type\".")
                # Collect the SNP indices to which the model should be applied.
                poss = []
                for pos in child:
                    if pos.tag != "pos":
                        raise ValueError("Unexpected XML element \"ParametrizedModel." + child.tag + "." + pos.tag + "\".")
                    try:
                        new_pos = int(pos.text)
                    except:
                        raise ValueError("The value of the XML element \"ParametrizedModel." + child.tag + ".pos\" must be convertible to an int between 0 and the size of the model - 1.")
                    if new_pos < 0 or new_pos >= self.size:
                        raise ValueError("The value of the XML element \"ParametrizedModel." + child.tag + ".pos\" must be convertible to an int between 0 and the size of the model - 1.")
                    poss.append(new_pos)
                sorted_poss = tuple(sorted(set(poss)))
                if len(sorted_poss) < len(poss):
                    raise ValueError("The XML element \"ParametrizedModel." + child.tag + "\" contains duplicate SNPs.")
                # Add the parsed marginal model to the dict of marginal models.
                if child.tag == "MarginalModel":
                    if len(sorted_poss) > 1:
                        raise ValueError("The XML element \"ParametrizedModel." + child.tag + "\" contains multiple SNPs.")
                    if sorted_poss[0] in self.marginal_models:
                        raise ValueError("Multiple definition of marginal model for SNP " + str(sorted_poss[0]) + ".")
                    self.marginal_models[sorted_poss[0]] = self._make_marginal_model(alpha, model_type, sorted_poss[0])
                # Add the parsed interaction model to the dict of interaction models.
                if child.tag == "InteractionModel":
                    if len(sorted_poss) < 2:
                        raise ValueError("The XML element \"ParametrizedModel." + child.tag + "\" contains less than two SNPs.")
                    if sorted_poss in self.interaction_models:
                        raise ValueError("Multiple definition of interaction model for SNP set " + str(sorted_poss) + ".")
                    self.interaction_models[sorted_poss] = self._make_interaction_model(alpha, model_type, list(sorted_poss))
        
    def __call__(self, gen_at_snp_set):
        """Generates the phenotype for a genotype at a SNP set of appropriate size.
        
        Args:
            gen_at_snp_set (numpy.array): A one-dimensional numpy.array of length self.size containing integers between 0 and 2,
                which represents the genotype of an individual at a set of disease SNPs.
            
        Returns:
            float/int: The generated phenotype, i.e. a float if self.phenotype equals "quantitative", 
                and an integer otherwise from {0,1}, otherwise.
        """
        
        # Compute the joint parameter.
        risk = self.baseline_alpha
        for item in self.marginal_models.items():
            risk *= item[1](gen_at_snp_set)
        for item in self.interaction_models.items():
            risk *= item[1](gen_at_snp_set)
            
        # Generate a quantitative phenotype.
        if self.phenotype == "quantitative":
            return np.random.normal(loc=risk, scale=self.stdev)    
        
        # Generate a categorical phenotype.
        prob = risk / (1.0 + risk)
        return np.random.choice(a=2, p = [1 - prob, prob])
            
    def _make_marginal_model(self, alpha, model_type, pos):
        """Generates a function representing a marginal model.
        
        Args:
            alpha (float): A non-negative float representing the alpha of the marginal model.
            model_type (str): A string representing the model type. Valid choices: "additive", "dominant", or "recessive".
            pos (int): An integer in range(self.size) representing the index of the SNP to which the marginal model should be applied.
            
        Returns:
            callable: A function that, given a one-dimensional numpy.array of length self.size, returns the alpha corresponding to the selected model.
        """
        
        def marginal_model(gen_at_snp_set):
            if model_type == "dominant":
                if gen_at_snp_set[pos] == 0:
                    return 1
                else:
                    return alpha
            elif model_type == "recessive":
                if gen_at_snp_set[pos] == 2:
                    return alpha
                else:
                    return 1
            elif model_type == "additive":
                if gen_at_snp_set[pos] == 0:
                    return 1
                else:
                    return (alpha * float(gen_at_snp_set[pos])) / 2.0
            else:
                msg = "Unsupported model type " + model_type + " for marginal model. "
                msg += "Choices: \"dominant\", \"recessive\", \"additive\"."
                raise ValueError(msg)
        return marginal_model
    
    def _make_interaction_model(self, alpha, model_type, poss):
        """Generates a function representing an interaction model.
        
        Args:
            alpha (float): The risk of the interaction model.
            model_type (str): The model type. Valid choices: "exponential", "multiplicative", "joint-dominant", or "joint-recessive".
            poss (list of int): A list of pairwise different integers from range(self.size) with length between 2 and self.size whose entries represent
                the indices of the SNPs to which the interaction model should be applied.
            
        Returns:
            callable: A function that, given a one-dimensional numpy.array of length self.size, returns the alpha corresponding to the selected model.
        """
        
        def interaction_model(gen_at_snp_set):
            if model_type == "multiplicative":
                return alpha ** np.sum(gen_at_snp_set[poss])
            elif model_type == "exponential":
                return alpha ** np.prod(gen_at_snp_set[poss])
            elif model_type == "joint-dominant":
                if np.prod(gen_at_snp_set[poss]) == 0:
                    return 1
                else:
                    return alpha
            elif model_type == "joint-recessive":
                if np.prod(gen_at_snp_set[poss] == 2) == 0:
                    return 1
                else:
                    return alpha
            else:
                msg = "Unsupported model type " + model_type + " for interaction model. "
                msg += "Choices: \"joint-dominant\", \"joint-recessive\", \"multiplicative\", \"exponential\"."
                raise ValueError(msg)
        return interaction_model