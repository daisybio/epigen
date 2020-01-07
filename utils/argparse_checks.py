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
from builtins import isinstance

"""Contains checks for arguments parsed by argparse."""

import argparse
import collections.abc as abc

def check_length(argname):
    """Ensures that at least two arguments are provided.
    
    Args:
        argname (str): Name of the argparse argument.
    """
    class CheckLength(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if not 2<=len(values):
                msg="The argument \"{f}\" requires at least 2 arguments.".format(f=argname)
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return CheckLength

def check_interval(argname):
    """Ensures that the provided arguments specify a sub-interval of [0,1].
    
    Args:
        argname (str): Name of the argparse argument.
    """
    class CheckInterval(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if not (0 <= values[0] < values[1] <= 1):
                msg="The argument \"{f}\" requires arguments that specify a sub-interval of [0,1].".format(f=argname)
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return CheckInterval

def check_positive(argname):
    """Ensures that the provided argument is positive.
    
    Args:
        argname (str): Name of the argparse argument.
    """
    class CheckPositive(argparse.Action):
        def __call__(self, parser, args, value, option_string=None):
            if isinstance(value, abc.Iterable):
                for item in value:
                    if not (0 < item):
                        msg="The argument \"{f}\" requires positive arguments.".format(f=argname)
                        raise argparse.ArgumentTypeError(msg)
            elif not (0 < value):
                msg="The argument \"{f}\" requires a positive argument.".format(f=argname)
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, value)
    return CheckPositive

def check_non_negative(argname):
    """Ensures that the provided argument is non-negative.
    
    Args:
        argname (str): Name of the argparse argument.
    """
    class CheckNonNegative(argparse.Action):
        def __call__(self, parser, args, value, option_string=None):
            if isinstance(value, abc.Iterable):
                for item in value:
                    if not (0 <= item):
                        msg="The argument \"{f}\" requires non-negative arguments.".format(f=argname)
                        raise argparse.ArgumentTypeError(msg)
            elif not (0 <= value):
                msg="The argument \"{f}\" requires a non-negative argument.".format(f=argname)
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, value)
    return CheckNonNegative