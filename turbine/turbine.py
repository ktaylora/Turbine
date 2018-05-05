# !/usr/bin/env python3

"""A shell interface to the Turbine 'R' package

Uses argparse to handle runtime arguments passed by the user and then wraps calls to the 'Rscript' interface
commonly deployed on unix systems (Linux / macOS)
"""

import argparse
import subprocess

class RLauncher:
    def __init__(self, *args):
        """Creates an 'R' launcher task with optional parameters specifying the\
        path of our local Rscript binary and the script and arguments we intend to\
        run

        Optional arguments:
        :param r_script_path: full path to the local Rscript binary
        :param src_path: full path to the R source script we intend to run
        """
        self._r_script_path = None
        self._src_path = None

    @property
    def script_path(self):
        """Get the full path to the R script binary path we are going to use for running our R code"""
        return self._script_path

    @script_path.setter
    def script_path(self, *args):
        """Set the full r_script_path we are going to use
        Optional arguments:
        :arg: full path to the local Rscript binary
        """
        self._script_path = args[0]

    @property
    def workspace_path(self):
        """Get the full path to the R script we are going to run for this session"""
        return self._src_path

    @workspace_path.setter
    def workspace_path(self, *args):
        """Set the full path to the R script we are going to run for this session"""

    def run(self):
        """ call Rscript on a user-specified script path and poll the run """
        pass

if __name__ == "__main__":
    pass
