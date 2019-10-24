# !/usr/bin/env python2

"""A shell interface to the Turbine 'R' package

Uses argparse to handle runtime arguments passed by the user and then wraps calls to the 'Rscript' interface
commonly deployed on unix systems (Linux / macOS)
"""

import argparse
import subprocess

from gee_asset_manager.batch_remover import delete
from gee_asset_manager.batch_uploader import upload
from gee_asset_manager.config import setup_logging


class ee_ingest:
    def __init__(self, *args):
        try:
            self._status = None
            self._asset_id = args[0]
            self._filename = args[1]
        except Exception as e:
            raise e
        # call-out built-in method for our upload task
        self.run_upload_task()

    def run_upload_task(self, *args):
        upload(
            user=args.user,
            source_path=self._filename,
            destination_path=self._asset_id,
            metadata_path=None,
            multipart_upload=None,
            nodata_value=None,
            bucket_name=None,
            band_names=None,
        )


class RLauncher:
    def __init__(self, *args):
        """Creates an 'R' launcher task with optional parameters specifying \
        the path of our local Rscript binary and the script and arguments we \
        intend to run

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
