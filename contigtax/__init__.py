# -*- coding: utf-8 -*-
import pkg_resources  # part of setuptools
__version__ = pkg_resources.require("contigtax")[0].version


def diamond_legacy():
    """
    Diamond versions older than 0.9.19 required the taxonmap to be specified
    at the aligner stage.

    :return: string to be passed to diamond aligner
    """
    # Check diamond version
    from distutils.version import StrictVersion
    import subprocess
    r = subprocess.run(["diamond", "version"], stdout=subprocess.PIPE)
    v = r.stdout.decode().rstrip().rsplit()[-1]
    # Older versions required you to supply a taxonmap at the aligner stage
    if StrictVersion(v) < StrictVersion("0.9.19"):
        return True
    return False


"""Main package for contigtax"""
name = "contigtax"
