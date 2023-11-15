# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: py:light,ipynb
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.4
# ---


"""
@author: VR Feb 2021
"""

import sys
import os
from datetime import datetime

sys.path.append(os.path.dirname(__file__))

from . import modules
from . import scripts

from . import version

now = datetime.now()

# define version
version, release_date = version.versionstrg()

__version__ = version
__release__ = release_date
