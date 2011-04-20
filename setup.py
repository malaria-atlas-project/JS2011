# Author: Pete gething
# Date: 20th April 2011
# License: Creative Commons BY-NC-SA
####################################

from setuptools import setup
from numpy.distutils.misc_util import Configuration
import os
config = Configuration('JS2011',parent_package=None,top_path=None)

#config.add_extension(name='variograms.directions',sources=['map_utils/variograms/directions.f'])

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(  version="0.1",
            description="The Malaria Atlas Project's utility functions.",
            author="Peter Gething and Anand Patil", 
            author_email="map@map.ox.ac.uk",
            url="www.map.ox.ac.uk",
            packages=['JS2011','JS2011/EXTRACT','JS2011/SIMULATE'],
            license="Public domain",
            **(config.todict()))


