"""
Downloaded from https://gist.github.com/jsanz/64d08d8ff2d07f619d45e17087d07573
on the 21/09/2021

Small script to fix geometries of the first file argument
using the native QGIS processing algorithm. You may need
to adjust the path to your installation.

Note from Carlota Segura-Garcia: a specific environment to run this script must be created.
I am using the environment "conda-qgis". I created this new conda environment following
https://stackoverflow.com/questions/35622661/import-qgis-modules-into-python-anaconda/67447061#67447061
together with
https://veillecarto2-0.fr/2019/12/11/tutorial-create-an-anaconda-environment-for-pyqgis-development/


Tutorials for PyQGIS and standalone Python scripts calling QGIS
https://docs.qgis.org/testing/en/docs/pyqgis_developer_cookbook/raster.html
https://www.geodose.com/p/pyqgis.html
https://docs.qgis.org/2.18/en/docs/pyqgis_developer_cookbook/intro.html#using-pyqgis-in-standalone-scripts


The key part of this script was extracted from:
    https://gist.github.com/jsanz/64d08d8ff2d07f619d45e17087d07573

"""

import sys
# sys.path.append('/usr/share/qgis/python/plugins')
# sys.path.append('C:/Users/scat8298/AppData/Roaming/QGIS/QGIS3/profiles/default/python/plugins')
# sys.path.append('C:/OSGeo4W64/apps/Python37')
# sys.path.append('C:/OSGEO4~1/apps/qgis-ltr/python/plugins')
# sys.path.append('C:/Users/scat8298/.conda/envs/conda-qgis/Library/python/qgis')
#import qgis.core


"""
Imports
"""

# For standalone Python script, you have to set QGIS Environment path. Code below
# should be added first.
# Information extracted from 
# https://gis.stackexchange.com/questions/348140/qgis-3-10-python-ide-application-path-not-initialized

from qgis.core import (
    QgsApplication,
    QgsProcessingFeedback,
    QgsVectorLayer
)

# Create a reference to the QgsApplication, setting the second argument to False
# disables the GUI
qgs = QgsApplication([], False)
# Supply path to qgis installation location
# The second argument controls whether the default paths are used
QgsApplication.setPrefixPath('C:/Users/scat8298/Miniconda3/envs/conda-qgis/Library', True)
# UPDATE 15/02/2023: Maybe now I have to use this one: C:\Users\scat8298\Miniconda3\envs\conda-qgis
# Older one was: C:/Users/scat8298/.conda/envs/conda-qgis/Library

# Load data providers and layer registry
qgs.initQgis()



# Here starts the script downloaded from 
# https://gist.github.com/jsanz/64d08d8ff2d07f619d45e17087d07573
# Where I have adapted the path to my installation

from qgis.analysis import QgsNativeAlgorithms
QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())

sys.path.append('C:/Users/scat8298/Miniconda3/envs/conda-qgis/Library/python/plugins/processing')
# UPDATE 15/02/2023: Maybe now I have to use this one: C:\Users\scat8298\Miniconda3\envs\conda-qgis
# Old one: C:/Users/scat8298/.conda/envs/conda-qgis/Library/python/plugins/processing

import processing
from processing.core.Processing import Processing

Processing.initialize()



# Getting the file paths
in_file = sys.argv[1]
out_file = sys.argv[2]

# Running the algorithm
params = {
    'INPUT': QgsVectorLayer(in_file, 'layer1', 'ogr'),
    'OUTPUT': out_file
}
feedback = QgsProcessingFeedback()

res = processing.run("native:fixgeometries", params)#, feedback=feedback)


# TODO: just trying to see if this manages to close the files
# When the script is complete, call exitQgis() to remove the provider and layer registries from memory
# qgs.exitQgis()