"""
Downloaded from https://gist.github.com/jsanz/64d08d8ff2d07f619d45e17087d07573
on the 21/09/2021

Small script to fix geometries of the first file argument
using the native QGIS processing algorithm. You may need
to adjust the path to your installation.
"""

import sys
# sys.path.append('/usr/share/qgis/python/plugins')
# sys.path.append('C:/Users/scat8298/AppData/Roaming/QGIS/QGIS3/profiles/default/python/plugins')
# sys.path.append('C:/OSGeo4W64/apps/Python37')
# sys.path.append('C:/OSGEO4~1/apps/qgis-ltr/python/plugins')
# sys.path.append('C:/Users/scat8298/.conda/envs/conda-qgis/Library/python/qgis')
#import qgis.core

# GOOD TUTORIAL: https://docs.qgis.org/testing/en/docs/pyqgis_developer_cookbook/raster.html
# https://www.geodose.com/p/pyqgis.html

# This needs to go first (https://gis.stackexchange.com/questions/348140/qgis-3-10-python-ide-application-path-not-initialized)


# THIS IS A WORKING VERSION!!!!!!!
# TRY TO SEE IF I COULD CALL IT FROM ANOTHER PYTHON SCRIPT
# should I use qgis environment? or can I use the geo_env...?
# what about starting a new environment with geo_env and qgis? and the changes that I have made

# https://gist.github.com/jsanz/64d08d8ff2d07f619d45e17087d07573
# Setting new conda environment with qgis: https://stackoverflow.com/questions/35622661/import-qgis-modules-into-python-anaconda/67447061#67447061
# together with 
# https://veillecarto2-0.fr/2019/12/11/tutorial-create-an-anaconda-environment-for-pyqgis-development/

"""
Imports
"""

from qgis.core import (
    QgsApplication,
    QgsProcessingFeedback,
    QgsVectorLayer
)



print("Initializing QGIS...")
qgs = QgsApplication([], False)
QgsApplication.setPrefixPath('C:/Users/scat8298/.conda/envs/conda-qgis/Library', True)

qgs.initQgis()
from qgis.analysis import QgsNativeAlgorithms
QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())


sys.path.append('C:/Users/scat8298/.conda/envs/conda-qgis/Library/python/plugins/processing')


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

print("Running the fix geometries algorithm...")
res = processing.run("native:fixgeometries", params)#, feedback=feedback)

print("Done!")
