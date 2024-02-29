# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 17:52:37 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk


Dictionary containing the crrespondence between the CSV files and the respective
shp file of municipalities. This is because in the IBGE's website (as of 08/07/2020)
the municipalities' shapefiles are not available for all years and only from 
2000 onwards. As well, the shapefiles are not updated every year.

https://www.ibge.gov.br/en/geosciences/territorial-organization/territorial-meshes/18890-municipal-mesh.html?=&t=downloads

(UPDATED on 17-08-2023 from the funciton in ./_chapter01/scripts/municipalities_population-density.py)

"""


def municipality_csv_to_shp():

    correspondence = {
            1985 : 1980,
            1986 : 1980,
            1987 : 1980,
            1988 : 1980,
            1989 : 1980,
            1990 : 1980,
            1991 : 1991,
            1992 : 1991,
            1993 : 1991,
            1994 : 1991,
            1995 : 1991,
            1996 : 1991,
            1997 : 2000,
            1998 : 2000,
            1999 : 2000,
            2000 : 2000,
            2001 : 2001,
            2002 : 2001,
            2003 : 2001,
            2004 : 2001,
            2005 : 2005,
            2006 : 2005,
            2007 : 2007,
            2008 : 2007,
            2009 : 2007,
            2010 : 2010,
            2011 : 2010,
            2012 : 2010,
            2013 : 2013,
            2014 : 2013,
            2015 : 2015,
            2016 : 2016,
            2017 : 2017,
            2018 : 2018,
            2019 : 2019,
            2020 : 2020
            }
    
    return correspondence



