# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 14:45:35 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Function to generate biplots with arows (PCA-style)

"""

"""
Imports
"""
import seaborn as sns


"""
Functions
"""
def biplot(scores, components, col_names, scores_colours, components_labels, palette, axes, scale = 1, inverse = False):
    
    xs = scores['{}1'.format(col_names)]
    ys = scores['{}2'.format(col_names)]
    cs = scores[scores_colours]
    
    coefx = components['{}1'.format(col_names)]
    coefy = components['{}2'.format(col_names)]
    coefl = components[components_labels]
    
    if inverse == True:
        xs = -1 * xs
        ys = -1 * ys
        
        coefx = -1 * coefx
        coefy = -1 * coefy
    
    # Creating the scatter plot of PCA scores
    sns.scatterplot(x = xs, 
                    y = ys,
                    hue = cs, 
                    palette = palette, ax = axes, s = 5, 
                    edgecolor = 'none', legend = False)
    
    # Adding a grid
    # axes.grid()
    axes.axvline(x = 0, c = 'grey', linestyle = '--', linewidth = 0.5)
    axes.axhline(y = 0, c = 'grey', linestyle = '--', linewidth = 0.5)
    
    # Adding an arrow per variable in pca_comp
    for i in range(len(components)):
        axes.arrow(0, 0, coefx[i]/scale, coefy[i]/scale,
                  color = 'r', alpha = 0.5)
        axes.text(coefx[i]*1.05/scale, coefy[i]*1.05/scale, 
                 coefl[i], color = 'g', ha = 'center', va = 'center')
