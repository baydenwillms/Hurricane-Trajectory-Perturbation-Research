# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 01:28:36 2021

@author: Bayden Willms
"""

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def MapInitPlot(lonlim,latlim):
    """ Plot the Gulf of Mexico a trajectorie and save plot to file
        usage is MapPlotTrajectory(lons,lats)
        where lons,lats are 1D numpay arrays holding the longitudes and latitudes of the trajectories in degrees
              u,v are 1D numpy arrays holding the velocity components 
              stride is an integer to control the drawing of vectors
                     0 does not draw any vector
                     1 draws every vector
                     2 draws every other vector
       the image is saved to file 'driftertrajectory.png'
    """
    fig = plt.figure()                     # create a figure
#    projObj = ccrs.Mercator()              # define a type of projection and store its parameters
    projObj = ccrs.PlateCarree()           # define PlateCarree projection
    loncenter = 0.5* (lonlim[0]+lonlim[1])
    latcenter = 0.5* (latlim[0]+latlim[1])
#    projObj = ccrs.Orthographic(loncenter,latcenter)
    ax = plt.axes(projection=projObj)      # generate an axis given the projObj
#   ax.set_global()
    ax.set_extent(lonlim+latlim, crs=projObj)           # define extent of the map
#    ax.stock_img()
    ax.coastlines(resolution='50m')        # draw the coastlines overriding the default 110m resolution
    ax.add_feature(cfeature.LAND)          # add land feature using the default land color
    ax.add_feature(cfeature.OCEAN)         # add ocean feature using the default ocean color
    ax.gridlines(draw_labels=True)         # add grid lines with labels (either gridlines or x,y ticks
#   ax.set_xticks(np.arange(lonlim[0],lonlim[1],5.0))
#   ax.set_yticks(np.arange(latlim[0],latlim[1],5.0))
    fig.canvas.draw()
    plt.pause(0.1)
    return fig,ax, projObj

#def MapPlotTrajectory(lons,lats, projObj, ax,fig):
def MapPlotTrajectory(lons, lats, projObj, ax, fig, colorstr='r', linewidthint=2):
    """ Plot the  trajectorie 
        usage is MapPlotTrajectory(lons,lats, projObj, ax,fig)
        where lons,lats are 1D numpay arrays holding the longitudes and latitudes of the trajectories in degrees
    """
  # modifying to plot both trajectories
    line, = ax.plot(lons,lats,color=colorstr, linewidth=linewidthint, transform=projObj)        # plots trajectory as a red line
#    fig.canvas.draw()
    plt.pause(0.1)
    return line

def MapPlotNewTrajectory(new_lons, new_lats, projObj, ax, fig):
    """ Plot the  trajectorie 
        usage is MapPlotTrajectory(lons,lats, projObj, ax,fig)
        where lons,lats are 1D numpay arrays holding the longitudes and latitudes of the trajectories in degrees
    """
  # modifying to plot both trajectories
    line, = ax.plot(new_lons,new_lats,'b+',transform=projObj)        # plots trajectory as a blue X
#    fig.canvas.draw()
    plt.pause(0.1)
    return line

#edited for 4-6 to include trajectory perturbations
###################################################################
def MapPlotTrajectoryPTB(new_PTB_coords, projObj, ax, fig):
    """ Plot the  trajectorie 
        usage is MapPlotTrajectory(lons,lats, projObj, ax,fig)
        where lons,lats are 1D numpay arrays holding the longitudes and latitudes of the trajectories in degrees
    """
  # modifying to plot both trajectories
    line = ax.plot(new_PTB_coords[0,:], new_PTB_coords[1,:],'b',transform=projObj)        # plots trajectory as a blue X
#    fig.canvas.draw()
    plt.pause(0.1)
    return line

###################################################################
def MapPlotTrajectoryAnimate(lons,lats, projObj, ax,fig):
    """ Plot the Gulf of Mexico trajectorie of a single drifter
        usage is MapPlotTrajectory(lons,lats,projObj, ax,fig)
        where lons,lats are 1D numpay arrays holding the longitudes and latitudes in degrees
              projObj is a projection object returned by cartopy.projection()
              ax is an axis object
              fig is a figure object where the map will be drawn
    """
    tail, = ax.plot(lons,lats,'b',transform=projObj)
    lastPosition, = ax.plot(lons[0],lats[0],'r+',transform=projObj)        # plots trajectory as a red line
    fig.canvas.draw()
    plt.pause(0.1)
    ntail = 200
    for j in range(len(lons)):
      jbeg = max(0,j-ntail)
      # tail.set_xdata(lons[jbeg:j])
      # tail.set_ydata(lats[jbeg:j])
      lastPosition.set_xdata(lons[j])
      lastPosition.set_ydata(lats[j])
      fig.canvas.draw()
      plt.pause(0.1)




#MapPlotTrajectoryAnimate(lons,lats, projObj, ax,fig)
