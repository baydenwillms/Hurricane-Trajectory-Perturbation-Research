# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 21:27:59 2021

@author: Bayden Willms
"""


import numpy as np
#modifying EBTsupport to allow for printing of the actual date on each velocity plot
import EBTsupport_edited as ebte
#modifying trajectoryplotcartopy to map new trajectory vs previous trajectory
import trajectory_modified as mplt
import matplotlib.pyplot as plt
import datetime as datetime
import matplotlib.dates as pltdte
import random

#we are going to be making the output and figures more readable
#we will change the lon and lat limits of the trajectories
#we will add a feature to plot multiple hurricanes on one plot

############################################################  
#edited for 3-30 to include perturbations  

def PlotHurricanePTB(hurricane, latlim, lonlim, beg, end, new_PTB_coords):
    # Plot all hurricanes in hurricaneList
    fig,ax,projObj = mplt.MapInitPlot(lonlim,latlim)
        
    lon = hurricane["lon"][beg:end+1]
    lat = hurricane["lat"][beg:end+1]
    #line = mplt.MapPlotTrajectory(hurricane["lon"], hurricane["lat"], projObj, ax,fig) #changed for 3-9 to include time interval

    for i in range(new_PTB_coords.shape[0]):
        linePTB = mplt.MapPlotTrajectoryPTB(new_PTB_coords[i, :, :], projObj, ax, fig)
        
    line = mplt.MapPlotTrajectory(lon, lat, projObj, ax,fig, linewidthint=2, colorstr='r')
    
    fig.savefig(hurricane["name"]+hurricane["year"]+"TrajectoryPerturbations.png")
    plt.show()
    
#######################################################################

def hurricanedate2mdates(datearray):
    """ Returns the date in matplotlib format for plotting. Usage is
    mdates = dates2datenumbers(yyyy,mm,dd)
    where YYYY is the year, MM is the month and DD is the day, they must be ints
    """
# transform dates to matplotlib.date dates
    mdates = np.empty(len(datearray))
    for i,date in enumerate(datearray):    
        datelist=date.split('/')
        year  = int(datelist[2])
        month = int(datelist[0])
        day = int(datelist[1])
        hh = int(datelist[3])
        ddate = datetime.datetime(year, month, day, hour=hh)
        mdates[i] = pltdte.date2num(ddate)    
    
    return mdates
    
def getNewIndex(hurricane):
    print(hurricane["name"], 'occurred from', hurricane["time"][0], 'until', hurricane["time"][-1])
    flag = 1
#    flag = int(input("Enter '1' to skip/use all of the hurricane's availability. Enter '0' to enter a time interval."))  #to use the entire time interval
    if(flag == 1):
        flag = print("using the entire hurricane's time series.")
        time1index,time2index = 0,len(hurricane["time"])-1
    else: 
 
    #for 3-11 modified this function to take hour inputs from the user
        print("Enter the time interval of the storm you would like to evaluate, in MM/DD/YYYY/HR format")
        time1 = input("Beginning of time interval: ")
        time2 = input("End of time interval: ")
          
        N = len(hurricane["time"])
        for i in range(N):
            if(hurricane["time"][i] == time1):
                time1index = i
                break
        for i in range(N):
            if(hurricane["time"][i] == time2):
                time2index = i
                break
           
    return time1index, time2index
    
def printDates(hurricane):  #added 3-9, used to help user select time interval
    print('')
    print(hurricane["name"], 'occurred from', hurricane["time"][0], 'until', hurricane["time"][-1])
    print('')
    
def FindHurricane(hurricaneList, name, year):
    # Returns hurricane or None
    for hurricane in hurricaneList: #running through hurricanelist
        if (hurricane['name'].upper() == name.upper() and hurricane['year'] == year):  
            return(hurricane)
    return(None)

def ChooseHurricane(hurricaneList):
    # Get user input and return chosen hurricane
    print('Please choose a storm for analysis.')
    while True:
        name = input('Name: ')
        year = input('Year: ')
        hurricane = FindHurricane(hurricaneList, name, year)

        if hurricane is not None:  # if the hurricane is found
            return(hurricane)  # then return that hurricane to be plotted
        else:
            print('Could not find that hurricane, please check your spelling and the year, and try again.')

def PlotHurricanes(hurricaneList, latlim, lonlim):
    # Plot all hurricanes in hurricaneList
    fig,ax,projObj = mplt.MapInitPlot(lonlim,latlim)

    lineList = []
    for j, hurricane in enumerate(hurricaneList):
        line = mplt.MapPlotTrajectory(hurricane["lon"],hurricane["lat"], projObj, ax,fig)
        print(hurricane["name"], hurricane["year"], hurricane["lon"].min(), hurricane["lon"].max(), hurricane["lat"].min(), hurricane["lat"].max())
        lineList.append(line)
        if j > 20:
            break

    #    keyinp = input("Hit any key to do next hurricane")
    plt.show()

#######################################################
def PlotHurricane(hurricane, latlim, lonlim, Vel, new_coords, beg, end):
    # Plot all hurricanes in hurricaneList
    fig,ax,projObj = mplt.MapInitPlot(lonlim,latlim)

        
    lon = hurricane["lon"][beg:end+1]
    lat = hurricane["lat"][beg:end+1]
    #line = mplt.MapPlotTrajectory(hurricane["lon"], hurricane["lat"], projObj, ax,fig) #changed for 3-9 to include time interval
    line, = mplt.MapPlotTrajectory(lon, lat, projObj, ax,fig)
    
    ax.quiver(lon[:-1], lat[:-1], Vel[0,:], Vel[1,:], transform=projObj)
    
  #  line_new_coords = mplt.MapPlotNewTrajectory(new_coords[0,:], new_coords[1,:], projObj, ax, fig)
    print(hurricane["name"], hurricane["year"], hurricane["lon"].min(), hurricane["lon"].max(), hurricane["lat"].min(), hurricane["lat"].max())

    #    keyinp = input("Hit any key to do next hurricane")
    plt.show()
    
    
#added 2-16-2021
def new_trajectory(lonlatinit, dt, REarth, Vel, beg, end):

    deg_2_rad = np.pi/180
    rad_2_deg = 180/np.pi
   # Ntimes = len(lons)
    
    #next line added for 3-11 to solve the index boundary error
    Ntimes = end - beg    
    
    #initial position
    new_coords = np.zeros([2, Ntimes])
    new_coords[:,0] = lonlatinit[:] * deg_2_rad    
    
    for i in range(1, Ntimes): #everything in loop is in radians
        
        #equations 7 and 8, new trajectory (lambda, theta)
        new_coords[1,i] = new_coords[1,i-1] + Vel[1,i-1]*dt / REarth
        temp = np.cos((new_coords[1,i-1] + new_coords[1,i])* 0.5)
        new_coords[0,i] = new_coords[0,i-1] + Vel[0,i-1]* dt / (REarth * temp)

    new_coords *= rad_2_deg #converting to degrees
    
    # print('')
    # print('new_coords size = ', new_coords.shape)
    # print('')
    
    return new_coords

#edited for 3-30-2021 to add perturbations
def PTB_trajectories(lonlatinit, dt, REarth, beg, end, VelPTB):
    """ Calculates trajectories associated with perturbed velocities 
    """
    print('Calculating Trajectory from Hurricane Velocity')
    
    
    #next line added for 3-11 to solve the index boundary error
    Ntimes = VelPTB.shape[-1]   # length of time-series for velocity
    Nperts = VelPTB.shape[0]    # number of realizations

    new_PTB_coords = np.zeros([Nperts, 2, Ntimes + 1]) # create storage for pertubed trajectories
    
    for z in range(Nperts):
        new_PTB_coords[z, :, :] = new_trajectory(lonlatinit, dt, REarth, VelPTB[z, :, :] , 0, Ntimes + 1)

    return new_PTB_coords

#changed for 3-16 to include dates in the plots
def velocity_sphere(lons, lats, dt, REarth):
    print('')
    print("Finite Difference Calculation of Velocity")
    deg_2_rad = np.pi/180
    Ntimes = len(lons)
    
    Vel = np.zeros([2, Ntimes - 1])
    Vel[1,:] = REarth * deg_2_rad * (lats[1:] - lats[0:-1])/dt
    cs = np.cos(deg_2_rad * 0.5 * (lats[0:-1] + lats[1:]))
    Vel[0,:] = REarth * deg_2_rad * (lons[1:] - lons[0:-1])/dt * cs

    return Vel

#modified 3-17 to create our array (Ii + alpha(Xi))
def vel_stats(Vel):
    
    VelMean = Vel.mean(axis=1) # compute time-mean
    VelStd = Vel.std(axis=1)   # compute time-stddev
    print('Mean velocities: ', VelMean)
    print('Standard Deviation: ', VelStd)

#   Calculate anomalies and scale them by stddev    
    VelAnom = np.zeros_like(Vel)  #gives an array of the same shape as Vel 
    VelAnom[0,:] = (Vel[0,:] - VelMean[0]) / VelStd[0]
    VelAnom[1,:] = (Vel[1,:] - VelMean[1]) / VelStd[1]
    
    #plotted anomaly
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.plot(VelAnom[0,:], 'r', VelAnom[1,:], 'k')
    # plt.title('Velocity Anomaly graph')
    # plt.xlabel('Entire Duration of the Storm')
    # plt.ylabel('unitless')
    # fig.show()

    
#Singular value decomposition of Vel
    SVD = np.linalg.svd(VelAnom)
    U = SVD[0] #left eigenvector
    S = SVD[1] # singular values
    W = SVD[2] 
#    print('U: ', U.shape)
#    print('S: ', S.shape)
#    print('W: ', W.shape)
    
    #checking our SVD
    Sd = np.diag(S)
    
    US = np.matmul(U,Sd) 
    USb = np.zeros_like(Vel)
    USb[:,0:2] = US[:,0:2]
    USW = np.matmul(USb, W)
    
    max_error = abs(VelAnom-USW).max()
    print("The maximum error in the SVD calculation: ", max_error)
    
    #edit for 3-23:
    ########################################################################
    Nperts = 30
    Ntimes = Vel.shape[1]
    VelPTB = np.zeros([Nperts, 2, Ntimes])
    alpha = 0.15 #size of the kick
    for r in range(Nperts):    
        PTB = np.zeros_like(U)  #creating the perturbations matrix
        
        Xi1 = random.uniform(-1,1)  #our random perturbation
        Xi2 = random.uniform(-1,1)
        PTB[0,0] = alpha * Xi1 + 1
        PTB[1,1] = alpha * Xi2 + 1  
        
        PTBs = np.matmul(PTB, Sd)
        
        US = np.matmul(U,PTBs) 
        USb = np.zeros_like(Vel)
        USb[:,0:2] = US[:,0:2]
        USW = np.matmul(USb, W)
   
        VelPTB[r, 0,:] = USW[0,:] * VelStd[0] + VelMean[0]
        VelPTB[r, 1,:] = USW[1,:] * VelStd[1] + VelMean[1]
        
    return VelPTB
  
def VelPTBPlot(hurricane, Vel, times, VelPTB=None):

# transform our hurricane dates to matplotlib.date dates
    mdates = hurricanedate2mdates(times)
    
    #some code for reproducing figure 5
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if VelPTB is not None:
       Nperts = VelPTB.shape[0]
       for r in range(Nperts): 
           ax.plot(mdates,VelPTB[r,0,:], 'g', mdates,VelPTB[r, 1,:], 'orange') #meridional and zonal velocities

    ax.plot(mdates,Vel[0,:], 'r', mdates,Vel[1,:], 'r', linewidth=2) #zonal + merid velocities BT
    ax.set_title('Velocity Peturbations')
    ax.set_xlabel('dates')
    ax.set_ylabel('Perturbation Forward Speed (m/s)')
    #format the x-axis ticks and their labels for dates
    myFmt = pltdte.DateFormatter("%Y-%m-%d")
    ax.xaxis.set_major_formatter(myFmt)
    fig.autofmt_xdate()
    
    fig.savefig(hurricane["name"]+hurricane["year"]+'VelocityPerturbations.png')
    fig.show()
       
def main():
    # Read the drifter data. The data is stored in a 2D numpy arrays which holds
    # the following data
    # drifterID timeInDays  Lat Lon   U V
    # readcarthefile define the indices for (time, lat, lon) as (itim, ilat, ilon)
    filename = 'ebtrk_atlc_1988_2018.txt'
    
    #the section of + signs below was edited 2-25-21 in an attempt to show dates on the velocity plots
   #+++++++++++++++++++++++++++++++++++++++++++++++
    
  #original:  #hurricaneList = ebt.ReadEBTFile(filename)
    hurricaneList = ebte.ReadEBTFile(filename)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # ######################################
    # Define map borders in degrees
    # modified on 4-14 to make map limits smaller
    latlim = [5, 60]            # define south/north limits of map
    lonlim = [-105, -15]         # define west/east limits of map
    print('')

  # PlotHurricanes(hurricaneList, latlim, lonlim)     #use if you want to plot the trajectory
    
    dt = 6 * 3600    #timestep in seconds
    REarth = 6371.e3  #radius of Earth in meters 

    nextHurricane = 1
    while nextHurricane > 0:    
        hurricane = ChooseHurricane(hurricaneList)

        printDates(hurricane) #added for 3-9, displays dates of a hurricane, beg and end
    
        beg, end = getNewIndex(hurricane)
        
        times = hurricane["time"][beg:end]
    
        lon = hurricane["lon"][beg:end+1]
        lat = hurricane["lat"][beg:end+1]
    
    # Calculate the hurricane velocity via Finite Differences    
        Vel = velocity_sphere(lon, lat, dt, REarth)
        print("Minimum: ", Vel.min(axis=1))
        print("Maximum: ", Vel.max(axis=1))
        print("Mean: ", Vel.mean(axis=1))
        
     # Velocity perturbations (via SVD decomposition)
        VelPTB = vel_stats(Vel)                   # svd decomposition + perturbations
        VelPTBPlot(hurricane, Vel, times, VelPTB) # plot vel perturbations
     # Compute and plot perturbed trajectories
        lonlatinit = np.array([lon[beg],lat[beg]]) # initial position
        new_PTB_coords = PTB_trajectories(lonlatinit, dt, REarth, beg, end, VelPTB)
        PlotHurricanePTB(hurricane, latlim, lonlim, beg, end, new_PTB_coords)
        nextHurricane = int( input("Enter a + number to do another hurricane:") )
main()
