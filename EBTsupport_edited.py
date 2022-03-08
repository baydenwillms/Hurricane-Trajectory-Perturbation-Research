# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 08:54:05 2021

@author: Bayden Willms, Mohamed Iskandarani
"""

import numpy as np
daysPerMonth    =np.array([31,28,31,30,31,30,31,31,30,31,30,31])
daysPerMonthLeap=np.array([31,29,31,30,31,30,31,31,30,31,30,31])

# list of parameters to identify the values held in the columns of the drifters matrix
icode = 0   # code column index
iname = 1   # hurricane name
idate = 2   # Code for mmddhh
iyear = 3   # year
ilat = 4   # latitude column index
ilon = 5   # longitude column index
iwmx = 6   # max wind
islp = 7   # sea level pressure

def timeinseconds(timestring):
    """ Returns the time provided in timestring in the form hh:mm:ss in second
    usage is 
    timeinsec = timeinseconds(timestring)
    """
    timelist = timestring.split(':')
    timeinsec = float(timelist[0])*3600.0 + float(timelist[1])* 60.0 + float(timelist[2])
    return timeinsec

def dayofyear(datestring):
    """ Returns the day provided in datestring in the form yyyy/mm/dd in days
    from the beginning of the year. Usage is: 
    day = dayofyear(datestring)
    """
    datelist=datestring.split('/')
    year = int(datelist[0])
    mm = int(datelist[1])
    dd = int(datelist[2])
    if (year%4==0):
       day=daysPerMonthLeap[0:mm-1].sum()+dd
    else:
       day=daysPerMonth[0:mm-1].sum()+dd
#   print('dayofyear=',dayofyear)
    return float(day)

def ReadEBTFile(fname):
    """ Read a specific carthe file
        returns arrays containing the days, times, latd and lond of drifter
        observations usage is
        ddata=ReadCartheFile(fname) where
        fname is a string containing the path to the file
    """ 
    fileObj = open(fname,'r');             # open file
    lines = fileObj.readlines();           # read the lines
    fileObj.close()                        # close the file

    nbofHeaderLines = 0                 # number of header lines
    nbofLines = len(lines)              # number of items in the 
    ndata = nbofLines - nbofHeaderLines  # number of data in the file
    print("file has ", nbofLines, " lines and ", ndata," data lines")

    hurricaneList = []    
    linedata = lines[nbofHeaderLines].split() # split line into a word list
    hcoderef = linedata[icode]
    hnameref = linedata[iname]
    hyearref = linedata[iyear]       
    wmxlist = []
    lonlist = []
    latlist = []
    slplist = []
    datelist = []
    latmax = -99999
    lonmax = -99999
    lonmin =  99999
    latmin =  99999
    for j in range(ndata):
#    for j in range(66):
       linedata = lines[j+nbofHeaderLines].split() # split line into a word list
       hcode = linedata[icode]
       hname = linedata[iname]
       hyear = linedata[iyear]       

       hmdt  = linedata[idate]
       hmm = hmdt[0:2]
       hdd = hmdt[2:4]
       hhh = str(hmdt[4:6]) 
  #original:     #datestring = hyear+'-'+hmm+'-'+hdd
       datestring = hmm+'/'+hdd+'/'+hyear + '/' + hhh   #Changed this to include hours 3-10
       hlat  = float(linedata[ilat])
       hlon  =-float(linedata[ilon])
       if (hlon >  180.0):
           hlon-= 360.0
       if (hlon < -180.0):
           hlon+= 360.0
       hwmx  = float(linedata[iwmx])
       hslp  = float(linedata[islp])
       #day = dayofyear(datestring)
       #date = day + hhh/24.0
       date = datestring
       
#edited to allow for the date to be viewed in the plots
           
#       print('j:',j,hname,hyear,date,hlon,hlat,hwmx,hslp)
       
       if (hcode!=hcoderef):
           hurricane = {"code":hcoderef, "name":hnameref, "year":hyearref, "lat":np.array(latlist),"lon":np.array(lonlist),"wmx":np.array(wmxlist),"slp":np.array(slplist),"time":np.array(datelist)}
           print(hurricane["name"],hurricane["year"],len(hurricane["lat"]),hurricane["wmx"].max(),hurricane["slp"].min())
           hurricaneList.append(hurricane)
           latmax = max(latmax,hurricane["lat"].max())
           lonmax = max(lonmax,hurricane["lon"].max())
           latmin = min(latmin,hurricane["lat"].min())
           lonmin = min(lonmin,hurricane["lon"].min())
           
           wmxlist.clear()
           lonlist.clear()
           latlist.clear()
           slplist.clear()
           datelist.clear()
           hcoderef = hcode
           hnameref = hname
           hyearref = hyear

           
       lonlist.append(hlon)
       latlist.append(hlat)
       wmxlist.append(hwmx)
       slplist.append(hslp)
       datelist.append(date)
       
    hurricane = {"code":hcoderef, "name":hnameref, "year":hyearref, "lat":np.array(latlist),"lon":np.array(lonlist),"wmx":np.array(wmxlist),"slp":np.array(slplist),"time":np.array(datelist)}
    print(hurricane["name"],hurricane["year"],len(hurricane["lat"]),hurricane["wmx"].max(),hurricane["slp"].min())
    hurricaneList.append(hurricane)
    latmax = max(latmax,hurricane["lat"].max())
    lonmax = max(lonmax,hurricane["lon"].max())
    latmin = min(latmin,hurricane["lat"].min())
    lonmin = min(lonmin,hurricane["lon"].min())    
    
    print("Read data for ",len(hurricaneList)," hurricanes",lonmin,lonmax,latmin,latmax)

    return hurricaneList

def CountTrajectoryLength(drifters,drifterid,ibeg):
    """ Find the last index of a trajectory of a given drifter whose start trajectory is at row ibeg
    Usage is iend =CountTrajectoryLength(drifters,drifterid,ibeg) where
    ibeg is the row where the trajectory history starts
    drifterid is the drifter ID
    drifters is the matrix of drifter data
    iend is the last row of the trajectory + 1
    """
    iend=ibeg                            # initialize trajectory end row to current row
    while (drifters[iend,0]==drifterid): # repeat until drifter-tag changes
       iend = iend + 1                   # increment row
       if (iend == drifters.shape[0]):   # stop before the end of the rows
          break
    return iend

def GetDrifterTrajectorySize(drifterlistbeg,drifterlistend,drifters):
    """ Count the number of drifters and set the lists giving the indices of
the start and end of the trajectories of each drifter. Usage is
    ndrifters,longesttrajectory = GetDrifterTrajectorySize(drifterlistbeg,drifterlistend,drifters)
    """

    ndrifters = 0            # initialize drifter count
    iend = 0                 # initialize the end of the first trajectory
    longesttrajectory=0      # initialize the max-length of trajectories
    while (iend < drifters.shape[0]):
       ndrifters = ndrifters + 1       # increment drifter count
       ibeg = iend                     # update the start of the NEXT trajectory
       drifterid = drifters[ibeg,0]    # get the new drifterID
       iend = CountTrajectoryLength(drifters,drifterid,ibeg) # find the last row of the trajectory
       drifterlistbeg.append(ibeg)     # append the start information to trajectory start-list
       drifterlistend.append(iend)     # append the end information to trajectory end-list
       print(drifterid,ibeg, iend, iend-ibeg) # print some info
       if (iend-ibeg > longesttrajectory): # update the longuest trajectory
          longesttrajectory = iend-ibeg

    return ndrifters, longesttrajectory # returns the number of drifters and longest trajectory
