from shapely.geometry import LineString, Point
import numpy as np

def line2pts(geom, dl=None):
    """Given an input Linestring geom, generate points at fixed interval dl in line units
    if no dl is specified, 1000 equidistant samples are created.
    code verified for geom of type: shapely Linestring
    
    Useful for extracting profile data from raster
    """
    #Extract list of (x,y) tuples at nodes
    nodes = [(x,y) for (x,y) in geom.coords]
    #print "%i nodes" % len(nodes)

    #Point spacing in map units
    if dl is None:
        nsteps=1000
        dl = geom.length/nsteps

    #Initialize empty lists
    l = [] #length along line
    mX = [] #x point along line
    mY = [] #y point along line

    #Add first point to output lists
    l += [0]
    x = nodes[0][0]
    y = nodes[0][1]
    mX += [x]
    mY += [y]

    #Remainder
    rem_l = 0
    #Previous length (initially 0)
    last_l = l[-1]
    
    #Loop through each line segment in the feature
    for i in range(0,len(nodes)-1):
        x1, y1 = nodes[i]
        x2, y2 = nodes[i+1]
      
        #Total length of segment
        tl = np.sqrt((x2-x1)**2 + (y2-y1)**2)

        #Number of dl steps we can fit in this segment
        #This returns floor 
        steps = int((tl+rem_l)/dl)

        if steps > 0:
            dx = ((x2-x1)/tl)*dl
            dy = ((y2-y1)/tl)*dl
            rem_x = rem_l*(dx/dl)
            rem_y = rem_l*(dy/dl)
            
            #Loop through each step and append to lists
            for n in range(1, steps+1):
                l += [last_l + (dl*n)]
                #Remove the existing remainder
                x = x1 + (dx*n) - rem_x
                y = y1 + (dy*n) - rem_y
                mX += [x]
                mY += [y]

            #Note: could just build up arrays of pX, pY for entire line, then do single z extraction
            #Update the remainder
            rem_l += tl - (steps * dl)
            last_l = l[-1]
        else:
            rem_l += tl 

    return l, mX, mY

def make_perp_transects(xpts, ypts, L, DEM_res):
    '''
    This function creates perpindicular transects of length L
    at certain locations along a channel centerline defined by xpts y pts
    (Pair with line2pts function above to get equispaced points along centerline)
    
    Parameters
    ------------
    xpts: list of points defining center of transect x locations
    ypts: list of points defining center of transect y locations
    L: length of transect
    DEM_res: resolution to extract topography at each transect(set to DEM res for finest interpolation)

    Output
    --------
    transects: a list of tuples defining the points along each transect of form (x, y, d) where
    x is the x location in native coordinate system units
    y is the y location in native coordinate system units
    d is the distance along the transect
    
    '''
    
    
    transects = [] #list of tuples containing x,y arrays of points within each transect
    slopes = []
    
    for i in range(1, len(xpts)-1): #n transects is two less than n points since using midpoint
            
        #calculate slope of transect for surrounding points
        S = (ypts[i+1] - ypts[i-1]) / (xpts[i+1] - xpts[i-1])
        St = -1/S
    
        #create distances to shift x and y pts away from line vertex
        x_shift = L / (2 * np.sqrt(1 + (St)**2)) 
        y_shift = St * x_shift
        
        #create endpoints of transect by moving points away from second point on centerline vector
        x1, y1 = (xpts[i] + (x_shift)), (ypts[i] + (y_shift)) 
        x2, y2 = (xpts[i] - (x_shift)), (ypts[i] - (y_shift))
        
        
        #create shapely line from this (uses shapely LineString and Point objects)
        # if statement keeps 0 point of each transect on left bank
        if ypts[i+1] > ypts[i-1]: #if river is flowing in positive y direction
            line = LineString([Point(x2, y2), Point(x1, y1)])
        elif ypts[i+1] < ypts[i-1]: #if river is flowing in negative y direction
            line = LineString([Point(x1, y1), Point(x2, y2)])
            
        #put intermediate points on line at resolution of DEM using line2pts function
        trans_l, trans_x, trans_y = line2pts(line, DEM_res)
        transects.append((trans_x, trans_y, trans_l))
    
    return transects
