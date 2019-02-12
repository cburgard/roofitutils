#!/bin/evn python

inf = float("inf")
nan = float("nan")

def findcontours(points,values,smooth):
    """find the contours in a 2d graph"""
    from numpy import array
    from math import isnan
    keys = sorted(points.keys())

    xvals = [ p[0] for p in keys ]
    yvals = [ p[1] for p in keys ]
    zvals = [ points[p] for p in keys ]

    npoints = 100
    
    from scipy.interpolate import griddata
    from numpy import mgrid
    grid_x, grid_y = mgrid[min(xvals):max(xvals):complex(npoints), min(yvals):max(yvals):complex(npoints)]

    grid_z = griddata(array(keys),array(zvals),(grid_x, grid_y), method='linear')

    nllmin = inf
    for i in range(0,npoints):
        for j in range(0,npoints):
            zval = grid_z[i][j]
            if zval < nllmin:
                nllmin = zval
                xval = grid_x[i][j]
                yval = grid_y[i][j]
    minimum = (xval,yval)
                
    from skimage import measure
    allcontours = []
    allimgcontours = []
    for v in values:
        imgcontours = measure.find_contours(grid_z, v+nllmin)
        allimgcontours.append(imgcontours)
        contours = []
        for c in imgcontours:
            realc = []
            for i,j in c:
                if isnan(i) or isnan(j): continue
                realc.append((i/npoints * (max(xvals)-min(xvals)) + min(xvals),j/npoints * (max(yvals)-min(yvals)) + min(yvals)))
            if len(realc) > 0:
                if smooth:
                    contours.append(smoothgraph(realc))
                else:
                    contours.append(realc)                    
                    
        allcontours.append(contours)

    # Display the image and plot all contours found
#    import matplotlib.pyplot as plt
#    fig, ax = plt.subplots()
#    ax.imshow(grid_z, interpolation='nearest', cmap=plt.cm.gray)
#    for contours in allimgcontours:
#        for n, contour in enumerate(contours):
#            ax.plot(contour[:, 1], contour[:, 0], linewidth=2)
#        
#    ax.axis('image')
#    ax.set_xticks([])
#    ax.set_yticks([])
#    plt.show() 

    return allcontours,minimum
    
def findcrossings(points,nllval):
    """find the minimum point of a 1d graph and the crossing points with a horizontal line at a given value"""
    from scipy.interpolate import PchipInterpolator as interpolate
    xvals = [ x for x,y in points ]
    yvals = [ y for x,y in points ]
    i0 = 0
    ymin = inf
    for i in range(0,len(xvals)):
        if yvals[i] < ymin:
            ymin = yvals[i]
            i0 = i
    x0 = xvals[i0]
    xl,xr = min(xvals),max(xvals)
    r = abs(xr - xl)
    xll = xl - 0.1*r
    xrr = xr + 0.1*r
    interp = interpolate(xvals, yvals, extrapolate=True)
    from scipy.optimize import ridder as solve
    up = nan
    down = nan
    try:
        if   interp(xl ) > nllval:      down = solve(lambda x : interp(x) - nllval,xl ,x0)
        elif interp(xll) > nllval:      down = solve(lambda x : interp(x) - nllval,xll,x0)
        if   interp(xr ) > nllval:      up   = solve(lambda x : interp(x) - nllval,x0,xr )
        elif interp(xrr) > nllval:      up   = solve(lambda x : interp(x) - nllval,x0,xrr)
        r = (x0,x0-down,up-x0)
        return r
    except ValueError as err:
        print("unable to determine crossings")
        return (x0,nan,nan)


def findminimum(points):
    """find the interpolated minimum of a 1d graph"""
    from scipy.interpolate import PchipInterpolator as interpolate
    from scipy.optimize import minimize
    from numpy import array
    xvals = sorted(points.keys())
    yvals = [ points[x] for x in xvals ]
    interp = interpolate(xvals, yvals, extrapolate=True)
    minimum = minimize(lambda v:interp(v[0]), array([min(xvals)]))
    return minimum.fun

def smoothgraph(graph):
    """smooth a graph by taking the center of each edge and constructing a new graph from those"""
    newgraph = []
    from math import isnan
    lastx,lasty = graph[-1]
    for x,y in graph:
        if not isnan(lastx) and not isnan(lasty):
            newx = 0.5*(x+lastx)
            newy = 0.5*(y+lasty)
            newgraph.append((newx,newy))
            lastx,lasty = x,y
    newgraph.append(newgraph[0])
    return newgraph       
