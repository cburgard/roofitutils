#!/bin/evn python

inf = float("inf")
nan = float("nan")

def findmergecontours(points1, points2, values, smooth):
    """find the contours in a 2d graph"""
    from numpy import array
    from math import isnan
    keys1 = sorted(points1.keys())
    keys2 = sorted(points2.keys())

    xvals1 = [ p[0] for p in keys1 ]
    yvals1 = [ p[1] for p in keys1 ]
    zvals1 = [ points1[p] for p in keys1 ]


    xvals2 = [ p[0] for p in keys2 ]
    yvals2 = [ p[1] for p in keys2 ]
    zvals2 = [ points2[p] for p in keys2 ]

    npoints = 50

    from scipy.interpolate import griddata
    from numpy import mgrid
    grid1_x, grid1_y = mgrid[min(min(xvals1),min(xvals2)):max(max(xvals2),max(xvals1)):complex(npoints), min(min(yvals1),min(yvals2)):max(max(yvals2),max(yvals1)):complex(npoints)]
    grid2_x, grid2_y = mgrid[min(min(xvals1),min(xvals2)):max(max(xvals2),max(xvals1)):complex(npoints), min(min(yvals1),min(yvals2)):max(max(yvals2),max(yvals1)):complex(npoints)]
#    grid1_x, grid1_y = mgrid[min(min(xvals1),min(xvals2)):max(max(xvals2),max(xvals1)):complex(npoints), 0.0:max(max(yvals2),max(yvals1)):complex(npoints)]
#    grid2_x, grid2_y = mgrid[min(min(xvals1),min(xvals2)):max(max(xvals2),max(xvals1)):complex(npoints), 0.0:max(max(yvals2),max(yvals1)):complex(npoints)]

    print min(min(xvals1),min(xvals2))
    print max(max(xvals2),max(xvals1))
    print min(min(yvals1),min(yvals2))
    print max(max(yvals2),max(yvals1))

    grid1_z = griddata(array(keys1),array(zvals1),(grid1_x, grid1_y), method='cubic')
    grid2_z = griddata(array(keys2),array(zvals2),(grid2_x, grid2_y), method='cubic')

    #find minima from the grid
    nllmin1_grid = inf
    minimum1_grid = (nan,nan)
    for i in range(0,npoints):
        for j in range(0,npoints):
            zval = grid1_z[i][j]
            if zval < nllmin1_grid:
                nllmin1_grid = zval
                minimum1_grid = (grid1_x[i][j],grid1_y[i][j])

    nllmin2_grid = inf
    minimum2_grid = (nan,nan)
    for i in range(0,npoints):
        for j in range(0,npoints):
            zval = grid2_z[i][j]
            if zval < nllmin2_grid:
                nllmin2_grid = zval
                minimum2_grid = (grid2_x[i][j],grid2_y[i][j])

#perform the evelope on the grid
    count = 0
    for i in range(0,npoints):
        for j in range(0,npoints):
	    if grid1_x[i][j] == grid2_x[i][j] and grid1_y[i][j] == grid2_y[i][j]:
                 if grid2_z[i][j] != nan and grid2_z[i][j] != nan:
                    grid1_z[i][j] = min(grid2_z[i][j],grid1_z[i][j])
#		    grid1_z[i][j] = grid2_z[i][j]
#		    grid1_z[i][j] = grid1_z[i][j]
		    if -0.001 < grid1_x[i][j] < 0.001 and -0.001 <  grid1_y[i][j] < 0.001:
			print grid1_z[i][j]
		    count = count + 1

    print count
    #find nllmin in the enveloped grid
    minimum = (nan,nan)
    nllmin = inf
    for i in range(0,npoints):
        for j in range(0,npoints):
            zval = grid1_z[i][j]
            if zval < nllmin:
                nllmin = zval
                minimum = (grid1_x[i][j],grid1_y[i][j])

 #   print nllmin
#    nllmin = 11112234.9491
#    nllmin = 11112234.9635
#    nllmin = 11112236.5508

#    nllmin = 11112236.5508
#    nllmin = 11112234.9515
    nllmin = 2156657.18102
    from skimage import measure
    allcontours = []
    allimgcontours = []
    for v in values:
        imgcontours = measure.find_contours(grid1_z,v + nllmin)
        allimgcontours.append(imgcontours)
        contours = []
        for c in imgcontours:
            realc = []
            for i,j in c:
                if isnan(i) or isnan(j): continue
                realc.append((i/npoints * (max(max(xvals1),max(xvals2))-min(min(xvals1),min(xvals2))) + min(min(xvals1),min(xvals2)),j/npoints * (max(max(yvals1),max(yvals2))-min(min(yvals1),min(yvals2))) + min(min(yvals1),min(yvals2))))
#                realc.append((i/npoints * (max(xvals1)-min(xvals1)) + min(xvals1),j/npoints * (max(yvals1)-min(yvals1)) + min(yvals1)))
            if len(realc) > 0:
                if smooth:
                    contours.append(smoothgraph(realc))
                else:
                    contours.append(realc)                    
                    
        allcontours.append(contours)


    # Display the image and plot all contours found
#    import matplotlib.pyplot as plt
#    fig, ax = plt.subplots()
#    ax.imshow(grid2_z, interpolation='nearest', cmap=plt.cm.coolwarm)
#    ax.set_zlim(11112230, 11112245)
#    for contours in allimgcontours:
#        for n, contour in enumerate(contours):
#            ax.plot(contour[:, 1], contour[:, 0], linewidth=2)
       
#    ax.axis('image')
#    ax.set_xticks([])
#    ax.set_yticks#([])
#    plt.show() 

#    from mpl_toolkits.mplot3d import Axes3D
#    import matplotlib.pyplot as plt
#    from matplotlib import cm
#    from matplotlib.ticker import LinearLocator, FormatStrFormatter
#    import numpy as np

#    fig = plt.figure()
#    ax = fig.gca(projection='3d')

    # Plot the surface.
#    surf = ax.plot_surface(grid1_x, grid1_y, grid1_z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
#    surf1 = ax.plot_surface(grid2_x, grid2_y, grid2_z, cmap=cm.Spectral, linewidth=0, antialiased=False)


    # Customize the z axis.
#    ax.set_zlim(11112230, 11112245)
#    ax.zaxis.set_major_locator(LinearLocator(10))
#    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
#    surf1.set_clim(11112230,11112245)
#  
#    surf.set_clim(11112230,11112245)
# Add a color bar which maps values to colors.
#    fig.colorbar(surf, shrink=0.5, aspect=5)
#    plt.show()


    return allcontours,minimum
 

def findcontours(points,values,smooth):
    """find the contours in a 2d graph"""
    from numpy import array
    from math import isnan
    keys = sorted(points.keys())

    xvals = [ p[0] for p in keys ]
    yvals = [ p[1] for p in keys ]
    zvals = [ points[p] for p in keys ]
    xyvals = [ [p[0],p[1]] for p in keys ]
    npoints = 50
    from scipy.interpolate import griddata
    from numpy import mgrid
    grid_x, grid_y = mgrid[min(xvals):max(xvals):complex(npoints), min(yvals):max(yvals):complex(npoints)]
#    grid_x, grid_y = mgrid[-0.20:0.5:complex(npoints), 0.0:0.5:complex(npoints)]
    grid_z = griddata(array(keys),array(zvals),(grid_x, grid_y), method='linear')

#    nllmin = 11112234.4895
#    nllmin = 2156649.49679 #hww
#    nllmin = 11112236.5508 #GenP
#   nllmin = 11112234.9515
#    nllmin = 2156657.16678 
#    nllmin = 2156657.18102
#    nllmin = 2156649.49679 #hgamgam
#    nllmin = 2156649.49679 #htautau
#    nllmin = 2156657.18102 #5xs
#    nllmin = 11112234.526  #CgCy
#    nllmin = 2156657.18102 #5xs old
    nllmin = 2156657.16678 #5xs new 
#    minimum = (nan,nan)
#    for i in range(0,npoints):
#        if zvals[i] < nllmin:
#            nllmin = zvals[i]
#            minimum = (xvals[i],yvals[i])

    nllmin_grid = inf
    minimum_grid = inf
    for i in range(0,npoints):
        for j in range(0,npoints):
            zval = grid_z[i][j]
            if zval < nllmin_grid:
                nllmin_grid = zval
                minimum_grid = (grid_x[i][j],grid_y[i][j])

#    for i in range(0,npoints):
#	for j in range(0,npoints):
#	    grid_z[i][j] = grid_z[i][j] - nllmin_grid + nllmin
#    x_sm = inf
#    y_sm = inf
#    z_sm = inf
#    for i in range(0,npoints):
#	for j in range(0,npoints):
#	    if abs(grid_x[i][j]) < 1e-02 and abs(grid_y[i][j]) < 1e-02:
#		if abs(grid_x[i][j]) < x_sm and  abs(grid_y[i][j]) < y_sm: 
#		    x_sm = grid_x[i][j]
#		    y_sm = grid_y[i][j]
#		    z_sm = grid_z[i][j]

#    print x_sm
#    print y_sm
#    print z_sm
#    for i in range(0,npoints):
#        for j in range(0,npoints):
#            if grid_x[i][j] == 0.0 and grid_y[i][j] == 0.0:
#		print "hello"
#		print grid_z[i][j]

    print nllmin
#    print minimum 
#    nllmin = 11112234.9491
#    nllmin = 11112234.9635


#    nllmin = 11112234.949976524
#    nllmin = 11112232.932974193


#    nllmin = 11112236.5508
#    nllmin = 11112236.629944898
    from skimage import measure
    allcontours = []
    allimgcontours = []
    for v in values:
        imgcontours = measure.find_contours(grid_z, v + nllmin)
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
       
#    ax.axis('image')
#    ax.set_xticks([])
#    ax.set_yticks([])
#    plt.show() 

    return allcontours,minimum_grid
    
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

