
inf = float("inf")
nan = float("nan")

def disp2dcontour(grid_z,allimgcontours):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.imshow(grid_z, interpolation='nearest', cmap=plt.cm.coolwarm)
    for contours in allimgcontours:
        for n, contour in enumerate(contours):
            ax.plot(contour[:, 1], contour[:, 0], linewidth=2)

    ax.axis('image')
    ax.set_xticks([])
    ax.set_yticks([])
    plt.show()

def disp3dcontour(gridx,gridy,gridz):
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter
    import numpy as np

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # Plot the surface.
    surf = ax.plot_surface(gridx, gridy, gridz, cmap=cm.coolwarm, linewidth=0, antialiased=False)

    # Customize the z axis.
    ax.set_zlim(11112230, 11112245)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    surf.set_clim(11112230,11112245)
    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()

def minfromscans(xvals,yvals,zvals):
    minimum = (inf,inf)
    nllmin = inf
    for i in range(0,len(zvals)):
        if nllmin > zvals[i]:
            minimum = (xvals[i],yvals[i])
            nllmin = zvals[i]
    return minimum,nllmin


def findmergecontours(points1,points2,values,smooth,npoints):
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

    from scipy.interpolate import griddata
    from numpy import mgrid
    grid1_x, grid1_y = mgrid[min(min(xvals1),min(xvals2)):max(max(xvals2),max(xvals1)):complex(npoints), min(min(yvals1),min(yvals2)):max(max(yvals2),max(yvals1)):complex(npoints)]
    grid2_x, grid2_y = mgrid[min(min(xvals1),min(xvals2)):max(max(xvals2),max(xvals1)):complex(npoints), min(min(yvals1),min(yvals2)):max(max(yvals2),max(yvals1)):complex(npoints)]

    grid1_z = griddata(array(keys1),array(zvals1),(grid1_x, grid1_y), method='cubic')
    grid2_z = griddata(array(keys2),array(zvals2),(grid2_x, grid2_y), method='cubic')

    minimum1,nllmin1 = minfromscans(xvals1,yvals1,zvals1)
    minimum2,nllmin2 = minfromscans(xvals2,yvals2,zvals2)

    if nllmin1 < nllmin2 : nllmin, minimum = nllmin1, minimum1
    else                 : nllmin, minimum = nllmin2, minimum2

#perform the evelope on the grid
    for i in range(0,npoints):
        for j in range(0,npoints):
            if grid1_x[i][j] == grid2_x[i][j] and grid1_y[i][j] == grid2_y[i][j]:
                if grid2_z[i][j] != nan and grid2_z[i][j] != nan:
                    grid1_z[i][j] = min(grid2_z[i][j],grid1_z[i][j])

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
                realc.append(((i+0.5)/npoints * (max(max(xvals1),max(xvals2))-min(min(xvals1),min(xvals2))) + min(min(xvals1),min(xvals2)),(j+0.5)/npoints * (max(max(yvals1),max(yvals2))-min(min(yvals1),min(yvals2))) + min(min(yvals1),min(yvals2))))
            if len(realc) > 0:
                if smooth:
                    contours.append(smoothgraph(realc))
                else:
                    contours.append(realc)

        allcontours.append(contours)

    # Display the image and plot all contours found
#    disp2dcontour(allimgcontours)

    return allcontours,minimum

def griddata_gp(xyvals,zvals,grids, options={}):
    from sklearn.gaussian_process import GaussianProcessRegressor
    gp = GaussianProcessRegressor()
    gp.fit(X=xyvals, y=zvals)
    from numpy import column_stack
    grid = column_stack([g.flatten() for g in grids])
    interp = gp.predict(grid)
    return [ [ interp[i+len(grids[0])*j] for i in range(0,len(grids[0])) ] for j in range(0,len(grids[1])) ]

def find_contours_root(xvals,yvals,grid_z,thresholds,smooth,npoints):
    import sys
    from math import isnan
    sys.argv = []
    import ROOT
    from array import array
    ROOT.gROOT.SetBatch(True)
    th2 = ROOT.TH2F("hist","hist",npoints,min(xvals),max(xvals),npoints,min(yvals),max(yvals))
    th2.SetDirectory(0)

    for i in range(0,npoints):
        for j in range(0,npoints):
            th2.SetBinContent(i+1,j+1,grid_z[i][j])

    c = ROOT.TCanvas("c","c",400,400)
    c.cd()
    th2.SetContour(len(thresholds),array("d",thresholds))
    th2.Draw("CONTLIST")
    c.Update()

    root_contours = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
    allcontours = []
    for ic in range(0,root_contours.GetEntries()):
        graphlist = root_contours.At(ic)
        contours = []
        for ig in range(0,graphlist.GetEntries()):
            g = graphlist.At(ig)
            contour = []
            for i in range(g.GetN()):
                x = g.GetX()[i]
                y = g.GetY()[i]
                if not isnan(x) and not isnan(y):
                    contour.append([x,y])
            if len(contour) > 0:
                contours.append(contour)
        if len(contours) > 0:
            allcontours.append(contours)
    return allcontours


def find_contours_skimage(xvals,yvals,grid_z,thresholds,smooth,npoints):
    from math import isnan
    from skimage import measure
    allcontours = []
    for v in thresholds:
        imgcontours = measure.find_contours(grid_z, v)
        contours = []
        for c in imgcontours:
            realc = []
            for i,j in c:
                if isnan(i) or isnan(j): continue
                realc.append(((i+0.5)/npoints * (max(xvals)-min(xvals)) + min(xvals),(j+0.5)/npoints * (max(yvals)-min(yvals)) + min(yvals)))
            if len(realc) > 0:
                if smooth:
                    contours.append(smoothgraph(realc))
                else:
                    contours.append(realc)
        allcontours.append(contours)
    return allcontours

def findcontours(points,values,smooth,npoints,algorithm="ROOT"):
    """find the contours in a 2d graph"""
    from numpy import array
    from math import isnan
    keys = sorted(points.keys())
    xvals = [ p[0] for p in keys ]
    yvals = [ p[1] for p in keys ]
    zvals = [ points[p] for p in keys ]
    xyvals = [ [p[0],p[1]] for p in keys ]
    from scipy.interpolate import griddata
    from numpy import mgrid
    grid_x, grid_y = mgrid[min(xvals):max(xvals):complex(npoints), min(yvals):max(yvals):complex(npoints)]
    grid_z = griddata(array(keys),array(zvals),(grid_x, grid_y), method='linear')
#    grid_z = griddata_gp(array(keys),array(zvals),(grid_x, grid_y))

    minimum,nllmin = minfromscans(xvals,yvals,zvals)
    if algorithm == "ROOT":
        allcontours = find_contours_root(xvals,yvals,grid_z,[ nllmin + v for v in values ],smooth,npoints)
    elif algorithm == "skimage":
        allcontours = find_contours_skimage(xvals,yvals,grid_z,[ nllmin + v for v in values ],smooth,npoints)
    else:
        print("unknown contour finding algorithm '"+algorithm+"'")

    return allcontours,minimum

def findintervals(points,nllval):
    """find the minimum point of a 1d graph and the crossing points with a horizontal line at a given value. returns tuple of central value, lower error and upper error"""
    from scipy.interpolate import PchipInterpolator as interpolate
    from math import isnan
    xvals = [ x for x,y in points ]
    yvals = [ y for x,y in points ]
    xl,xr = min(xvals),max(xvals)
    r = abs(xr - xl)
    xll = xl - 0.1*r
    xrr = xr + 0.1*r
    n = 100
    step = (xrr-xll)/n
    interp = interpolate(xvals, yvals, extrapolate=True)
    from scipy.optimize import ridder as solve
    leftbound = None
    rightbound = None
    intervals = []
    for i in range(1,n):
        left  =  (xll + i    *step)
        right =  (xll + (i-1)*step)
        if (interp(left)-nllval)*(interp(right)-nllval) < 0:
            xcoord = solve(lambda x : interp(x) - nllval,left,right)
            if interp(left) < nllval:
                rightbound = None
                leftbound = xcoord
            else:
                rightbound = xcoord
            if leftbound != None and rightbound != None:
                intervals.append((leftbound,rightbound))
                leftbound = None
                rightbound=None
    return intervals

def findallminima(points,nllmin,minthreshold=0.05):
    """find all the local minima of the nll curve that are no further than the threshold away from the global minimum"""
    intervals = []
    interval = []
    for x,y in points:
        if y<nllmin+minthreshold:
            interval.append((x,y))
        elif len(interval) > 0:
            intervals.append(interval)
            interval = []
    minima = []
    for interval in intervals:
        localmin_y = inf
        localmin_x = None
        for x,y in interval:
            if y<localmin_y:
                localmin_y = y
                localmin_x = x
        minima.append(localmin_x)
    return minima

def findcrossings(points,nllval,minthreshold=0.05):
    """find the minimum point of a 1d graph and the crossing points with a horizontal line at a given value. returns tuple of central value, lower error and upper error"""
    from scipy.interpolate import PchipInterpolator as interpolate
    xvals = [ x for x,y in points ]
    yvals = [ y for x,y in points ]
    i0 = 0
    ymin = inf
    for i in range(0,len(xvals)):
        if yvals[i] < ymin or (yvals[i] < ymin+minthreshold and abs(xvals[i]) < abs(xvals[i0])):
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
    minimum = minimize(lambda v:interp(v[0]), array([min(xvals)]), bounds=[[min(xvals),max(xvals)]])
    miny = minimum.fun
    for y in yvals:
        if y<miny:
            miny=y
    return miny

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

