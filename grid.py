# coding: utf-8

import numpy as np

def get_cube(x,y,z,l,map,level,type, xmin=0.,xmax=1.,ymin=0.,ymax=1.,zmin=0.,zmax=1.):
    """
    return a cube if type=="3d" or a slice if type=="2d"
    """

    def prolong(x,dl,type):
        """
        increase the resolution of the map by a direct injection prolongation methode
        """

        n=2**dl
        if type=="2d":
            return np.repeat(np.repeat(x,n,axis=0),n,axis=1)*n
        elif type=="3d":
            return np.repeat(np.repeat(np.repeat(x,n,axis=0),n,axis=1),n,axis=2)
        else:
            print("type should be 2d or 3d")

    def get_cube_level(x,y,z,l,map,level,lmax, type):
        """
        project the data of current level on a grid
        """

        if level==lmax:
            mask=np.where(l>=level)[0]
            w=map[mask] * np.power(0.5,3*(l[mask]-lmax))
        else:
            mask=np.where(l==level)[0]
            w=map[mask]

        bins_1d=np.arange(0.,1.+0.5**level,0.5**level)

        if type=="2d":
            """
            project a map
            """
            x=x[mask]
            y=y[mask]
            r= np.array((x,y)).transpose()
            bin_edges = [bins_1d,bins_1d]

        elif type=="3d":
            """
            project a cube
            """
            x=x[mask]
            y=y[mask]
            z=z[mask]
            r= np.array( (x,y,z)) .transpose()
            bin_edges = [bins_1d,bins_1d,bins_1d]

        else:
            print("type should be 2d or 3d")

        h,_=np.histogramdd(r, weights=w, bins=bin_edges)
        return h


    def recursive_cube(x,y,z,l,map,level,lmax,type):
        """
        get the grid recursively
        """
        if level<lmax:
            cur_map = prolong(get_cube_level(x,y,z,l,map,level,lmax,type),lmax-level,type)
            cur_map += recursive_cube(x,y,z,l,map,level+1,lmax,type)
            return cur_map
        else:
            return get_cube_level(x,y,z,l,map,level,lmax,type)

    #control level min
    lmin=np.min(l)
    if level<lmin:
        level=lmin
        print("level min = %d"%level)

    #control level max
    lmax=np.max(l)
    if level>lmax:
        level=lmax
        print("level max = %d"%level)

    ngrid=2**level

    #get coordinate in grid space (0,ngrid)
    _xmin=int(xmin*ngrid)
    _xmax=int(xmax*ngrid)
    _ymin=int(ymin*ngrid)
    _ymax=int(ymax*ngrid)
    _zmin=int(zmin*ngrid)
    _zmax=int(zmax*ngrid)

    #get coordinate in grid space (0,1)
    xmin= _xmin/ngrid
    xmax= _xmax/ngrid
    ymin= _ymin/ngrid
    ymax= _ymax/ngrid
    zmin= _zmin/ngrid
    zmax= _zmax/ngrid

    #get only data of interest
    mask=np.where(  (x>=xmin) & (x<xmax) &
                    (y>=ymin) & (y<ymax) &
                    (z>=zmin) & (z<zmax) )
    x = x[mask]
    y = y[mask]
    z = z[mask]
    l = l[mask]
    map = map[mask]

    #projection
    cube=recursive_cube(x,y,z,l,map,lmin,level,type)

    if type=="2d":
        return cube[_xmin:_xmax,_ymin:_ymax]/(_zmax-_zmin)
    if type=="3d":
        return cube[_xmin:_xmax,_ymin:_ymax,_zmin:_zmax]
