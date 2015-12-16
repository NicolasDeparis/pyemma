import os
import numpy as np
import field

class Grid:
	def __init__(self,number,folder):
		"""			
			Create a grid object from EMMA output 					
		"""
		
		self.number=number
		
		path = "%s%05d/"%(folder,number)
		
		for cur_folder in  os.listdir(path):
			if "grid_" in cur_folder:
				
				key=cur_folder.strip("grid_").replace(".","_")
				val= field.Field(folder,number,cur_folder)
				setattr(self,key,val)
		
def get_cube(x,y,z,l,map,level,type, xmin=0.,xmax=1.,ymin=0.,ymax=1.,zmin=0.,zmax=1.):

    def project(x,dl,type):
        n=2**dl
        if type=="2d":
            return np.repeat(np.repeat(x,n,axis=0),n,axis=1)*n
        elif type=="3d":
            return np.repeat(np.repeat(np.repeat(x,n,axis=0),n,axis=1),n,axis=2)*n
        else:
            print("type should be 2d or 3d")

    def get_cube_level(x,y,z,l,map,level,lmax, type):
        if level==lmax:
            mask=np.where(l>=level)[0]
            dv_lmax=0.5**(3*lmax)
            dv_lcur=0.5**(3*l[mask])
            w=map[mask] * dv_lcur/dv_lmax
        else:
            mask=np.where(l==level)[0]
            w=map[mask]

        if type=="2d":            
            x=x[mask]
            y=y[mask]
            r= np.array( (x,y)) .transpose()
            bins_1d=np.arange(0.,1.+0.5**level,0.5**level)
            bin_edges = [bins_1d,bins_1d]
            h,_=np.histogramdd(r, weights=w,bins=bin_edges)
            return h
        elif type=="3d":
            x=x[mask]
            y=y[mask]
            z=z[mask]
            r= np.array( (x,y,z)) .transpose()
            bins_1d=np.arange(0.,1.+0.5**level,0.5**level)        
            bin_edges = [bins_1d,bins_1d,bins_1d]
            h,_=np.histogramdd(r, weights=w,bins=bin_edges)
            return h
        else:
            print("type should be 2d or 3d")

    def recursive_cube(x,y,z,l,map,level,lmax,type):
        if level<lmax:            
            cur_map = project(get_cube_level(x,y,z,l,map,level,lmax,type),lmax-level,type)
            cur_map += recursive_cube(x,y,z,l,map,level+1,lmax,type)
            return cur_map
        else:
            return get_cube_level(x,y,z,l,map,level,lmax,type)

    mask=np.where(  (x>=xmin) & (x<xmax) &
                    (y>=ymin) & (y<ymax) &
                    (z>=zmin) & (z<zmax) )

    x = x[mask]
    y = y[mask]
    z = z[mask]
    l = l[mask]
    map = map[mask]

    lmin=min(np.min(l),level)
    print(lmin)
    lmax=np.max(l)
    if level>lmax:
        level=lmax
        print("level max = %d"%level)
    
    return recursive_cube(x,y,z,l,map,lmin,level,type)/(2**level)
