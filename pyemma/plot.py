import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from mayavi import mlab
from scipy.interpolate import griddata
from scipy.ndimage.interpolation import map_coordinates



from convert import *
from amr import *
import physique


def hist2d(X_all, Y_all):
	n=1024	
	xedges=np.logspace(np.min(X_all),np.max(X_all),n)
	yedges=np.logspace(np.min(Y_all),np.max(Y_all),n)
	
	x=np.log10(X_all)
	y=np.log10(Y_all)
	
	H,yedges,xedges=np.histogram2d(y,x,bins=n)
	X, Y = np.meshgrid(xedges, yedges)
	H = np.ma.masked_where(H==0,H) 

	plt.pcolormesh(X, Y, np.log(H))
	
	cb=plt.colorbar()
	cb.set_label('log10 of number of occurence')

	plt.xlim(xedges[0], xedges[-1])
	plt.ylim(yedges[0], yedges[-1])

		
	

def getSlice(filename,level,force=0, nproc=0, field="density", xmin=0, xmax=-1, ymin=0, ymax=-1, zmin=0, zmax=-1, log=True):
	N = pow(2,level)
	if xmax == -1 :
		xmax = N
	if ymax == -1 :
		ymax = N
	if zmax == -1 :
		zmax = N
		
	if xmin<0:
		xmin=0	
	if xmax>1:
		xmax=1	

	if ymin<0:
		ymin=0	
	if ymax>1:
		ymax=1
		
	if zmin<0:
		zmin=0	
	if zmax>1:
		zmax=1		
	
		
	if nproc==0:
		nproc=IO.getNproc(filename)
		
	convert.oct2grid(filename,level, force, nproc, field, xmin, xmax, ymin, ymax, zmin, zmax, proj=3)

	data = cube(filename.replace("grid.","cube"+ str(level)+ ".") + "." + field)

	print "Z=%.4f"%data.getZ()
	data = data.getData()

	if log:
		data = np.log10(data)
	
	data = np.max(data,axis=0)
	#data = np.mean(data,axis=0)
	return data

def slice(filename,level,force=0, nproc=0, field="density", xmin=0, xmax=-1, ymin=0, ymax=-1, zmin=0, zmax=-1, log=False):
	if nproc==0:
		nproc=IO.getNproc(filename)
	data = getSlice(filename,level,force, nproc, field, xmin, xmax, ymin, ymax, zmin, zmax, log)
	
	print np.mean(data)
	print np.max(data)
	
	N = pow(2,level)
	if xmax == -1 :
		xmax = N
	if ymax == -1 :
		ymax = N
	if zmax == -1 :
		zmax = N
	plt.clf()
	plt.imshow(data, interpolation='nearest',extent=(xmin,xmax,ymin,ymax),origin='lower' )
	plt.colorbar()
	plt.show(block=False)

def diffslice(filename1,filename2,level,force=0, nproc=0, field="density", xmin=0, xmax=-1, ymin=0, ymax=-1, zmin=0, zmax=-1, log=False):
	if nproc==0:
		nproc=IO.getNproc(filename1)
	data1 = getSlice(filename1,level,force, nproc, field, xmin, xmax, ymin, ymax, zmin, zmax, log)
	data2 = getSlice(filename2,level,force, nproc, field, xmin, xmax, ymin, ymax, zmin, zmax, log)

	data=np.subtract(data1,data2)

	"""	
	if log:
		data = np.log10(np.abs(data))
	"""
	plt.clf()
	plt.imshow(data, interpolation='nearest',extent=(xmin,-xmax,ymin,-ymax),origin='lower' )
	plt.colorbar()
	plt.show(block=False)


def readMovie(file):
	m=[]
	
	f = open(file, "rb")
	x 	= np.fromfile(f, dtype=np.int32  ,count=1)[0]
	y 	= np.fromfile(f, dtype=np.int32  ,count=1)[0]
	a 	= np.fromfile(f, dtype=np.float32  ,count=1)[0]			
	m.append(np.fromfile(f, dtype=np.float32  ,count=x*y))		#pot
	m.append(np.fromfile(f, dtype=np.float32  ,count=x*y))		#den
	m.append(np.fromfile(f, dtype=np.float32  ,count=x*y))		#X
	m.append(np.fromfile(f, dtype=np.float32  ,count=x*y))		#temp	
	f.close()

	return x,y,a,m

def movie(folder ="data/movie/", field=1):

	files = np.sort(os.listdir(folder))
	i,j=0,0
	
	nsub = 1
	fig = plt.figure(figsize=(10,10))
	
	if field == 0:
		fieldname = "gdata.p"
	if field == 1:
		fieldname = "field.d"
	if field == 2:
		fieldname = "rfield.xion"
	if field == 3:
		fieldname = "rfield.temp"
	
	try:
		os.mkdir(folder+"../img")
	except  OSError:
		print "folder exists"
	
	for file in files:		
		i+=1
		if  i%nsub == 0 :
			j+=1		
			image_file = "img/"+file.replace("movie",fieldname)+".png"
			
			if not os.path.isfile(image_file):

				print "%d / %d"%(j,len(files)/nsub)

				x,y,a,m = readMovie(folder+file)

				t = physique.a2t(a)
				
				data1 = m[field].reshape(x,y)
	
				
				plt.axis("off")
				plt.text( 2, 10, 'size = 2 Mpc.h^-1')
				plt.text( 2, 20, 'z = %s '%"{:.2f}".format(1./a-1.))
				plt.text( 2, 30, 't = %s Myrs'%"{:.2f}".format(t/1e6))
				
				fig=plt.imshow(np.log(data1[::8,::8]), interpolation='nearest')
				fig.axes.get_xaxis().set_visible(False)
				fig.axes.get_yaxis().set_visible(False)
				
				plt.savefig(folder+"../"+image_file, bbox_inches='tight', pad_inches=0)
				plt.clf()


def part(parts, level=1):
	s=pow(2.0,level)

	plt.plot(parts.x*s,parts.y*s,'.')
	#fig = plt.figure()
	#ax = fig.add_subplot(111, projection='3d')
	#ax.scatter(parts.x*s,parts.y*s,parts.z*s,'.')
	plt.show(block=False)

def interp3D(x, y, z, v, xi, yi, zi, **kwargs):
    """Sample a 3D array "v" with pixel corner locations at "x","y","z" at the
    points in "xi", "yi", "zi" using linear interpolation. Additional kwargs
    are passed on to ``scipy.ndimage.map_coordinates``."""
    def index_coords(corner_locs, interp_locs):
        index = np.arange(len(corner_locs))
        if np.all(np.diff(corner_locs) < 0):
            corner_locs, index = corner_locs[::-1], index[::-1]
        return np.interp(interp_locs, corner_locs, index)

    orig_shape = np.asarray(xi).shape
    xi, yi, zi = np.atleast_1d(xi, yi, zi)
    for arr in [xi, yi, zi]:
        arr.shape = -1

    output = np.empty(xi.shape, dtype=float)
    coords = [index_coords(*item) for item in zip([x, y, z], [xi, yi, zi])]

    map_coordinates(v, coords, order=1, output=output, **kwargs)

    return output.reshape(orig_shape)

def p3D(file,level,field,nproc=0,force=0,xmin=0, xmax=-1, ymin=0, ymax=-1, zmin=0, zmax=-1, log=False):
	if nproc==0:
		nproc=IO.getNproc(file)
		
	oct2grid(file,level, force, nproc, field, xmin, xmax, ymin, ymax, zmin, zmax, proj=0)
	data = cube(file.replace("grid.","cube"+ str(level)+ ".") + "." + field)
	data = data.getData()
	
	
	n=128
	x1 = np.linspace(0,1,n)
	y1 = np.linspace(0,1,n)
	z1 = np.linspace(0,1,n)
	
	n=256
	x2 = np.linspace(0,1,n)
	y2 = np.linspace(0,1,n)
	z2 = np.linspace(0,1,n)
	x2,y2,z2=np.meshgrid(x2,y2,z2)
	
	data = interp3D(x1,y1,z1, data, x2,y2,z2)
	
	contours=[1e4]
	mlab.contour3d(data,opacity=0.1, contours=contours)
	mlab.savefig('test_256.obj')

	mlab.show()
	
