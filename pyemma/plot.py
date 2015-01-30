import matplotlib.pylab as plt
import numpy as np

from convert import *
from amr import *
import physique

def getSlice(filename,level,force=0, nproc=0, field="density", xmin=0, xmax=-1, ymin=0, ymax=-1, zmin=0, zmax=-1, log=False):
	N = pow(2,level)
	if xmax == -1 :
		xmax = N
	if ymax == -1 :
		ymax = N
	if zmax == -1 :
		zmax = N
		
	if nproc==0:
		nproc=IO.getNproc(filename)
		
	oct2grid(filename,level, force, nproc, field, xmin, xmax, ymin, ymax, zmin, zmax, proj=3)

	data = cube(filename.replace("grid","cube"+ str(level)) + "." + field)

	print data.geta()
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

	plt.clf()
	plt.imshow(data, interpolation='nearest',extent=(xmin,xmax,ymin,ymax) )
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
	plt.imshow(data, interpolation='nearest',extent=(xmin,xmax,ymin,ymax) )
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

def movie(folder ="data/movie/", field=2):

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
				plt.text( 2, 10, 'size = 4 Mpc.h^-1')
				plt.text( 2, 20, 'z = %s '%"{:.2f}".format(1./a-1.))
				plt.text( 2, 30, 't = %s Myrs'%"{:.2f}".format(t/1e6))
				
				fig=plt.imshow(np.log(data1), interpolation='nearest')
				fig.axes.get_xaxis().set_visible(False)
				fig.axes.get_yaxis().set_visible(False)
				
				plt.savefig(image_file, bbox_inches='tight', pad_inches=0)
				plt.clf()

				
