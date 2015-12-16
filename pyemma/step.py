import numpy as np
import matplotlib.pylab as plt
import mpl_toolkits.mplot3d.axes3d as p3

import sfr
import convert
import amr
import param
import plot
import part
import physique
import sys,os

class step:
	def __init__(self,number,folder="data", injection_type=None):
		
		self.number=number
		format_number = "{:05d}".format(number)

		self.run_folder = folder
		self.step_folder = ("%s%s/"%(folder,format_number))

		self.grid_folder=("%sgrid/"%(self.step_folder))
		self.part_folder=("%spart/"%(self.step_folder))
		self.star_folder=("%sstar/"%(self.step_folder))
		
		self.grid_name=("grid.%s"%(format_number))
		self.part_name=("part.%s"%(format_number))
		self.star_name=("star.%s"%(format_number))
		self.halo_name=("halo.%s"%(format_number))
		
		self.grid_path=("%s%s"%(self.grid_folder,self.grid_name))
		self.part_path=("%s%s"%(self.part_folder,self.part_name))
		self.star_path=("%s%s"%(self.star_folder,self.star_name))
		self.halo_path=("%s%s"%(self.part_folder,self.halo_name))		
#______________________________________________________________________#

		self.param=param.Param(self.step_folder,number)		
#______________________________________________________________________#

		try:
			self.a = part.geta(self.part_path+".p00000")		
			print self.a
			self.t = physique.a2t(self.a)
			self.z = physique.a2z(self.a)
		except IOError:
			try:
				self.a = amr.geta(self.grid_path+".p00000")
				print amr.geta(self.grid_path+".p00000")
				self.t = physique.a2t(self.a)
				self.z = physique.a2z(self.a)
			except IOError:
				print "can't determine the current time"
				self.a=-1
				self.t=-1
				self.z=-1
				
			
#______________________________________________________________________#
		
		try:
			files = os.listdir(self.grid_folder)
			name = self.grid_name
		except  OSError:
			print "can't found grid trying part"		
			try:			
				files = os.listdir(self.part_folder)
				name = self.part_name
				print "part found"
			except  OSError:
				print "can't found part : abort"
				sys.exit(0)
			
		self.nproc =0
		for file in files:
			if name in file :
				self.nproc += 1
				
		self.injection_type = injection_type

########################################################################

		
	def pdf(self,field, **kwargs):
		AMR=amr.allOct(self.grid_path,field)
		rho=AMR.getMap()
		l=AMR.getLevel()
		
		dx = np.power(0.5,l)
		dv = np.power(dx,3.)
		
		Nbins= 512
	
		rho_min = np.min(rho)
		rho_max = np.max(rho)
		bin_edges = np.logspace(np.log10(rho_min),np.log10(rho_max),Nbins+1)
		y,bin_edges=np.histogram(rho, bins=bin_edges, weights=dv)
		x=(bin_edges[1:]+bin_edges[:-1])/2
		print x,y
		
		dx = np.diff(bin_edges)
		plt.loglog(x,y/dx,'.', **kwargs)
	
		plt.axvline(x=25, ymin=0, ymax =1, color='k')

	def sfr(self, **kwarg):
		sfr.fromSnap(self, **kwarg)		
		plt.ylim(-5,0)
		sfr.observation()

	def xion(self):
		x,y = self.param.avg.getxy("z","mean_xion")
		
		mask =  np.where( (y>0) & (y<1) )
		
		x=x[mask]
		y=y[mask]
		
		plt.semilogy(x,y, label="x")
		plt.semilogy(x,1.-y, label="1-x")
			
		"""
		plt.plot(x,y, label="x")
		plt.plot(x,1.-y, label="1-x")
		"""
		plt.legend()
		
		er=[[0.5],[1]]
		plt.errorbar(6., 0.01, xerr=er, ls='none', fmt='ko')
		plt.errorbar(6., 0.5, xerr=er, ls='none', fmt='ko')
	
		plt.xlim(4,10)
		plt.xlabel('z')
		plt.ylabel(r'$ionisation fraction$')

	def bilan(self):
		plt.figure()	
		
		src_int = self.param.run.get()["src_int_or_fesc"]
		overdensity_cond = self.param.run.get()["overdensity_cond"]
		efficiency = self.param.run.get()["efficiency"]
		

		
		#plt.subplot(121)
		self.sfr()
		plt.xlim(0,14)
		
		x=2
		y=-4
		dy = -0.25
		plt.text(2, y, 'src_int=%s'%src_int)
		plt.text(2, y+dy, 'overdensity_cond=%s'%overdensity_cond)
		plt.text(2, y+2*dy, 'efficiency=%s'%efficiency)
		
		#plt.subplot(122)
		plt.figure()
		self.xion()

		

		
########################################################################

def freefall(step):
	def get_thresh(p):
		cond = float(step.param.run.get()["overdensity_cond"])
		info = step.param.info.get()	
		thresh = cond *info["ob"]/info["om"]
		print thresh		
		return thresh


	def get_field(step, field):		
		AMR=amr.allOct(step.grid_path,field)
		m=AMR.getMap()
		l=AMR.getLevel()	
		return l,m
		
		
		

	#plt.figure()
	cst = physique.Constantes()
	thresh=get_thresh(step.param)

	
	l,rho_g= get_field(step, "field.d") 
	#rho_g_phy = physique.code2d(step,rho_g/thresh)
	rho_g_phy = physique.code2d(step,rho_g)
	
	
	l,rho_m= get_field(step, "den")
	rho_m += 1.
	rho_m_phy = physique.code2d(step,rho_m/thresh)
	#rho_m_phy = physique.code2d(step,rho_m)
	print np.mean(rho_g)
	print np.min(rho_g)
	print np.max(rho_g)
	
	
	l,rho_p= get_field(step, "density")
	rho_p_phy = physique.code2d(step,rho_p)
	print rho_p


	dx = np.power(0.5,l)
	dv = np.power(dx,3.)
	

	mask = np.where((rho_p!=0) & (l==np.max(l)-2) )
	mask = np.where((rho_p!=0)) 
	#mask = np.where( (rho_m-rho_g) )
	
	#plot.hist2d(rho_g[mask] ,rho_m[mask])#, weights= dv[mask])
	#plot.hist2d(rho_g[mask] ,rho_m[mask]-rho_g[mask])#, weights= dv[mask])
	#plot.hist2d(rho_g ,rho_m)#, weights= dv[mask])

	
	"""
	l,p= get_field(step, "field.p")	
	cs = np.sqrt(5./3.*p/rho_g);
	tj = physique.code2t(step,cs/dx)
	"""
	


	tff = np.sqrt(3*np.pi/(32*cst.G*(rho_m_phy))) # s

	#plot.hist2d(rho_g, tff/cst.yr )
	
	#mask = np.where(l==np.max(l)-3)
	#plot.hist2d(rho_g[mask], rho_m[mask], dv[mask] )
	#plot.hist2d(rho_g, rho_m, dv )
	
	eps = float(step.param.run.get()["efficiency"])
	sfr = eps*rho_g_phy/tff; print sfr # kg/m3/s
	sfr*= (cst.Parsec*1e6)**3; print sfr #Mpc3
	sfr*=cst.yr; print sfr #yr
	sfr/=cst.M0; print sfr #Mo
	#sfr*=tj/tff

	print np.mean(sfr)
	print np.min(sfr)
	print np.max(sfr)
	
	#l0 sfr 6.65163048e-05
	mask = np.where(sfr)
	plot.hist2d(rho_g[mask],sfr[mask])#,weights=dv[mask])
	
	"""
	Vbox= (step.param.info.get()["box_size_Mpc/h"] * (step.param.info.get()["H0"]/100) )**3
	x,y = plot.hist1d(rho_g,sfr*(dv*Vbox))
	y = np.cumsum( y[::-1])[::-1]
	plt.loglog(x,y)
	"""
	
	
	#x,y = plot.hist1d(rho_m, dv )
	#plt.loglog(x,y)
	
	
	'''
	thresh_sfr = 5.-06
		
	thresh_sfr_proj = np.where(np.diff(np.sign(y-thresh_sfr)))[0]
	xline = x[thresh_sfr_proj]
	print xline
	plt.axvline(x=xline, color='k')
	'''
	
	"""
	thresh_proj = np.where(np.diff(np.sign(x-thresh)))[0]
	yline = y[thresh_proj+1]
	print yline
	"""

	#plt.axhline(y=yline, color='k')
	#plt.axvline(x=np.log10(thresh), color='k')
	plt.axvline(x=thresh, color='k')
	
	
	#plt.axhline(y=yline, color='r', linewidth=3)
	#plt.axvline(x=xline, color='r', linewidth=3)
	
	
	#plt.axvline(x=thresh, color='r', linewidth=3)
	



	plt.xlabel(r'log10 gaz overdensity')
	
	#plt.ylabel(r'log10 free fall time (yr)')
	#plt.ylabel(r'log10 SFR ')
	plt.ylabel(r'log10 collisionless matter overdensity')
	
	
def test_sfr():
	
	rho = [3.9391664, 7.27806527, 24.11086891]
	#rho = [5.10295449e-06,	4.81801098e-05,	0.00024406]
	x= range(4)[1:]

	plt.xscale('log')
	plt.yscale('log')
	
	plt.plot(x,rho)
	
		
	#rho = np.log10(rho)
	#x = np.log10(x)
	A= np.polyfit(x,rho,2)
	print A
	x2= np.arange(np.min(x),np.max(x),1e-2)
	#y2= A[0]*x2+A[1]
	y2= A[0]*x2**2 + A[1]*x2+ A[2]
	plt.plot(x2,y2)
	

