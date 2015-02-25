import numpy as np
import matplotlib.pylab as plt
import mpl_toolkits.mplot3d.axes3d as p3

import convert
import amr
import IO
import param
import plot
import part

class step:
	def __init__(self,number,folder="data"):
		self.grid_path=	("%s%s/grid/grid.%s"%(folder,"{:05d}".format(number),"{:05d}".format(number) )) 
		self.part_path=	("%s%s/part/part.%s"%(folder,"{:05d}".format(number),"{:05d}".format(number) )) 
		self.star_path=	("%s%s/star/star.%s"%(folder,"{:05d}".format(number),"{:05d}".format(number) )) 
		
		self.param=param.Param(folder)
		
	def grid(self):
		return self.grid_path
	def part(self):
		return self.part_path
	def star(self):
		return self.star_path
					
	def t(self):
		return part.geta(self.part()+".p00000")
		
	def Diag(self,fieldX="field.d",fieldY="rfield.temp",type="hist"):
		
		X_all=amr.allOct(self.grid_path, fieldX)
		Y_all=amr.allOct(self.grid_path, fieldY)
		
		plt.figure()
		if type =="hist":			
			plot.hist2d(X_all.getMap(), Y_all.getMap())
			plt.xlabel('log10 of %s'%fieldX)
			plt.ylabel('log10 of %s'%fieldY)
		elif type=="dot":
			x=X_all.getMap()
			y=Y_all.getMap()
			
			for i in range(X_all.lmin(),int(X_all.lmax()+1)):
				mask=np.where(X_all.getLevel()==i)
				plt.loglog(x[mask],y[mask],'.',label='level%d'%i)
			
			plt.legend()
			plt.xlabel('%s'%fieldX)
			plt.ylabel('%s'%fieldY)
		plt.show()


	def diagXm_Xv(self):
		X_all=amr.allOct(self.grid_path, "field.d").getMap()
		Y_all=amr.allOct(self.grid_path, "rfield.xion").getMap()
		
		plt.figure()
		plt.xlabel('log10 of Xion_v')
		plt.ylabel('log10 of Xion_m')
		
		X_all+=1.
	
		tmp=Y_all
		Y_all = X_all*Y_all 
		X_all=tmp
		plot.hist2d(X_all, Y_all)
		

	def phaseDiag(self):		
		X_all=amr.allOct(self.grid_path, "field.d").getMap()
		Y_all=amr.allOct(self.grid_path, "rfield.xion").getMap()
		
		plt.figure()
		plt.xlabel(r'$log10 \left( 1+\delta \right)$')
		plt.ylabel(r'$log10\left(\left(1+\delta\right) \left( \frac{{x_{HII}}^2}{x_{HI}} \right) \right)$')
					
		mask = np.where(Y_all!=1.)
		X_all = X_all[mask]
		Y_all = Y_all[mask]
		X_all+=1.
		Y_all = X_all * np.power(Y_all,2)/(1.-Y_all)
		
	
		plot.hist2d(X_all, Y_all)
		
		#plt.xlim(0,4.5)
		#plt.ylim(-10,10)
		
	def phaseMass(self):
		AMR=amr.allOct(self.grid_path, "field.d")
		rho=AMR.getMap()
		l=AMR.getLevel()
		
		dx = np.power(2.,-l)
		dv = np.power(dx,3.)
		
		nstar,a,star = part.read(self.star_path)

		m=rho*dv
		
		"""
		for i in range(AMR.lmin(),int(AMR.lmax()+1)):
			mask=np.where(l==i)
			plt.loglog(rho[mask],m[mask],'.',label='level%d'%i)
		plt.legend()
		"""
		
		plot.hist2d(rho,m)
		
		m1 = np.log10(np.unique(star.mass)[0])
		plt.axhline(y=m1)

		rho1 = 86
		plt.axvline(x=np.log10(rho1))
