
# coding: utf-8

# In[139]:

import numpy as np
import pickle
import os
from scipy import spatial


# In[3]:

class Fof:
    def __init__(self,folder,stepnum):
        self.exec_folder = "/home/deparis/Emma/utils/pfof/"
        
        self.folder=folder
        self.stepnum=stepnum               
        self.path="%s%05d/halo/"%(self.folder,self.stepnum)
        
        self._isLoaded_masst = False
        self._isLoaded_struct = False
        
    def __getattr__(self, name):
                
        if name == 'R200':            
            self.get_R200()
            return self.__getattribute__(name)
        
        elif name == 'lmin':            
            self.getLmin()
            return self.__getattribute__(name)

        elif name == 'ob':            
            self.getCosmo()
            return self.__getattribute__(name)

        elif name == 'om':
            self.getCosmo()
            return self.__getattribute__(name)
        
        elif name == 'nfoftot':            
            self.getNfofTot()
            return self.__getattribute__(name)
                
        elif name == 'x':            
            self.read_masst()
            return self.__getattribute__(name)
        
        elif name == 'y':
            self.read_masst()
            return self.__getattribute__(name)
        
        elif name == 'z':            
            self.read_masst()
            return self.__getattribute__(name)
        
        elif name == 'idx':
            self.read_masst()
            return self.__getattribute__(name)
        
        elif name == 'npart':
            self.read_masst()
            return self.__getattribute__(name)
        
        elif name == 'part_n':
            self.read_masst()
            return self.__getattribute__(name)
              
        elif name == 'part_pos':
            self.read_struct()
            return self.__getattribute__(name)
        
        elif name == 'part_vel':
            self.read_struct()
            return self.__getattribute__(name)
        
        elif name == 'part_id':
            self.read_struct()
            return self.__getattribute__(name)
        
        elif name == 'nproc_sim':
            self.getNproc_sim()
            return self.__getattribute__(name)
        
        elif name == 'ncpu':
            self.getNcpu()
            return self.__getattribute__(name)
        
        elif name == 'inertia_eig_val':
            self.get_inertia_tensor()
            return self.__getattribute__(name)
        elif name == 'inertia_eig_vec':
            self.get_inertia_tensor()
            return self.__getattribute__(name)
        
        elif name == 'part_mass_fine':
            self.get_part_mass_fine()
            return self.__getattribute__(name)
        
        
        else:
            print("ERROR : Can't automatically determine %s"%name)
            raise AttributeError
            return
            
#         if name == 'cells':
#             print("ERROR : call get_cells() first")
#             return
        
#         if name == 'cells_fine':
#             print("ERROR : call get_cells_fine() first")
#             return
            
    def write_infosim(self):
        self.getLmin()
        with open("%sinfosim.txt"%self.exec_folder,"w") as f:
            f.write("Npro %d\n"%self.nproc_sim)
            f.write("Ndim 3\n")
            f.write("lmin %d\n"%self.lmin)
                    
    def write_fofin(self):        
        with open("%sfof.in"%self.exec_folder,"w") as f:
            f.write("%s%05d/halo/halo\n"%(self.folder,self.stepnum))
            f.write("EMM\n")
            f.write("%s\n"%self.folder)
            f.write("dummyname\n")
            f.write("dummyname\n")
            f.write("0 !group size should be set to zero\n")
            f.write("0.2 !percolation parameter\n")
            f.write("%d !snapnum\n"%self.stepnum)
            f.write("10 !min mass\n")
            f.write("100000000 !mass max\n")
            f.write(".false. !star ?\n")
            f.write(".false. ! metal ?\n")
            f.write(".false. !output apres lecture\n")
            f.write(".true.  !detection des structures\n")
            f.write(".false. !lire depuis cube\n")
            f.write(".true. !timings\n")

    def gen(self,nproc_fof):
        
        self.ncpu = nproc_fof
        self.write_fofin()
        self.write_infosim()
        
        curdir=os.getcwd()
        os.chdir(self.exec_folder)
        
        folder_name ="%s%05d/halo/"%(self.folder, self.stepnum) 
        try : 
            os.mkdir(folder_name)
        except OSError:
            print( "Can't open %s"%folder_name)
            pass
        
        commande= "mpirun -np %d ./fof"%(nproc_fof)
        print ("executing %s"%commande)
        os.system(commande)
        os.chdir(curdir)        
        
    def getNcpu(self):    
        """
        count the number of halo_masst* files in path
        """
        try:     
            self.ncpu=0
            files=os.listdir(self.path)            
            for f in files:
                if "halo_masst" in f:
                    self.ncpu+=1            
        except FileNotFoundError:
            print("ERROR : can't determine nproc in \n%s"%(self.path))
    
    def getNproc_sim(self): 
        """
        Find coarse level by reading param.run
        """
        from pyemma import io
        info = io.Info(self.folder+"../")
        self.nproc_sim=info.nproc
        
    def getLmin(self):    
        """
        Find coarse level by reading param.run
        """
        from pyemma import io
        runparam = io.RunParam(self.folder+"../")
        self.lmin=runparam.level_coarse
        
    def getCosmo(self):
        """
        Find cosmo by reading param.info
        """
        from pyemma import io
        info = io.Info(self.path+"../../../")
        self.ob=info.ob
        self.om=info.om
        self.unit_mass=info.unit_mass
                    
    def getNfofTot(self):
        """
        get the total number oh halos by reading halo_masst headers
        """
        print("Getting nfoftot")
        nfoftot=0
        for icpu in range(self.ncpu):            
            with open("%s%05d/halo/halo_masst_%05d"%(self.folder,self.stepnum, icpu)) as f:
                dummy=np.fromfile(f,dtype=np.int32,count=1)[0]
                nfoftot+=np.fromfile(f,dtype=np.int32,count=1)[0]
        self.nfoftot=nfoftot
        
####################################################################################################################
# READERS
####################################################################################################################
        
    def read_masst(self):
        print("Reading fof masst")
        
        self.idx=np.zeros(self.nfoftot)
        self.npart=np.zeros(self.nfoftot)
        self.x=np.zeros(self.nfoftot)
        self.y=np.zeros(self.nfoftot)
        self.z=np.zeros(self.nfoftot)

        offset=0
        for icpu in range(self.ncpu):
            filename= "%s%05d/halo/halo_masst_%05d"%(self.folder,self.stepnum, icpu)
            with open(filename) as f:
                dummy=np.fromfile(f,dtype=np.int32,count=1)[0]
                nfof=np.fromfile(f,dtype=np.int32,count=1)[0]
                dummy=np.fromfile(f,dtype=np.int32,count=1)[0]        

                for i in range(nfof):
                    dummy=np.fromfile(f,dtype=np.int32,count=1)[0]
                    self.idx[i+offset]=np.fromfile(f,dtype=np.int64,count=1)[0]
                    self.npart[i+offset]=np.fromfile(f,dtype=np.int32,count=1)[0]
                    self.x[i+offset]=np.fromfile(f,dtype=np.float32,count=1)[0]
                    self.y[i+offset]=np.fromfile(f,dtype=np.float32,count=1)[0]
                    self.z[i+offset]=np.fromfile(f,dtype=np.float32,count=1)[0] 
                    dummy=np.fromfile(f,dtype=np.int32,count=1)[0]

                offset+=nfof
    
    def read_struct(self):
        print("Reading fof struct")
        
        nhalo=0

        for icpu in range(self.ncpu):
            filename= "%s%05d/halo/halo_strct_%05d"%(self.folder,self.stepnum, icpu)
            with open(filename) as f:
                dummy=np.fromfile(f,dtype=np.int32,count=1)
                nhalo+=np.fromfile(f,dtype=np.int32,count=1)[0]

        self.part_n = np.zeros(nhalo, dtype=np.int32)
        self.part_pos = np.empty(nhalo, dtype=np.object)
        self.part_vel = np.empty(nhalo, dtype=np.object)
        self.part_id = np.empty(nhalo, dtype=np.object)

        cur_halo=0
        for icpu in range(self.ncpu):
            filename= "%s%05d/halo/halo_strct_%05d"%(self.folder,self.stepnum, icpu)
            with open(filename) as f:
                dummy=np.fromfile(f,dtype=np.int32,count=1)
                cur_nhalo=np.fromfile(f,dtype=np.int32,count=1)[0]
                dummy=np.fromfile(f,dtype=np.int32,count=1)

                for i in range(cur_nhalo):
                    dummy=np.fromfile(f,dtype=np.int32,count=1)
                    self.part_n[cur_halo]=np.fromfile(f,dtype=np.int32,count=1)[0] 
                    dummy=np.fromfile(f,dtype=np.int32,count=1)

                    dummy=np.fromfile(f,dtype=np.int32,count=1)
                    self.part_pos[cur_halo]=np.fromfile(f,dtype=np.float32,count=3*self.part_n[cur_halo])
                    dummy=np.fromfile(f,dtype=np.int32,count=1)

                    dummy=np.fromfile(f,dtype=np.int32,count=1)
                    self.part_vel[cur_halo]=np.fromfile(f,dtype=np.float32,count=3*self.part_n[cur_halo])
                    dummy=np.fromfile(f,dtype=np.int32,count=1)

                    dummy=np.fromfile(f,dtype=np.int32,count=1)
                    self.part_id[cur_halo]=np.fromfile(f,dtype=np.int64,count=self.part_n[cur_halo])
                    dummy=np.fromfile(f,dtype=np.int32,count=1)

                    cur_halo+=1

####################################################################################################################
# Finding
####################################################################################################################

    def get_R200(self, force=0):        
        name = self.path+"R200"
        if os.path.isfile(name) and not force:
            print("Reading %s"%name)
            with open(name, 'rb') as input:                
                self.R200 = pickle.load(input)
            return

        nptot=2**(3*self.lmin)
        part_mass=(1.-self.ob/self.om)/(nptot)

        self.R200 = np.zeros(self.nfoftot, dtype=np.float32)
        for iHalo in range(self.nfoftot):                
            self.R200[iHalo] = np.power(3.*self.npart[iHalo]*part_mass/(200.*4.*np.pi),1./3.)

        with open(name, 'wb') as output:
            pickle.dump(self.R200, output,-1)

    def get_star(self,star, force=0):        
        self._get_Part(star,"stars",force=force)
        
    def get_part(self,part, force=0):        
        self._get_Part(part,"part", force=force)
        
    def _get_Part(self,part, type,force=0):
        """
        get part in R200
        """
                
        name = self.path+type
        if os.path.isfile(name) and not force:
            with open(name, 'rb') as input:
                print("Reading %s"%name)
                setattr(self,type,pickle.load(input))
            return 
        
        x=part.x.data
        y=part.y.data
        z=part.z.data
        
        tree = spatial.cKDTree( np.transpose( [x,y,z] )) 

        n=self.nfoftot
        part=np.empty(n, dtype=np.object)

        for i in range(n):
            x=self.x[i]
            y=self.y[i]
            z=self.z[i]
            r=self.R200[i]

            part[i]=tree.query_ball_point((x,y,z), r)

            # Boundary conditions                
#             bound=[]
#             if x-r<0:
#                 bound.append(tree.query_ball_point((x+1,y,z), r))
#             if x+r>1:
#                 bound.append(tree.query_ball_point((x-1,y,z), r))
#             if y-r<0:
#                 bound.append(tree.query_ball_point((x,y+1,z), r))
#             if y+r>1:
#                 bound.append(tree.query_ball_point((x,y-1,z), r))
#             if z-r<0:
#                 bound.append(tree.query_ball_point((x,y,z+1), r))
#             if z+r>1:
#                 bound.append(tree.query_ball_point((x,y,z-1), r))

#                 if len(bound):
#                     print(bound)
#                     part[i].append(bound)

        with open(name, 'wb') as output:
            pickle.dump(part, output,-1)

        setattr(self,type,part)
    
    def get_star_fine(self,grid, stars, force=0):  
        """
        get halo stars the fine way
        """
        
        name = self.path+"stars_fine"
        if os.path.isfile(name) and not force:
            print("Reading %s"%name)
            with open(name, 'rb') as input:
                self.stars_fine = pickle.load(input)
            return

        x=stars.x.data
        y=stars.y.data
        z=stars.z.data
        tree = spatial.cKDTree( np.transpose( [x,y,z] ))

        self.stars_fine=np.empty(self.nfoftot, dtype=np.object)
        for i in range(self.nfoftot):
            cells = self.cells_fine[i]

            l=grid.l.data[cells]

            dx=np.power(0.5,l)

            x=grid.x.data[cells]+dx/2
            y=grid.y.data[cells]+dx/2
            z=grid.z.data[cells]+dx/2  

            r=dx/2 *np.sqrt(3)            

            search=[]
            for j in range(len(cells)):
                search.append( tree.query_ball_point( (x[j],y[j],z[j]), r[j] ) )

            self.stars_fine[i]=np.int32(np.unique(np.concatenate(search)))

        with open(name, 'wb') as output:
            pickle.dump(self.stars_fine, output,-1)

    def get_cells(self,grid, force=0):
        """
        get grid cells in R200
        """
        
        name = self.path+"cells"
        if os.path.isfile(name) and not force:
            print("Reading %s"%name)
            with open(name, 'rb') as input:
                self.cells = pickle.load(input)            
            return
        
        #get the center of cells
        l=grid.l.data
        dx=np.power(0.5,l+1)
        x=grid.x.data+dx
        y=grid.y.data+dx
        z=grid.z.data+dx

        tree = spatial.cKDTree( np.transpose( [x,y,z] ))
        self.cells=np.zeros(self.nfoftot,dtype=np.object)

        for halo_num in range(self.nfoftot):
            xc=self.x[halo_num]
            yc=self.y[halo_num]
            zc=self.z[halo_num]
            r=self.R200[halo_num]
            self.cells[halo_num]=tree.query_ball_point((xc,yc,zc), r)

            # Boundary conditions
#             if xc-r<0:
#                 self.cells[halo_num].append(tree.query_ball_point((xc+1,yc,zc), r))
#             if xc+r>1:
#                 self.cells[halo_num].append(tree.query_ball_point((xc-1,yc,zc), r))
#             if yc-r<0:
#                 self.cells[halo_num].append(tree.query_ball_point((xc,yc+1,zc), r))
#             if yc+r>1:
#                 self.cells[halo_num].append(tree.query_ball_point((xc,yc-1,zc), r))
#             if zc-r<0:
#                 self.cells[halo_num].append(tree.query_ball_point((xc,yc,zc+1), r))
#             if zc+r>1:
#                 self.cells[halo_num].append(tree.query_ball_point((xc,yc,zc-1), r))

        with open(name, 'wb') as output:
            pickle.dump(self.cells, output,-1)

    def get_cells_fine(self,grid, force=0):
        """
        get grid cells of halo by associating each part to the nearest cell
        """
        
        name = self.path+"cells_fine"
        if os.path.isfile(name) and not force:
            print("Reading %s"%name)
            with open(name, 'rb') as input:
                self.cells_fine = pickle.load(input)            
            return

        #get the center of cells
        l=grid.l.data
        dx=np.power(0.5,l+1)
        x=grid.x.data+dx
        y=grid.y.data+dx
        z=grid.z.data+dx

        tree = spatial.cKDTree( np.transpose( [x,y,z] ))

        self.cells_fine=np.empty(self.nfoftot,dtype=np.object)
        for halo_num in range(self.nfoftot):

            xp=self.part_pos[halo_num][0::3]
            yp=self.part_pos[halo_num][1::3]
            zp=self.part_pos[halo_num][2::3]

            self.cells_fine[halo_num]=np.unique(tree.query(np.transpose( [xp,yp,zp] ))[1])

        with open(name, 'wb') as output:
            pickle.dump(self.cells_fine, output,-1)

####################################################################################################################
# Mass 
####################################################################################################################

    def get_star_mass(self,star):
        self._get_Part_mass(star, type="star_mass")
        
    def get_part_mass(self,part):        
        self._get_Part_mass(part, type="part_mass")        
        
    def _get_Part_mass(self,part,type):
        """
        get part mass in Mo
        """    
        from pyemma import io
        info = io.Info(self.path+"../../../")
        
        if type == "star_mass":
            partID = self.stars
        if type == "part_mass":
            partID = self.part
            
        part_mass=np.zeros(self.nfoftot)
        for i in range(self.nfoftot):
            part_mass[i]=np.sum(part.mass.data[partID[i]])
        part_mass = part_mass/1.9891e30*info.unit_mass 
        setattr(self,type,part_mass)

    def get_part_mass_fine(self): 
        """
        get halo mass in Mo
        """

        nptot=2**(3*self.lmin)
        part_mass=(1.-self.ob/self.om)/(nptot)
                        
        self.part_mass_fine=self.npart*part_mass
        self.part_mass_fine=self.part_mass_fine/1.9891e30*self.unit_mass
                    
    def get_star_mass_fine(self,star):
        """
        get halo star mass in Mo
        """    
        from pyemma import io
        info = io.Info(self.path+"../../../")
        
        self.star_mass_fine=np.zeros(self.nfoftot)
        for i in range(self.nfoftot):
            stars=self.stars_fine[i]            
            self.star_mass_fine[i]=np.sum(star.mass.data[stars])
        self.star_mass_fine = self.star_mass_fine/1.9891e30*info.unit_mass
        
    def get_gas_mass(self,grid):
        """
        get halo gas mass in Mo in R200
        """

        from pyemma import io
        info = io.Info(self.path+"../../../")
        
        self.gas_mass=np.zeros(self.nfoftot)
        for halo_num in range(self.nfoftot):
            cells=self.cells[halo_num]
            l=grid.l.data[cells]
            d=grid.field_d.data[cells]
            self.gas_mass[halo_num] = np.sum( d* np.power(0.5,3*l) )

        self.gas_mass = self.gas_mass/1.9891e30*info.unit_mass

    def get_gas_mass_fine(self,grid, info):
        """
        get halo gas mass in Mo in fine methode
        """        
                
        self.gas_mass_fine=np.zeros(self.nfoftot)
        for halo_num in range(self.nfoftot):
            cells=self.cells_fine[halo_num]
            l=grid.l.data[cells]
            d=grid.field_d.data[cells]
            self.gas_mass_fine[halo_num] = np.sum( d* np.power(0.5,3*l) )

        self.gas_mass_fine = self.gas_mass_fine/1.9891e30*info.unit_mass                
        
####################################################################################################################
# Luminosity
####################################################################################################################
        
    def get_integ_egy(self,age,mass,tlife,force=0):
        import luminosity
        name = self.path+"int_egy"
        if os.path.isfile(name) and not force:
            with open(name, 'rb') as input:
                self.integ_egy = pickle.load(input)
            return
            
        self.integ_egy = np.zeros(self.nfoftot)
        for halo_num in range(self.nfoftot):
            stars = self.stars[halo_num]
            for i in self.stars[halo_num]:                    
                self.integ_egy[halo_num]+=mass[i] * luminosity.get_tot_egy(age[i],tlife)

        with open(name, 'wb') as output:
            pickle.dump(self.integ_egy, output,-1)
                            
    def get_luminosity_1600(self,cur_step):
        from pyemma import io,time,luminosity
        
        info = io.Info(self.path+"../../../")

        t=time.a2t_quad(cur_step.a, info.om, info.H0)
        
        luminosity.get_all_flux_1600(cur_step.star,t,info.unit_mass)

        flux = cur_step.star.flux_1600
        
        flux_tot=np.zeros(self.nfoftot)
        for i in range(self.nfoftot):
            flux_tot[i]=np.sum(flux[self.stars[i]])
        self.star_flux_1600=flux_tot
        self.mag_1600=luminosity.flux2mag(flux_tot[flux_tot!=0])
        
    def get_luminosity_UV(self,cur_step):
        from pyemma import io, time,luminosity
        info = io.Info(self.path+"../../../")
                    
        t=time.a2t_quad(cur_step.a, info.om, info.H0)
        
        luminosity.get_all_flux_UV(cur_step.star,t,info.unit_mass)

        flux = cur_step.star.flux_UV
        
        flux_tot=np.zeros(self.nfoftot)
        for i in range(self.nfoftot):
            flux_tot[i]=np.sum(flux[self.stars[i]])
        self.star_flux_UV=flux_tot
        self.mag_UV=luminosity.flux2mag(flux_tot[flux_tot!=0])

####################################################################################################################
####################################################################################################################

    def get_mean_vel(self, force=False):
        """
        get the mean velocity of halos by averaging the speed of their particles
        """
        
        name = self.path+"mean_vel"
        if os.path.isfile(name) and not force:
            print("Reading %s"%name)
            with open(name, 'rb') as input:
                self.mean_vel = pickle.load(input)
            return
            
        self.mean_vel = np.empty(self.nfoftot, dtype=np.object)
        for halo_num in range(self.nfoftot):            
            vx=np.mean(self.part_vel[halo_num][0::3])
            vy=np.mean(self.part_vel[halo_num][1::3])
            vz=np.mean(self.part_vel[halo_num][2::3])

            self.mean_vel[halo_num] = [vx,vy,vz]

        with open(name, 'wb') as output:
            pickle.dump(self.mean_vel, output,-1)

####################################################################################################################
# GEOMETRY
####################################################################################################################

    def get_inertia_tensor(self, force=False):
        """
        compute inertia tensors of halos and diagonalize them
        """
                
        name = self.path+"inertia"
        if os.path.isfile(name) and not force:
            print("Reading %s"%name)
            with open(name, 'rb') as input:
                self.inertia_eig_val = pickle.load(input)
                self.inertia_eig_vec = pickle.load(input)
            return
        
        import scipy.linalg

        self.inertia_eig_val = np.empty(self.nfoftot, dtype=np.object)
        self.inertia_eig_vec = np.empty(self.nfoftot, dtype=np.object)

        for halo_num in range(self.nfoftot):

            xc =self.x[halo_num]
            yc =self.y[halo_num]
            zc =self.z[halo_num]

            part_x = self.part_pos[halo_num][0::3] - xc
            part_y = self.part_pos[halo_num][1::3] - yc
            part_z = self.part_pos[halo_num][2::3] - zc

            #inertia matrix
            A = [part_x, part_y, part_z]
            B = np.transpose(A)
            I = np.dot(A,B)

            self.inertia_eig_val[halo_num], self.inertia_eig_vec[halo_num] = scipy.linalg.eig(I)
        
        with open(name, 'wb') as output:
            pickle.dump(self.inertia_eig_val, output,-1)
            pickle.dump(self.inertia_eig_vec, output,-1)

    def get_mass_center(self, cur_step, force=False):

#         cur_step.fof.get_cells_fine(cur_step.grid)
#         cur_step.fof.get_stars_fine(cur_step.grid,cur_step.star)

        self.mass_center_dm_x = np.zeros(self.nfoftot)
        self.mass_center_dm_y = np.zeros(self.nfoftot)
        self.mass_center_dm_z = np.zeros(self.nfoftot)

        self.mass_center_star_x= np.zeros(self.nfoftot)
        self.mass_center_star_y= np.zeros(self.nfoftot)
        self.mass_center_star_z= np.zeros(self.nfoftot)

        self.mass_center_gas_x= np.zeros(self.nfoftot)
        self.mass_center_gas_y= np.zeros(self.nfoftot)
        self.mass_center_gas_z= np.zeros(self.nfoftot)

        self.mass_center_baryon_x= np.zeros(self.nfoftot)
        self.mass_center_baryon_y= np.zeros(self.nfoftot)
        self.mass_center_baryon_z= np.zeros(self.nfoftot)        

        for halo_num in range(cur_step.fof.nfoftot):

            self.mass_center_dm_x[halo_num] = np.mean(self.part_pos[halo_num][0::3])
            self.mass_center_dm_y[halo_num] = np.mean(self.part_pos[halo_num][1::3])
            self.mass_center_dm_z[halo_num] = np.mean(self.part_pos[halo_num][2::3])

            gas_m_tot = 0
            _gas_x = 0
            _gas_y = 0
            _gas_z = 0
            cells_mask = self.cells_fine[halo_num]
            if len(cells_mask):
                dv = np.power(0.5,3*cur_step.grid.l.data[cells_mask])
                _gas_m = cur_step.grid.field_d.data[cells_mask] * dv
                _gas_x = np.sum(cur_step.grid.x.data[cells_mask]*_gas_m)
                _gas_y = np.sum(cur_step.grid.y.data[cells_mask]*_gas_m)
                _gas_z = np.sum(cur_step.grid.z.data[cells_mask]*_gas_m)

                gas_m_tot = np.sum(_gas_m)


            star_m_tot = 0
            _star_x = 0
            _star_y = 0
            _star_z = 0    
            star_mask = self.stars_fine[halo_num]
            if len(star_mask):
                _star_m = cur_step.star.mass.data[star_mask]
                _star_x = np.sum(cur_step.star.x.data[star_mask]*_star_m)
                _star_y = np.sum(cur_step.star.y.data[star_mask]*_star_m)
                _star_z = np.sum(cur_step.star.z.data[star_mask]*_star_m)

                star_m_tot = np.sum(_star_m)
            

            if (star_m_tot>0):    
                self.mass_center_star_x[halo_num]= _star_x / star_m_tot
                self.mass_center_star_y[halo_num]= _star_y / star_m_tot
                self.mass_center_star_z[halo_num]= _star_z / star_m_tot

            if (gas_m_tot>0):
                self.mass_center_gas_x[halo_num]= _gas_x / gas_m_tot
                self.mass_center_gas_y[halo_num]= _gas_y / gas_m_tot
                self.mass_center_gas_z[halo_num]= _gas_z / gas_m_tot                

            m_tot = star_m_tot + gas_m_tot
            if (m_tot>0):
                self.mass_center_baryon_x[halo_num]= (_gas_x + _star_x) / m_tot
                self.mass_center_baryon_y[halo_num]= (_gas_y + _star_y) / m_tot
                self.mass_center_baryon_z[halo_num]= (_gas_z + _star_z) / m_tot  
            
####################################################################################################################
# FLUX AT R200
####################################################################################################################

    def get_flux_r200(self,grid, type, fact=1, force=0):
        """
        type could be : hydro, rad, ref
        
        #read mean vel
        #get_cells
        """

        name = self.path+"flux_r200_%s"%type
        fact_name = ("%0.2f"%(fact)).replace(".","_")  
        if os.path.isfile(name) and not force:
            print("Reading %s"%name)
            with open(name, 'rb') as input:                      
                setattr(self,"flux_data_%s_%s"%(type,fact_name),pickle.load(input))
                setattr(self,"mean_flux_%s_%s"%(type,fact_name),pickle.load(input))
            return                
        
        #healpix shere
        import healpy as hp        
        n=4
        nside = 2**n
        x_healpix,y_healpix,z_healpix=hp.pix2vec(nside, range(hp.nside2npix(nside) ))
            
        mean_flux = np.zeros(self.nfoftot)        
        flux_data=np.empty(self.nfoftot, dtype=np.object)        

        skipped = 0
        for halo_num in range(self.nfoftot):
            if not halo_num%(1000):
                print(halo_num)

            #select halo in catalog
            xc =self.x[halo_num]
            yc =self.y[halo_num]
            zc =self.z[halo_num]
            R200=self.R200[halo_num]
            cells = self.cells[halo_num]

            # skip halo smaller than ncells_min
            ncells_min = 2
            if (len(cells)<ncells_min):
                skipped +=1
                continue

            # get cell size
            l_grid = grid.l.data[cells]
            dx = np.power(0.5,l_grid)

            # get halo cells, set cell-centred coordinate, center on the halo
            x_grid = (grid.x.data[cells] +dx/2. -xc)
            y_grid = (grid.y.data[cells] +dx/2. -yc)
            z_grid = (grid.z.data[cells] +dx/2. -zc)

            #gen tree
            tree = spatial.cKDTree(np.transpose((x_grid,y_grid,z_grid)))
            
            #get the cell nearest neighbourg of each healpix point
            idx = tree.query(np.transpose((x_healpix*fact*R200,y_healpix*fact*R200,z_healpix*fact*R200)))[1]
        
            # Select the flux
            if type == "hydro":
                # hydro flux
                fx=grid.field_d.data[cells]*(grid.field_u.data[cells]-self.mean_vel[halo_num][0])
                fy=grid.field_d.data[cells]*(grid.field_v.data[cells]-self.mean_vel[halo_num][1])
                fz=grid.field_d.data[cells]*(grid.field_w.data[cells]-self.mean_vel[halo_num][2])
            
            if type == "rad":
                # radiative flux
                fx=grid.rfield_fx0.data[cells]
                fy=grid.rfield_fy0.data[cells]
                fz=grid.rfield_fz0.data[cells]

            if type == "ref":
                #reference flux
                fx=-np.ones(len(cells))
                fy=-np.ones(len(cells))
                fz=-np.ones(len(cells))
                
            # scalar product
            flux_data[halo_num]=(x_healpix*fx[idx] +
                                 y_healpix*fy[idx] +
                                 z_healpix*fz[idx] )    *4*np.pi*(fact*R200)**2
            # mean flux
            mean_flux[halo_num]=np.mean(flux_data[halo_num])

        # mean_scal_rad[abs(mean_scal_rad) < abs(mean_scal_ref) ]=0
        
        print("skipped %d/%d = %.02f %%"%(skipped,self.nfoftot, skipped/self.nfoftot*100))
                                
        setattr(self,"flux_data_%s_%s"%(type,fact_name),flux_data)
        setattr(self,"mean_flux_%s_%s"%(type,fact_name),mean_flux)
        
        with open(name, 'wb') as output:
            pickle.dump(flux_data, output,-1)
            pickle.dump(mean_flux, output,-1)
        
        print("get_flux_r200 done")
        
####################################################################################################################
####################################################################################################################
    def get_instant_SFR(self,stars):
        from pyemma import io, time
        info = io.Info(self.path+"../../../")
        
        t=time.a2t_quad(stars.age.tsim, info.om, info.H0)
        
        self.instant_sfr=np.zeros(self.nfoftot)
        for halo_num in range(self.nfoftot):
            star_mask=self.stars[halo_num]
            age=t-stars.age.data[star_mask]
            thresh=1e7
            self.instant_sfr[halo_num]=np.sum( (stars.mass.data[star_mask])[age<thresh] )
            self.instant_sfr[halo_num]*=info.unit_mass/1.9891e30
            self.instant_sfr[halo_num]/=thresh
            
    def get_time_newest_star(self, stars):
        self.time_newest_star=np.zeros(self.nfoftot)        
        for halo_num in range(self.nfoftot):
            star_mask=self.stars[halo_num]
            if len(star_mask):
                self.time_newest_star[halo_num] = np.max(stars.age.data[star_mask])                
                
    def get_time_oldest_star(self, stars):
        self.time_oldest_star=np.zeros(self.nfoftot)
        for halo_num in range(self.nfoftot):
            star_mask=self.stars[halo_num]
            if len(star_mask):
                self.time_oldest_star[halo_num] = np.min(stars.age.data[star_mask])

