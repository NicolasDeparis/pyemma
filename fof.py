# coding: utf-8

import numpy as np
import pickle
# import cPickle as pickle
import os
from scipy import spatial

class Fof:
    def __init__(self,folder,stepnum,step):
        self.exec_folder = "/home/deparis/Emma/utils/pfof/"

        self.folder=folder
        self.stepnum=stepnum
        self.path="%s%05d/halo/"%(self.folder,self.stepnum)

        self._isLoaded_masst = False
        self._isLoaded_struct = False
        self.step=step

    def __getattr__(self, name):

        if name == 'R200':
            print("Getting %s"%name)
            self.get_R200()
            return self.__getattribute__(name)

        elif name == 'lmin':
            print("Getting %s"%name)
            self.getLmin()
            return self.__getattribute__(name)

        elif name == 'ob':
            print("Getting %s"%name)
            self.getCosmo()
            return self.__getattribute__(name)

        elif name == 'om':
            print("Getting %s"%name)
            self.getCosmo()
            return self.__getattribute__(name)

        elif name == 'nfoftot':
            print("Getting %s"%name)
            self.getNfofTot()
            return self.__getattribute__(name)

        elif name == 'x':
            print("Getting %s"%name)
            self.read_masst()
            return self.__getattribute__(name)

        elif name == 'y':
            print("Getting %s"%name)
            self.read_masst()
            return self.__getattribute__(name)

        elif name == 'z':
            print("Getting %s"%name)
            self.read_masst()
            return self.__getattribute__(name)

        elif name == 'idx':
            print("Getting %s"%name)
            self.read_masst()
            return self.__getattribute__(name)

        elif name == 'npart':
            print("Getting %s"%name)
            self.read_masst()
            return self.__getattribute__(name)

        elif name == 'part_n':
            print("Getting %s"%name)
            self.read_masst()
            return self.__getattribute__(name)

        elif name == 'part_pos':
            print("Getting %s"%name)
            self.read_struct()
            return self.__getattribute__(name)

        elif name == 'part_vel':
            print("Getting %s"%name)
            self.read_struct()
            return self.__getattribute__(name)

        elif name == 'part_id':
            print("Getting %s"%name)
            self.read_struct()
            return self.__getattribute__(name)

        elif name == 'nproc_sim':
            print("Getting %s"%name)
            self.getNproc_sim()
            return self.__getattribute__(name)

        elif name == 'ncpu':
            print("Getting %s"%name)
            self.getNcpu()
            return self.__getattribute__(name)

        elif name == 'inertia_eig_val':
            print("Getting %s"%name)
            self.get_inertia_tensor()
            return self.__getattribute__(name)

        elif name == 'inertia_eig_vec':
            print("Getting %s"%name)
            self.get_inertia_tensor()
            return self.__getattribute__(name)

        elif name == 'part_mass':
            print("Getting %s"%name)
            self.get_part_mass()
            return self.__getattribute__(name)

        elif name == 'part_mass_fine':
            print("Getting %s"%name)
            self.get_part_mass_fine()
            return self.__getattribute__(name)

        elif name == 'getStars':
            print("Getting %s"%name)
            self.get_getStars()
            return self.__getattribute__(name)

        elif name == 'stars':
            print("Getting %s"%name)
            self.get_stars()
            return self.__getattribute__(name)

        elif name == 'stars_fine':
            print("Getting %s"%name)
            self.get_stars_fine()
            return self.__getattribute__(name)

        elif name == 'part':
            print("Getting %s"%name)
            self.get_part()
            return self.__getattribute__(name)

        elif name == 'cells':
            print("Getting %s"%name)
            self.get_cells()
            return self.__getattribute__(name)

        elif name == 'cells_fine':
            print("Getting %s"%name)
            self.get_cells_fine()
            return self.__getattribute__(name)

        elif name == 'gas_mass':
            print("Getting %s"%name)
            self.get_gas_mass()
            return self.__getattribute__(name)

        elif name == 'gas_mass_fine':
            print("Getting %s"%name)
            self.get_gas_mass_fine()
            return self.__getattribute__(name)

        elif name == 'star_mass':
            print("Getting %s"%name)
            self.get_star_mass()
            return self.__getattribute__(name)

        elif name == 'star_mass_fine':
            print("Getting %s"%name)
            self.get_star_mass_fine()
            return self.__getattribute__(name)

        elif name == 'instant_SFR':
            print("Getting %s"%name)
            self.get_instant_SFR()
            return self.__getattribute__(name)

        elif name == 'instant_SFR_from_gas':
            print("Getting %s"%name)
            self.get_instant_SFR_from_gas()
            return self.__getattribute__(name)

        else:
            print("ERROR : Can't automatically determine %s"%name)
            raise AttributeError
            return

    def write_infosim(self):
        self.getLmin()
        with open("%sinfosim.txt"%self.exec_folder,"w") as f:
            # f.write("Npro %d\n"%self.nproc_sim)
            f.write("Npro %d\n"%1)
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

    def gen(self,nproc_fof, force=0):

        folder_name ="%s%05d/halo/"%(self.folder, self.stepnum)

        if not force:
            try:
                files=os.listdir(folder_name)
                for f in files:
                    if "halo_masst" in f:
                        print("fof already exist")
                        return 0
            except OSError:
                pass

        self.ncpu = nproc_fof
        self.write_fofin()
        self.write_infosim()

        self.convert()

        curdir=os.getcwd()
        os.chdir(self.exec_folder)

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

            files=os.listdir(self.path)
            self.ncpu=0
            for f in files:
                if "halo_masst" in f:
                    self.ncpu+=1
        except FileNotFoundError:
            print("ERROR : can't determine nproc in \n%s"%(self.path))

    def getNproc_sim(self):
        """
        return the number of cpu used by the simulation by reading param.run
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


    def convert(self):

        import h5py

        path = self.folder
        stepnum = self.stepnum



        folder ="%s%05d/"%(path, stepnum)
        f = h5py.File("%spart_x_%05d.h5"%(folder,stepnum), "r")
        tsim = f.attrs['a']
        data = f['data'][:]
        f.close()

        try:
            os.mkdir("%s/part_x"%folder)
        except FileExistsError:
            pass
        f=open("%s/part_x/x.%05d.p00000"%(folder,stepnum), "wb")
        f.write(np.array([len(data)], dtype=np.int32).tobytes())
        f.write(np.array([1], dtype=np.int32).tobytes())
        f.write(np.array([tsim], dtype=np.float32).tobytes())
        for i in range(7):
            f.write(np.array([1], dtype=np.int32).tobytes())
        f.write(np.array(data, dtype=np.float32).tobytes())
        f.close()

        f = h5py.File("%spart_y_%05d.h5"%(folder,stepnum), "r")
        tsim = f.attrs['a']
        data = f['data'][:]
        f.close()

        try:
            os.mkdir("%s/part_y"%folder)
        except FileExistsError:
            pass
        f=open("%spart_y/y.%05d.p00000"%(folder,stepnum), "wb")
        f.write(np.array([len(data)], dtype=np.int32).tobytes())
        f.write(np.array([1], dtype=np.int32).tobytes())
        f.write(np.array([tsim], dtype=np.float32).tobytes())
        for i in range(7):
            f.write(np.array([1], dtype=np.int32).tobytes())
        f.write(np.array(data, dtype=np.float32).tobytes())
        f.close()

        try:
            os.mkdir("%s/part_z"%folder)
        except FileExistsError:
            pass
        f = h5py.File("%spart_z_%05d.h5"%(folder,stepnum), "r")
        tsim = f.attrs['a']
        data = f['data'][:]
        f.close()
        f=open("%spart_z/z.%05d.p00000"%(folder,stepnum), "wb")
        f.write(np.array([len(data)], dtype=np.int32).tobytes())
        f.write(np.array([1], dtype=np.int32).tobytes())
        f.write(np.array([tsim], dtype=np.float32).tobytes())
        for i in range(7):
            f.write(np.array([1], dtype=np.int32).tobytes())
        f.write(np.array(data, dtype=np.float32).tobytes())
        f.close()

        f = h5py.File("%spart_vx_%05d.h5"%(folder,stepnum), "r")
        tsim = f.attrs['a']
        data = f['data'][:]
        f.close()

        try:
            os.mkdir("%s/part_vx"%folder)
        except FileExistsError:
            pass
        f=open("%spart_vx/vx.%05d.p00000"%(folder,stepnum), "wb")
        f.write(np.array([len(data)], dtype=np.int32).tobytes())
        f.write(np.array([1], dtype=np.int32).tobytes())
        f.write(np.array([tsim], dtype=np.float32).tobytes())
        for i in range(7):
            f.write(np.array([1], dtype=np.int32).tobytes())
        f.write(np.array(data, dtype=np.float32).tobytes())
        f.close()

        f = h5py.File("%spart_vy_%05d.h5"%(folder,stepnum), "r")
        tsim = f.attrs['a']
        data = f['data'][:]
        f.close()

        try:
            os.mkdir("%s/part_vy"%folder)
        except FileExistsError:
            pass
        f=open("%spart_vy/vy.%05d.p00000"%(folder,stepnum), "wb")
        f.write(np.array([len(data)], dtype=np.int32).tobytes())
        f.write(np.array([1], dtype=np.int32).tobytes())
        f.write(np.array([tsim], dtype=np.float32).tobytes())
        for i in range(7):
            f.write(np.array([1], dtype=np.int32).tobytes())
        f.write(np.array(data, dtype=np.float32).tobytes())
        f.close()

        f = h5py.File("%spart_vz_%05d.h5"%(folder,stepnum), "r")
        tsim = f.attrs['a']
        data = f['data'][:]
        f.close()

        try:
            os.mkdir("%s/part_vz"%folder)
        except FileExistsError:
            pass
        f=open("%spart_vz/vz.%05d.p00000"%(folder,stepnum), "wb")
        f.write(np.array([len(data)], dtype=np.int32).tobytes())
        f.write(np.array([1], dtype=np.int32).tobytes())
        f.write(np.array([tsim], dtype=np.float32).tobytes())
        for i in range(7):
            f.write(np.array([1], dtype=np.int32).tobytes())
        f.write(np.array(data, dtype=np.float32).tobytes())
        f.close()


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

    def get_stars(self, force=0, fact=1., Rmin=None):
        self._get_Sphere(self.step.star,"stars",force=force, fact=fact)

    def get_part(self, force=0, fact=1., Rmin=None):
        self._get_Sphere(self.step.part,"part", force=force, fact=fact)

    def get_cells(self, force=0, fact=1., Rmin=None):
        self._get_Sphere(self.step.grid,"cells", force=force, fact=fact)

    def _get_Sphere(self,part, type,force=0, fact=1., Rmin=None):
        """
        get part in R200
        """

        name = self.path+type
        if os.path.isfile(name) and not force:
            with open(name, 'rb') as input:
                print("Reading %s"%name)
                setattr(self,type,pickle.load(input))
            return

        try :
            x=part.x.data
            y=part.y.data
            z=part.z.data
        except AttributeError:
            return

        if type == "cells":
            #get the center of cells
            dx=np.power(0.5,part.l.data+1)
            x+=dx
            y+=dx
            z+=dx

        tree = spatial.cKDTree( np.transpose( [x,y,z] ))

        n=self.nfoftot
        part=np.empty(n, dtype=np.object)

        N_halo_smaller_than_Rmin=0
        for i_halo in range(n):
            x=self.x[i_halo]
            y=self.y[i_halo]
            z=self.z[i_halo]
            r=self.R200[i_halo]*fact

            if Rmin is not None and r<Rmin:
                N_halo_smaller_than_Rmin+=1
                r=Rmin

            # Boundary conditions
            xm= (x-r)<0
            xp= (x+r)>=1
            xo= not(xm) and not(xp)

            ym= (y-r)<0
            yp= (y+r)>=1
            yo= not(ym) and not(yp)

            zm= (z-r)<0
            zp= (z-r)>=1
            zo= not(zm) and not(zp)

            if xo and yo and zo:
                part[i_halo]=np.uint(tree.query_ball_point((x,y,z), r))
                continue
            else:
                c=[1,0,-1]
                for k in range(3):
                    for j in range(3):
                        for i in range(3):
                            p=(x+c[i],y+c[j],z+c[k])
                            part_tmp=np.uint(tree.query_ball_point(p,r))

                            if part[i_halo] is None:
                                part[i_halo]=part_tmp
                            else:
                                part[i_halo]=np.append(part[i_halo], part_tmp)

        if Rmin is not None:
            print("Halos smaller than Rmin : ",N_halo_smaller_than_Rmin, Rmin)

        with open(name, 'wb') as output:
            pickle.dump(part, output,-1)

        setattr(self,type,part)

    def get_stars_radius(self):
        """
        get distance of all stars in halo compare to the center of the halo
        """

        stars=self.step.stars
        self.stars_radius=np.zeros(self.nfoftot,dtype=np.object)

        for i in range(self.nfoftot):
            if self.getStars[i]:

                dx=self.x[i]- stars.x.data[self.stars[i]]
                dy=self.y[i]- stars.y.data[self.stars[i]]
                dz=self.z[i]- stars.z.data[self.stars[i]]

                self.stars_radius[i]=np.sqrt(np.power(dx,2)+np.power(dy,2)+np.power(dz,2))

        # def get_part_fine(self,force=0):
        #     """
        #     get halo stars the fine way
        #     """
        #     grid=self.step.grid
        #     stars=self.step.part
        #     name = self.path+"part_fine"
        #     if os.path.isfile(name) and not force:
        #         print("Reading %s"%name)
        #         with open(name, 'rb') as input:
        #             self.stars_fine = pickle.load(input)
        #         return
        #
        #     x=stars.x.data
        #     y=stars.y.data
        #     z=stars.z.data
        #     tree = spatial.cKDTree( np.transpose( [x,y,z] ))
        #
        # self.stars_fine=np.empty(self.nfoftot, dtype=np.object)
        # for i in range(self.nfoftot):
        #     cells = self.cells_fine[i]
        #
        #     l=grid.l.data[cells]
        #
        #     dx=np.power(0.5,l)
        #
        #     x=grid.x.data[cells]+dx/2
        #     y=grid.y.data[cells]+dx/2
        #     z=grid.z.data[cells]+dx/2
        #
        #     r=dx/2 *np.sqrt(3)
        #
        #     search=[]
        #     for j in range(len(cells)):
        #         search.append( tree.query_ball_point( (x[j],y[j],z[j]), r[j] ) )
        #
        #     self.stars_fine[i]=np.int32(np.unique(np.concatenate(search)))
        #
        # with open(name, 'wb') as output:
        #     pickle.dump(self.stars_fine, output,-1)


    def get_stars_fine(self,force=0):
        """
        get halo stars the fine way
        """
        grid=self.step.grid
        stars=self.step.star
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

    def get_getStars(self):
        """
        get halo with stars
        """

        self.getStars=np.zeros(self.nfoftot, dtype=np.bool)
        for i in range(self.nfoftot):
            if len(self.stars[i])>0:
                self.getStars[i]=True

    def get_getYoungStars(self,age_max):
        """
        Check if halo get stars of age lower than age_max
        """
        star=self.step.star

        n=self.nfoftot
        self.getYoungStars=np.zeros(n, dtype=np.bool)
        for i in range(self.nfoftot):
            if self.getStars[i]:

                # t=np.max(star.age.data)
                t=self.step.t
                youngest_stars = np.min(t- star.age.data[self.stars[i]])
                if youngest_stars <= age_max:
                    self.getYoungStars[i] = True

    # def get_cells(self,force=0,fact=1.):
    #     """
    #     get grid cells in R200
    #     """
    #
    #     grid=self.step.grid
    #     name = self.path+"cells"
    #     if os.path.isfile(name) and not force:
    #         print("Reading %s"%name)
    #         with open(name, 'rb') as input:
    #             self.cells = pickle.load(input)
    #         return
    #
    #     #get the center of cells
    #     l=grid.l.data
    #     dx=np.power(0.5,l+1)
    #     x=grid.x.data+dx
    #     y=grid.y.data+dx
    #     z=grid.z.data+dx
    #
    #     tree = spatial.cKDTree( np.transpose( [x,y,z] ))
    #     self.cells=np.zeros(self.nfoftot,dtype=np.object)
    #
    #     for halo_num in range(self.nfoftot):
    #         xc=self.x[halo_num]
    #         yc=self.y[halo_num]
    #         zc=self.z[halo_num]
    #         r=self.R200[halo_num]
    #         self.cells[halo_num]=tree.query_ball_point((xc,yc,zc), fact*r)
    #
    #         # Boundary conditions
    #         xm= (x-r)<0
    #         xp= (x+r)>=1
    #         xo= not(xm) and not(xp)
    #
    #         ym= (y-r)<0
    #         yp= (y+r)>=1
    #         yo= not(ym) and not(yp)
    #
    #         zm= (z-r)<0
    #         zp= (z-r)>=1
    #         zo= not(zm) and not(zp)
    #
    #         if xo and yo and zo:
    #             self.cells[halo_num]=np.uint(tree.query_ball_point((x,y,z), r))
    #             continue
    #         else:
    #             c=[1,0,-1]
    #             for k in range(3):
    #                 for j in range(3):
    #                     for i in range(3):
    #                         p=(x+c[i],y+c[j],z+c[k])
    #                         part_tmp=np.uint(tree.query_ball_point(p,r))
    #
    #                         if self.cells[halo_num] is None:
    #                             self.cells[halo_num]=part_tmp
    #                         else:
    #                             self.cells[halo_num]=np.append(self.cells[halo_num], part_tmp)
    #
    #     with open(name, 'wb') as output:
    #         pickle.dump(self.cells, output,-1)

    def get_cells_fine(self,force=0):
        """
        get grid cells of halo by associating each part to the nearest cell
        """

        grid=self.step.grid

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

    def get_star_mass(self):
        self._get_Part_mass(self.step.star, type="star_mass")

    def get_part_mass(self):
        self._get_Part_mass(self.step.part, type="part_mass")

    def _get_Part_mass(self,part,type):
        """
        get part mass in Mo
        """
        from pyemma import io
        info = io.Info(self.path+"../../../")

        part_mass=np.zeros(self.nfoftot)

        if type == "star_mass":
            partID = self.stars
            for i in range(self.nfoftot):
                part_mass[i]=np.sum(part.mass.data[partID[i]])

        if type == "part_mass":
            partID = self.part
            unit_mass=2**(-3*self.lmin) * (1.-info.ob/info.om)
            for i in range(self.nfoftot):
                part_mass[i]=unit_mass*len(partID[i])

        part_mass = part_mass/1.9891e30*info.unit_mass
        setattr(self,type,part_mass)

    def get_part_mass_fine(self):
        """
        get halo mass in Mo
        """

        nptot=2**(3*self.lmin)
        part_mass=(1.-self.ob/self.om)/(nptot)

        Mo = self.unit_mass/1.9891e30
        self.part_mass_fine=self.npart*part_mass*Mo

    def get_star_mass_fine(self):
        """
        get halo star mass in Mo
        """

        star=self.step.star

        from pyemma import io
        info = io.Info(self.path+"../../../")

        self.star_mass_fine=np.zeros(self.nfoftot)
        for i in range(self.nfoftot):
            stars=self.stars_fine[i]
            self.star_mass_fine[i]=np.sum(star.mass.data[stars])
        self.star_mass_fine = self.star_mass_fine/1.9891e30*info.unit_mass

    def get_gas_mass(self):
        """
        get halo gas mass in Mo in R200
        """

        grid=self.step.grid
        from pyemma import io
        info = io.Info(self.path+"../../../")

        self.gas_mass=np.zeros(self.nfoftot)
        for halo_num in range(self.nfoftot):
            cells=self.cells[halo_num]
            l=grid.l.data[cells]
            d=grid.field_d.data[cells]
            self.gas_mass[halo_num] = np.sum( d* np.power(0.5,3*l) )

        self.gas_mass = self.gas_mass/1.9891e30*info.unit_mass

    def get_gas_mass_fine(self):
        """
        get halo gas mass in Mo in fine methode
        """

        grid=self.step.grid
        from pyemma import io
        info = io.Info(self.path+"../../../")

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

        fesc = 0.4
        print("WARNING applying escape fraction of %f"%fesc)
        flux *= fesc

        flux_tot=np.zeros(self.nfoftot)
        for i in range(self.nfoftot):
            flux_tot[i]=np.sum(flux[self.stars[i]])
        self.star_flux_1600=flux_tot
        self.mag_1600=luminosity.flux2mag(flux_tot[flux_tot!=0])

    def get_ageLuminosity(self,cur_step):
        """
        compute the age of stars when photon arriving at R200 were emited
        """
        t=np.max(cur_step.star.age.data)
        from pyemma import io,time,luminosity
        info = io.Info(self.path+"../../../")
        run = io.RunParam(self.path+"../../../")

        self.ageLuminosity=np.zeros(self.nfoftot,dtype=np.object)

        for i in range(self.nfoftot):
            if self.getStars[i]:
                # dt = (self.R200[i] - self.stars_radius[i]) *cur_step.a*info.unit_l /(run.clight*299792458) /(365*24*3600)
                # dt = self.R200[i] *cur_step.a *info.unit_l /(run.clight*299792458) /(365*24*3600) # yr
                dt=0
                self.ageLuminosity[i]= t -cur_step.star.age.data[self.stars[i]] -dt # yr

    def get_luminosity_UV(self,cur_step):
        from pyemma import io,luminosity
        run = io.RunParam(self.path+"../../../")
        info = io.Info(self.path+"../../../")

        tlife = run.tlife_rad #yr
        E0 = run.src_int_or_fesc*run.fesc #phot/s/kg

        def model_flux_UV(age,tlife,E0):
            y=np.ones(len(age)) *E0
            y[age>tlife] *= np.power(age[age>tlife]/tlife ,-4.)
            y[age<0] = 0
            return y

        self.star_flux_UV=np.zeros(self.nfoftot)
        for i in range(self.nfoftot):
            if self.getStars[i]:

                flux=model_flux_UV(self.ageLuminosity[i],tlife,E0) # #phot/s/kg
                mass=cur_step.star.mass.data[self.stars[i]] *info.unit_mass # kg
                self.star_flux_UV[i]=np.sum(flux*mass) #phot/s

        # self.mag_UV=luminosity.flux2mag(self.star_flux_UV)

    def get_luminosity_UV_old(self,cur_step, dt=0):
        from pyemma import io, time,luminosity
        info = io.Info(self.path+"../../../")
        run = io.RunParam(self.path+"../../../")

        t=time.a2t_quad(cur_step.a, info.om, info.H0)

        luminosity.get_all_flux_UV_old(cur_step.star,t,info.unit_mass,run)

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
            # vx=np.mean(self.part_vel[halo_num][0::3])
            # vy=np.mean(self.part_vel[halo_num][1::3])
            # vz=np.mean(self.part_vel[halo_num][2::3])

            # vx=np.mean(self.step.part.vx.data[self.part_id[halo_num]])
            # vy=np.mean(self.step.part.vy.data[self.part_id[halo_num]])
            # vz=np.mean(self.step.part.vz.data[self.part_id[halo_num]])

            vx=0
            vy=0
            vz=0

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

    def get_flux_r200(self,type,fact=1,force=0,getonly=0, Rmin=None):
        """
        type could be : hydro, rad, ref

        #read mean vel
        #get_cells
        """
        grid=self.step.grid
        name = self.path+"flux_r200_%s"%type
        fact_name = ("%0.2f"%(fact)).replace(".","_")
        if os.path.isfile(name) and not force:
            print("Reading %s"%name)
            with open(name, 'rb') as input:
                setattr(self,"flux_data_%s_%s"%(type,fact_name),pickle.load(input))
                setattr(self,"tot_flux_%s_%s"%(type,fact_name),pickle.load(input))
            return


        #healpix shere
        import healpy as hp
        n=4
        nside = 2**n
        n_healpix = hp.nside2npix(nside)
        x_healpix,y_healpix,z_healpix=hp.pix2vec(nside, range(n_healpix))

        tot_flux = np.zeros(self.nfoftot)
        flux_data=np.empty(self.nfoftot, dtype=np.object)

        # self.cos = np.zeros(self.nfoftot, dtype=np.object)
        # self.norm = np.zeros(self.nfoftot, dtype=np.object)
        # self.cosstd = np.zeros(self.nfoftot)
        # self.cosmean = np.zeros(self.nfoftot)


        #get the center of cells
        l=grid.l.data
        dx=np.power(0.5,l+1)
        x=grid.x.data+dx
        y=grid.y.data+dx
        z=grid.z.data+dx

        tree = spatial.cKDTree( np.transpose( [x,y,z] ))

        skipped = 0
        N_halo_smaller_than_Rmin=0
        for halo_num in range(self.nfoftot):
        # for halo_num in range(100):
            if not halo_num%(1000):
                print(halo_num)

            #select halo in catalog
            xc =self.x[halo_num]
            yc =self.y[halo_num]
            zc =self.z[halo_num]
            R200=self.R200[halo_num]

            if Rmin is not None:
                if R200<Rmin:
                    N_halo_smaller_than_Rmin+=1
                    R200=Rmin


            # cells = self.cells[halo_num]
            #
            # # skip halo smaller than ncells_min
            # ncells_min = 2
            # if (len(cells)<ncells_min):
            #     skipped +=1
            #     continue
            #
            # # l_min = 11
            # # if (np.min(grid.l.data[cells])<l_min):
            # #     skipped +=1
            # #     continue


            # skip boundary conditions
            # if (xc+R200 >=1.) | (xc-R200 <0.) | \
            #    (yc+R200 >=1.) | (yc-R200 <0.) | \
            #    (zc+R200 >=1.) | (zc-R200 <0.):
            #     skipped +=1
            #     continue
            #
            # # get cell size
            # l_grid = grid.l.data[cells]
            # dx = np.power(0.5,l_grid)
            #
            # # get halo cells, set cell-centred coordinate, center on the halo
            # x_grid = (grid.x.data[cells] +dx/2. -xc)
            # y_grid = (grid.y.data[cells] +dx/2. -yc)
            # z_grid = (grid.z.data[cells] +dx/2. -zc)

            #gen tree
            # tree = spatial.cKDTree(np.transpose((x_grid,y_grid,z_grid)))

            #get the cell nearest neighbourg of each healpix point
            xp=xc+x_healpix*fact*R200
            yp=yc+y_healpix*fact*R200
            zp=zc+z_healpix*fact*R200
            idx = tree.query(np.transpose((xp,yp,zp)))[1]


            # Select the flux
            if type == "hydro":
                # hydro flux
                # fx=grid.field_d.data[cells]*(grid.field_u.data[cells]-self.mean_vel[halo_num][0])
                # fy=grid.field_d.data[cells]*(grid.field_v.data[cells]-self.mean_vel[halo_num][1])
                # fz=grid.field_d.data[cells]*(grid.field_w.data[cells]-self.mean_vel[halo_num][2])

                fx=grid.field_d.data[idx]*(grid.field_u.data[idx]-self.mean_vel[halo_num][0])
                fy=grid.field_d.data[idx]*(grid.field_v.data[idx]-self.mean_vel[halo_num][1])
                fz=grid.field_d.data[idx]*(grid.field_w.data[idx]-self.mean_vel[halo_num][2])

            if type == "rho":
                r=grid.field_d.data[idx]

            if type == "speed":
                # hydro flux
                # fx=grid.field_u.data[cells]-self.mean_vel[halo_num][0]
                # fy=grid.field_v.data[cells]-self.mean_vel[halo_num][1]
                # fz=grid.field_w.data[cells]-self.mean_vel[halo_num][2]

                fx=grid.field_u.data[idx]-self.mean_vel[halo_num][0]
                fy=grid.field_v.data[idx]-self.mean_vel[halo_num][1]
                fz=grid.field_w.data[idx]-self.mean_vel[halo_num][2]


            if type == "rad":
                # radiative flux
                # fx=grid.rfield_fx0.data[cells]
                # fy=grid.rfield_fy0.data[cells]
                # fz=grid.rfield_fz0.data[cells]
                fx=grid.rfield_fx0.data[idx]
                fy=grid.rfield_fy0.data[idx]
                fz=grid.rfield_fz0.data[idx]


            # if type == "ref":
            #     #reference flux
            #     fx=np.ones(len(cells))*3
            #     fy=np.zeros(len(cells))
            #     fz=np.zeros(len(cells))

            #surface element
            ds= 4*np.pi*(fact*R200)**2 / n_healpix

            if type == "speed":
                ds= 1. / n_healpix
            if type == "rho":
                ds=1.

            # scalar product
            # flux_data[halo_num]=(x_healpix*fx[idx] +
            #                      y_healpix*fy[idx] +
            #                      z_healpix*fz[idx] ) * ds


            if type == "rho":
                flux_data[halo_num]= r
            else:
                flux_data[halo_num]=(x_healpix*fx +
                                     y_healpix*fy +
                                     z_healpix*fz ) * ds

            # self.norm[halo_num]=np.sqrt( np.power(fx[idx],2) + np.power(fy[idx],2) + np.power(fz[idx],2) )

            # self.cos[halo_num]= flux_data[halo_num] / self.norm[halo_num]
            #
            # self.cosmean[halo_num]=np.mean(self.cos[halo_num])
            # self.cosstd[halo_num]=np.std(self.cos[halo_num])

            #consider only outflow
            if getonly >0:
                flux_data[halo_num][flux_data[halo_num]<0]=0

            #consider only inflow
            if getonly <0:
                flux_data[halo_num][flux_data[halo_num]>0]=0

            # tot flux
            tot_flux[halo_num]=np.sum(flux_data[halo_num])
            # tot_flux[halo_num]=np.max(flux_data[halo_num])

        print("skipped %d/%d = %.02f %%"%(skipped,self.nfoftot, skipped/self.nfoftot*100))

        if Rmin is not None:
            print("Halos smaller than Rmin %d/%d = %.02f %% Rmin=%e"%(N_halo_smaller_than_Rmin,self.nfoftot, N_halo_smaller_than_Rmin/self.nfoftot*100, Rmin))


        setattr(self,"flux_data_%s_%s"%(type,fact_name),flux_data)
        setattr(self,"tot_flux_%s_%s"%(type,fact_name),tot_flux)

        with open(name, 'wb') as output:
            pickle.dump(flux_data, output,-1)
            pickle.dump(tot_flux, output,-1)

        print("get_flux_r200 done")

####################################################################################################################
####################################################################################################################
    def get_instant_SFR(self, force=0, tshift=None):
        """
        compute the averaged SFR over the last 10Myr by getting the mass of recent stellar particles
        """

        try:
            ages=self.step.t - self.step.star.age.data
            masses=self.step.star.mass.data
        except AttributeError:
            #case where snap do not get stars
            self.instant_SFR=None
            return

        from pyemma import io, time
        info = io.Info(self.path+"../../../")

        self.instant_SFR=np.zeros(self.nfoftot)
        for halo_num in range(self.nfoftot):
            star_mask=self.stars[halo_num]

            age=ages[star_mask]

            #To compute escape fraction age have to be shift by R200/c
            if tshift is not None:
                age-=tshift[halo_num]

            mass=masses[star_mask]

            thresh=1e7# thresh should be lower than tlife_feedback to neglect mass return

            self.instant_SFR[halo_num]=np.sum( mass[age<thresh] )
            self.instant_SFR[halo_num]*=info.unit_mass/1.9891e30
            self.instant_SFR[halo_num]/=thresh



    def get_instant_SFR_from_gas(self):
        """
        compute theorical instant SFR using e*rho/tff
        """

        step=self.step

        from pyemma import io, time
        info = io.Info(self.path+"../../../")
        run = io.RunParam(self.path+"../../../")

        rho=step.grid.field_d.data *info.unit_mass /  np.power(step.a *info.unit_l,3, dtype=np.float128) # kg.m-3

        n_halo = self.nfoftot
        self.instant_SFR_from_gas = np.zeros(n_halo)
        for ihalo in range(n_halo):

            if len(self.cells[ihalo])>0:

                dx=np.power(0.5,step.grid.l.data[self.cells[ihalo]]) *step.a *info.unit_l # m
                dv=np.power(dx,3, dtype=np.float128) # m3

                rho_halo=rho[self.cells[ihalo]] # kg.m-3
                tff=np.sqrt(3*np.pi/(32*6.67384e-11*rho_halo)) # s

                mask = rho_halo>run.overdensity_cond *info.unit_mass /  np.power(step.a *info.unit_l,3, dtype=np.float128)

                sfr_all = run.efficiency*rho_halo/tff *dv

                sfr = np.sum( sfr_all[mask] ) # kg.s-1
                sfr *= 1./1.9891e30 *(365*24*3600) # Mo.yr-1

                self.instant_SFR_from_gas[ihalo]=sfr

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
