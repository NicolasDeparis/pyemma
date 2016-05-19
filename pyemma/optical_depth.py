# coding: utf-8
import numpy as np

"""
For more information about this file see Chardin et al 2015
http://adsabs.harvard.edu/abs/2015MNRAS.453.2943C
This is an adaptation of Chardin's code
"""

class OpticalDepth:
    def __init__(self):
        pass

    def _LOS (self, x, y, z, l, d, t, vz, xion) :
        """
        Line Of Sight
        Récupère l'information d'une ligne de visée
        """

        l_coarse=np.min(l)
        maskl = [l==l_coarse]
        x_rand = np.random.choice(x[maskl])
        y_rand = np.random.choice(y[maskl])
        mask = [x==x_rand, y==y_rand]
        mask2 = mask[0] & mask[1]

        z_new = list(z[mask2])
        l_new = list(l[mask2])
        d_new = list(d[mask2])
        t_new = list(t[mask2])
        vz_new = list(vz[mask2])
        xion_new = list(xion[mask2])

        l_new.sort(key=dict(zip(l_new, z_new)).get)
        d_new.sort(key=dict(zip(d_new, z_new)).get)
        t_new.sort(key=dict(zip(t_new, z_new)).get)
        vz_new.sort(key=dict(zip(vz_new, z_new)).get)
        xion_new.sort(key=dict(zip(xion_new, z_new)).get)

        return l_new, d_new, t_new, vz_new, xion_new


    def _THR (self, x_init, y_init, z_init, l_init, d_init, t_init, vz_init, xion_init, l_coarse) :
        """
        To Higher Resolution
        Permet de passer de l=7 à l=10 (pour avoir toute l'information sur une ligne de visée)
        """
        l, d, t, vz, xion = self._LOS(x_init, y_init, z_init, l_init, d_init, t_init, vz_init, xion_init)

        lmax = np.max(l_init)
        n=2**(lmax)
        d_new = np.zeros(n)
        t_new = np.zeros(n)
        vz_new = np.zeros(n)
        xion_new = np.zeros(n)

        idx=0
        for i in range(len(l)) :
            for j in range(int(pow(2,lmax-l[i]))) :
                d_new[idx] = d[i]
                t_new[idx] = t[i]
                vz_new[idx] = vz[i]
                xion_new[idx] = xion[i]

                idx+=1

        return d_new, t_new,vz_new, xion_new

    def get_LOS(self,step, nline):
        self.nline = nline

        self.LOS_d=np.zeros(nline, dtype=np.object)
        self.LOS_t=np.zeros(nline, dtype=np.object)
        self.LOS_vz=np.zeros(nline, dtype=np.object)
        self.LOS_xion=np.zeros(nline, dtype=np.object)

        lmin = np.min(step.grid.l.data)

        for i in range(self.nline) :
            self.LOS_d[i], self.LOS_t[i], self.LOS_vz[i], self.LOS_xion[i] = self._THR(
                step.grid.x.data,
                step.grid.y.data,
                step.grid.z.data,
                step.grid.l.data,
                step.grid.field_d.data,
                step.grid.rfield_temp.data,
                step.grid.field_w.data,
                step.grid.xion.data,
                lmin)

    def get_tau(self, run, step):

        oM = run.param.info.om
        oV = run.param.info.ov
        oB = run.param.info.ob
        aexp = step.a[0]

        H0 = run.param.info.H0
        h = H0/100
        taille_boite = run.param.info.box_size_hm1_Mpc/h

        lmin=run.param.info.level_min
        lmax=np.max(step.grid.l.data)

        nb_seg = int(2**lmax)

        c_light = 299792458
        mH = 1.672623e-27
        kB = 1.3806488e-23
        Mpc = 3.085677581e22

        # sigma_alpha = 4.48e-22 #theuns 1998

        sigma_thomson=6.625e-21 #theuns 1998

        sigma_thomson=2.35188329313e-22 # EMMA 1grp sigma E
        # sigma_thomson=1.82867991612e-22 # EMMA 1grp sigma I

        f = 0.41615 #oscillator strength
        lambda0 = 1215.6 *1e-10 # H lyman alpha transition, angstrom -> m
        sigma_alpha=np.sqrt(3.*np.pi*sigma_thomson/8.)*f*lambda0

        print(sigma_alpha)

        taille_seg = taille_boite * aexp / nb_seg  # Mpc
        H_z = H0 * np.sqrt((1-oM-oV)/(aexp**2) + oM/(aexp**3) + oV) # km s-1 Mpc-1 \\* 3.2407e-20 (s-1)
        distance = np.linspace(0,nb_seg-1,nb_seg)
        v_H = H_z * distance * taille_seg# km s-1
        v_H *= 1000 # m s-1
        K = sigma_alpha*c_light*taille_seg*Mpc/np.sqrt(np.pi) # m^4 s-1

        self.LOS_d *= run.param.info.unit_mass / np.power(run.param.info.unit_l*aexp,3) / mH # atome H m-3
        nHI = self.LOS_d * (1-self.LOS_xion) # atome H neutre m-3
        self.LOS_vz /= run.param.info.unit_v / aexp # m s-1

        self.teff=np.zeros(self.nline)
        self.tau_liste=np.zeros((self.nline,nb_seg))

        for k in range(self.nline):
            if (k+1)%(100) == 0 :   # Affiche toutes les 100 lignes de visées réalisées
                print(k+1, '-> Done..')

            vtot = v_H + self.LOS_vz[k]
            bHI = np.sqrt(2*kB*self.LOS_t[k]/mH) # m s-1

            # need a floor value to avoid divide by 0
            bHI[bHI==0]=1e-5

            for i in range(nb_seg) :
                gauss = (nHI[k]/bHI) * np.exp(-((v_H[i]-vtot)/bHI)**2)
                self.tau_liste[k][i] = K*np.sum(gauss)

            self.teff[k] = -np.log(np.mean(np.exp(-self.tau_liste[k])))
