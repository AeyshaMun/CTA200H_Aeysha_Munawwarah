a0Warning = 'WARNING! Time necessary for binary to shrink from a0 to current_a is shorter than tflightmax. Number of ejected stars may not be accurate. Consider increasing a0 in BMBH() call'

cutWarning = 'WARNING! Ejected star sample size is being cut down from its original size. Please ensure this is intentional'

metError = 'Error: invalid input for metallicity argument'

imfError = 'Error: invalid input for initial mass function'

BMBHMergerStateError = 'Error: "a_current" parameter or "tlb" parameter must be argument to BMBH() but not both. See speedystar.eject.BMBH() docstring for help'

metLowWarning = 'Warning: Some assigned metallicities are too low. SSE module within AMUSE only provides stellar evolution formulae for stars in the metallicity range [0.0001, 0.03]. Metallicites below 0.0001 have been set to 0.0001'

metHighWarning = 'Warning: Some assigned metallicities are too high. SSE module within AMUSE only provides stellar evolution formulae for stars in the metallicity range [0.0001, 0.03]. Metallicites above 0.03 have been set to 0.03'


#try:
import numpy as np
from tqdm import tqdm
import astropy.constants as const
from astropy import units as u
from astropy.table import Table
import speedystar

from .utils.mainsequence import t_MS
from math import ceil
#import emcee
import importlib
import sys
import time
import os

from .utils.imfmaster.imf import imf

class EjectionModel:
    '''
        Ejection model class
    '''

    _name = 'Unknown'

    def __init__(self, name):
                self._name      = name

    def sampler(self):
        '''
        Sampler of the ejection distribution
        '''
        raise NotImplementedError

    def _evo_star(self,mass,age=None,met=0.0):
        from amuse.community.sse.interface import SSE
        from amuse import datamodel
        from amuse.units import units
        #Evolve a star of a certain mass and metallicity until a certain age 
        #using the SSE module within AMUSE

        #Initialize
        stellar_evolution = SSE()
        #Adjust metallicity for new Zsun assumption
        stellar_evolution.parameters.metallicity = self.Zsun*10**(met)
        star      = datamodel.Particles(len(mass))
        star.mass = mass | units.MSun

        if age is not None:
            age = age.to('Myr').value | units.Myr

        #Evolve the star
        star = stellar_evolution.particles.add_particles(star)
        stellar_evolution.commit_particles()
        stellar_evolution.evolve_model(end_time = age)

        stellar_evolution.stop()

        return star

    def set_imf(self):
        '''
        Sets the initial mass function which governs ejected star masses
        ''' 

        #Single power law
        if self.imftype=='Salpeter':
            return imf.Salpeter(alpha=self.kappa,mmin=self.mmin,mmax=self.mmax)

        #Broken power law with break at m=mbreak
        elif self.imftype=='Kroupa':
            if (self.mmin < self.mbreak1) and (self.mmax <= self.mbreak1):
                return imf.Salpeter(alpha=self.kappa1, mmin=self.mmin, mmax=self.mmax) 
            elif (self.mmin < self.mbreak1) and (self.mmax > self.mbreak1) \
                    and (self.mmax <= self.mbreak2):
                return imf.Kroupa(mmin=self.mmin,mmax=self.mmax,
                                    p1=self.kappa1,p2=self.kappa2,p3=self.kappa2, 
                                    break1=self.mbreak1,break2=(self.mbreak1+self.mmax/2))
            elif (self.mmin < self.mbreak1) and (self.mmax > self.mbreak2) :
                return imf.Kroupa(mmin=self.mmin,mmax=self.mmax,
                                    p1=self.kappa1,p2=self.kappa2,p3=self.kappa, 
                                    break1=self.mbreak1,break2=self.mbreak2)

            elif (self.mmin >= self.mbreak1) and (self.mmin < self.mbreak2) \
                    and (self.mmax > self.mbreak1) and (self.mmax <= self.mbreak2):
                 return imf.Salpeter(alpha=self.kappa2, mmin=self.mmin, mmax=self.mmax)           

            elif (self.mmin >= self.mbreak1) and (self.mmin < self.mbreak2) \
                    and (self.mmax > self.mbreak1) and (self.mmax > self.mbreak2):
                return imf.Kroupa(mmin=self.mmin,mmax=self.mmax,
                                    p1=self.kappa2,p2=self.kappa2,p3=self.kappa, 
                                    break1=(self.mmin + self.mbreak2)/2, break2=self.mbreak2)
            elif (self.mmin >= self.mbreak2) and (self.mmax >= self.mbreak2):
                return imf.Salpeter(alpha=self.kappa, mmin=self.mmin, mmax=self.mmax)           
            #elif self.mmin < self.mbreak2:
            #    return imf.Kroupa(mmin=self.mmin,mmax=self.mmax,
            #                        p1=self.kappa2,p2=self.kappa2,p3=self.kappa, 
            #                        break1=(self.mmin+self.mbreak2)/2,break2=self.mbreak2)
            #else:
            #    return imf.Kroupa(mmin=self.mmin,mmax=self.mmax,
            #                        p1=self.kappa,p2=self.kappa,p3=self.kappa, 
            #                        break1=(self.mbreak2+self.mmax)/2,break2=self.mbreak2)

        #User-specified imf
        elif callable(self.imftype):
            return self.imftype

        else:
            raise ValueError(imfError)

    def _inverse_cumulative_mp_old(self, x):
        #Inverse sampling of simple power law IMF between 0.1 and 100 Msun
        #Depreciated
        mmin = 0.1
        mmax = 100.
        if self.kappa==-1:
            return mmin*(mmax/mmin)**x
        else:
            return (  (mmax**(1.+self.kappa) - mmin**(1.+self.kappa))*x 
                    + mmin**(1.+self.kappa) )**(1./(1.+self.kappa)) 

    def _inverse_cumulative_mp_old2(self, x):
        #Inverse sampling of a Kroupa IMF between 0.1 and 100 Msun
        #Depreciated

        '''
            Inverse of the cumulative function of the Kroupa IMF f(mp) as a function of x.

            F(mp) = int_[0.1, 100] f(mp)

            returns mp such that

            F(mp) = x

        '''

        x = np.atleast_1d(x)

        F_cut = 0.729162 # F(x) at the breaking point
        total_int = 6.98626 # total integral of the broken power law between 0.1 and 100

        #F_cut = 0. # F(x) at the breaking point
        #total_int = 1.8941 # total integral of the broken power law between 0.1 and 100


        idx_cut = x<F_cut
        result = np.zeros(x.shape)

        result[idx_cut] = ((0.1)**(-0.3)-x[idx_cut]*0.3*0.5*total_int)**(1./(-0.3))
        result[~idx_cut] = (-total_int * (x[~idx_cut] - F_cut)*(1.3) + (0.5)**(-1.3))**(-1./1.3)
        return result


    def _inverse_cumulative_mp(self, IMF, x):
        '''
            Inverse of the CDF of a single power law IMF as a function of x
        '''

        #from .utils.imfmaster.imf import imf

        #testimf = imf.Salpeter(alpha=self.kappa,mmin=0.1,mmax=100)
        result = imf.inverse_imf(x,mmin=IMF.mmin,mmax=IMF.mmax,massfunc=IMF)
        return result#*u.Msun

    def _inverse_cumulative_mp_broken(self, x):
        '''
            Inverse sampling function for a Kroupa IMF. Depreciated.
            Inverse of the CDF of a broken power law IMF as a function of x
            By default, breaks and low/intermediate regime slopes are set as 
            in Kroupa, high-end slope can be an argument to __init__().
            A caveat: if m_range[0] is less than break1 
            (as is the case by default), this class throws an error.
            You can avoid that by saying p1=p2 and setting break1 to anywhere 
            in between mmin and break2, as in the default class call here

        '''
        from imfmaster.imf import imf #%or wherever else you installed it to

        myimf = imf.Kroupa(mmin=0.1,mmax=100,p1=1.3,p2=1.3,p3=self.kappa,
                            break1=0.2,break2=0.5)

        #Uncomment to get a proper Kroupa IMF
        #myimf = imf.Kroupa(mmin=0.05,mmax=900,p1=0.3,p2=1.3,p3=self.kappa,
                            #break1=0.08,break2=0.5)

        result = imf.inverse_imf(x,massfunc=myimf)
        return result

    def sample_met(self, n):
        '''
            Assigns a total stellar metallicity to each mock ejected star.
        '''        
        
        #Sample metallicities

        #If self.Met is a float, all stars have the same metallicity
        if isinstance(self.Met,float):
            met = np.full(n, self.Met)
        elif isinstance(self.Met,int):
            met = np.full(n, float(self.Met))
        #If self.Met is 'uniform', metallicities are sampled uniformly in the
        # range [-0.25, 0.25] as in Evans et al. 2022 (MNRAS, 517, 3469)
        elif self.Met=='uniform':
            met = np.round(np.random.uniform(-0.25, 0.25, len(mp)), decimals=2)
        elif callable(self.Met):
            met = np.round(self.Met(size=n, *self.metargs),decimals=2)

        else:
            raise ValueError(metError)

        if any(self.Zsun*10.**met < 0.0001):
            print(metLowWarning)
            
            idx = self.Zsun*10.**met < 0.0001
            met[np.where(idx)[0]] = np.log10(0.0001 / self.Zsun)

        if any(self.Zsun*10.**met > 0.03):
            print(metHighWarning)
            
            idx = self.Zsun*10.**met > 0.03
            met[np.where(idx)[0]] = np.log10(0.03 / self.Zsun)            

        return met

class Hills(EjectionModel):

    '''
    HVS ejection model from Rossi+ 2017. Isotropic ejection from 3 pc from GC 
    and with a mass/velocity distribution based on MC. Can generate an 
    ejection sample using a Monte Carlo approach based on inverse transform 
    sampling. See also Marchetti+ 2018, Evans+2021, Evans+2022.
    '''

    T_MW = 13.8*u.Gyr # MW maximum lifetime from Planck2015

    #Initial ejection radius from GC
    centralr = 3*u.pc

    def __init__(self, M_BH = 4e6*u.Msun, Met = np.random.uniform, Zsun = 0.0142,
                metargs = (-0.25,+0.25), tflightmax = 100.*u.Myr, cut_factor = 1,
                name_modifier = None, alpha=-1, gamma=0, rate=1e-4/u.yr, 
                imftype = 'Salpeter', kappa1 = 0.3, kappa2 = 1.3, kappa = 2.3,
                mmin = 0.1, mmax = 100, mbreak1 = 0.08, mbreak2 = 0.5,
                v_range = [500, 5e4]*u.km/u.s, m_range = [0.5, 1e3]*u.Msun,
                amuseflag = False):

        '''
        Parameters
        ----------
        name_modifier : str
            Add a modifier to the class name

        M_BH: astropy Msun quantity
            Change mass of Sgr A*

        Zsun : float
            Total metallicity of the sun. Default is 0.0142 (Asplund+2009) 
            0.02 is another common choice (Anders & Grevasse 1989)

        Met: Float or callable function
            If float, all ejected stars have total metallicity 
            xi = log10(Z/Zsun) = Met.
            If callable, mtallicities are assigned by calling the function.
            Default is uniform distribution.
            Arguments to this callable can be provided in order as metargs.
            E.g. if you want metallicites to be a Gaussian distribution centred
            on 0.0 with a standard deviation of 1.0, call ejection model with 
            (Met = np.random.normal, metargs=[0.0, 1.0])

        metargs: tuple or listlike 
            Arguments, in order, to the supplied metallicity distribution, 
            see above.

         tflightmax : float
            Stars with a flight time more than tflightmax are tossed out 
            and aren't propagated. Provides a lot of speedup, as long as you're
            certain stars with tflight<tflightmax are not relevant for 
            your science case

        cut_factor: integer
            Final sample will be cut down by a factor of cut_factor, 
            i.e. size = size / cut_factor. 
            E.g. set cut_factor to 10 means the returned sample will be 
            10 times smaller than it should be. Useful for situations where 
            you know you'll get way too many stars.

        alpha : float
            Exponent of the power-law for the distribution of the 
            semi-major axes in binaries. Default is -1

        gamma : float
            Exponent of the power-law for the distribution of the mass ratio 
            in binaries. Default is 0

        rate : astropy quantity (yr^-1)
            Assumed ejection rate from the SMBH

        imftype: 'Salpeter' or 'Kroupa' or callable.
            Functional form of initial mass function. If 'Salpeter', IMF is 
            single power law. If 'Kroupa', IMF is broken power law. Power law 
            slopes are NOT necessarily assumed to follow Salpter/Kroupa slopes.
            Slope is set by kappa argument (see below). 
            Alternatively, a complete IMF defined by the imf package can 
            be passed explicitly and stellar masses are drawn from it.

        kappa: float
            log-slope of the IMF (if imftype='Salpeter') or log-slope of the 
            m>mbreak2 regime of the IMF (if imftype=='Kroupa'). Irrelevant
            if imftype is a callable,

        mmin, mmax: floats
            Smallest and largest stellar mass formed by the IMF.
            Irrelevant if imftype is a callable.
            *NOT the same as m_range, see below.

        mbreak1, mbreak2: floats
            Stellar mass where the IMF changes slope if imftype=='Kroupa'. 
            Irrelevant if imftype=='Salpeter' or imftype is a callable.

        kappa1, kappa2: floats
            log-slopes of the IMF in the m<mbreak1 and mbreak1 < m < mbreak2 regimes.
            Irrelevant if imftype=='Salpeter' or imftype is a callable.


        v_range : list of astropy quantities (km s^-1)
            Stars drawn with initial velocities outside this range are
            immediately discarded

        m_range : list of astropy quantities (Msun)
            Stars sampled with masses outside this range are immediately discarded.
            Can speed up the code a lot since the lifetimes of the low-mass stars
            won't be calculated.

        amuse_flag : Boolean
            If True, stellar evolution and maximum lifetimes are handled
            using the AMUSE python package. Allows more detailed stellar evolution
            past the main sequence, but installation can be a headache
            If False, stellar evolution is handled through in-house functions.
            Cannot handle evolution past the end of the main sequence, which may
            or may not be an issue depending on use case.

        '''

        self.alpha      = alpha
        self.gamma      = gamma
        self.M_BH       = M_BH
        self.rate       = rate

        self.mmin       = mmin
        self.mmax       = mmax

        if mbreak1>=mbreak2:
            #raise ValueError(mBreakError)
            raise ValueError('Error: Lower mass break in the two-break power law IMF is set higher than the upper mass break.')
        if mmin >= mmax:
            raise ValueError('Error: mmin argument must be smaller than mmax argument')
        if mmin > m_range[0].value:
            print('WARNING -- IMF is not generating stars down to the lower mass cut')
        if mmax > m_range[1].value:
            print('WARNING -- IMF is not generating stars up to the upper mass cut')
        
        self.mbreak1    = mbreak1
        self.mbreak2    = mbreak2
        self.kappa1     = kappa1
        self.kappa2     = kappa2
        self.kappa      = kappa

        self.tflightmax = tflightmax
        self.cut_factor = cut_factor

        self.v_range    = v_range
        self.m_range    = m_range
        self.Met        = Met
        self.metargs    = metargs
        self.imftype    = imftype
        self.Zsun       = Zsun

        self.amuseflag  = amuseflag
        
        if(name_modifier is not None):
            self._name = 'Hills - ' + name_modifier
        else:
            self._name = 'Hills'

    def _inverse_cumulative_q(self, x, mp):
        '''
           Inverse of the CDF of a single power law mass ratio distribution 
            as a function of x
        '''

        qmin = 0.1
        qmax = 1.

        if self.gamma==-1:
            return qmin*(qmax/qmin)**x
        else:
            q = (  (qmax**(1.+self.gamma) - qmin**(1.+self.gamma))*x \
                    + qmin**(1.+self.gamma) )**(1./(1.+self.gamma))
            return q 

    def _inverse_cumulative_logP(self,x,Rmax,mtot):

        '''
           Inverse of the CDF of a single power law log-period 
            distribution as a function of x
        '''

        #Pmin is set by Roche lobe overflow, Pmax is (arbitrarily) 2000 Rsun
        Pmin = np.log10(np.sqrt((4*np.pi**2 * (2.5*Rmax)**3) \
                        / (const.G*mtot)).to('second').value)
        Pmax = np.log10(np.sqrt((4*np.pi**2 * (2e3*u.Rsun)**3) \
                        / (const.G*mtot)).to('second').value)

        if self.alpha==-1:
            return (10**(Pmin*(Pmax/Pmin)**x)) * u.second
        else:
            return (  10**( ((Pmax**(1.+self.alpha) - Pmin**(1.+self.alpha))*x \
                    + Pmin**(1.+self.alpha) )**(1./(1.+self.alpha))))*u.second

    def sampler(self):
        '''
        Samples from the ejection distribution to generate an ejection sample.
        The distribution mass and velocity is generated using a Monte Carlo 
        approach (Rossi 2014). The functions _inverse_* dictate the 
        parameters of the progenitor binary population. 
        The velocity vector is assumed radial.

        Returns
        -------
            r0, phi0, theta0, v0, phiv0, thetav0 : Quantity
                Initial phase space position in spherical coordinates, 
                centered on the GC
            m, tage, tflight : Quantity
                Stellar mass of the HVS, age at observation and tflight 
                between ejection and observation      
            stage, stagebefore : int
                Evolutionary stage (e.g. main sequence, red giant) of your HVS 
                *today* and at the moment of ejection.
                Stage conventions follow Hurley et al. (2000) 
                https://ui.adsabs.harvard.edu/abs/2000MNRAS.315..543H/abstract             
            n : int
                Size of the output ejection sample

        '''
    
        import time
        from math import ceil

        #try:
        from astropy import constants as const

        from .utils.hurley_stellar_evolution import get_t_BAGB, get_t_MS
        from .utils.hurley_stellar_evolution import Radius, get_TempRadL

        #print('snark2')
        #from amuse.community.sse.interface import SSE
        #from amuse import datamodel

        # Sample the binary properties q, a, mp using inverse sampling
        #Take only stars ejected within the last tflightmax Myr

        n = np.rint((self.tflightmax*self.rate).decompose()).value.astype(int)

        #Cut down sample size if required
        n = int(n/self.cut_factor)

        if(self.cut_factor>1):
            print(cutWarning)

        tflight = np.random.uniform(0., self.tflightmax.to('Myr').value, n) \
                    *u.Myr

        # Inverse sampling a primary mass and mass ratio
        uniform_for_mp, uniform_for_q = np.random.uniform(0, 1, (2, n))
        #mp = self._inverse_cumulative_mp2(uniform_for_mp)
        
        IMF = self.set_imf()
        mp = self._inverse_cumulative_mp(IMF=IMF, x=uniform_for_mp)

        #mp = self._inverse_cumulative_mpold3(uniform_for_mp)
        q  = self._inverse_cumulative_q(uniform_for_q, mp)

        mp = mp*u.Msun

        #Randomly designate one member of the binary as the HVS and 
        #the other as the companion (C)
        ur  = np.random.uniform(0,1,n)
        idx = ur>=0.5

        #mem is 1 if the ejected HVS is the heavier star of the binary, 
        #2 if the lighter.
        #Generally not super helpful but can be good for debugging
        mem_HVS       = np.zeros(n) 
        mem_HVS[idx]  = 1
        mem_HVS[~idx] = 2

        M_HVS, M_C    = np.zeros(n)*u.Msun, np.zeros(n)*u.Msun
        M_HVS[idx]    = mp[idx]
        M_HVS[~idx]   = mp[~idx]*q[~idx]
        M_C[idx]      = mp[idx]*q[idx]
        M_C[~idx]     = mp[~idx]

        idx        = (M_HVS > self.m_range[0]) & (M_HVS < self.m_range[1])
        tflight, q = tflight[idx], q[idx]
        m, mc, mem = M_HVS[idx], M_C[idx], mem_HVS[idx]

        n = idx.sum()

        #Solar metallicity
        #met = np.zeros(n) 

        #Uncomment this line to give each HVS a metallicity xi = log10(Z/Zsun) 
        #randomly distributed between -0.25 and +0.25 
        #met = np.round(np.random.uniform(-0.25, 0.25, n),decimals=2)       

        met = self.sample_met(n)

        if self.amuseflag:
            #Calculate stellar ages and remove stars that are remnants today
            #Calculating maximum lifetimes (taken here as the time from ZAMS 
            #to the beginning of the AGB branch) takes a while.
            #Calculating main sequence lifetimes is much faster, can save time by 
            #pre-selecting using the MS lifetime.
            #We cut out stars older than 1.4 times their main sequence lifetime
            #(they are definitely dead)
            #For the rest, we calculate proper maximum lifetimes

            #Calculate the main sequence lifetime of the heavier star in the 
            #binary using the Hurley+2000 relation. These relations assume 
            #Zsun=0.02, so might need to be adjusted if self.Zsun is different

            T_maxbig = get_t_MS(np.maximum(m,mc),met+np.log10(self.Zsun/0.02))
            T_maxbig[T_maxbig>self.T_MW] = self.T_MW

            #Binary hung around the GC for a random fraction e1 
            #of its maximum lifetime
            #("maximum lifetime" here is taken as 1.4*t_MS, see above)
            e1       = np.random.random(n)

            t_before = 1.4*T_maxbig*e1

            #Calculate time remaining for the binary
            t_rest   = 1.4*T_maxbig*(1-e1)

            #Calculate the HVSs maximum lifetime. This may or may not be the same 
            #as the binary's maximum lifetime. 
            #It's longer if the HVS is the smaller member
            T_maxHVS       = T_maxbig.copy()
            T_maxHVS[m<mc] = get_t_MS(m[m<mc],met[m<mc]+np.log10(self.Zsun/0.02))
            T_maxHVS[T_maxHVS>self.T_MW] = self.T_MW

            #Calculate HVS's age
            tage = np.zeros(len(tflight))*u.Myr
            tage = t_before + tflight

            #Calculate the maximum lifetime of each star
            #Remove HVS from sample if its age is >1.4 times 
            #its main sequence lifetime
            tage[tage>=1.4*T_maxHVS] = np.nan
            idx = (~np.isnan(tage))

            m, mc, tage          = m[idx], mc[idx], tage[idx]
            tflight, q, mem, met = tflight[idx], q[idx], mem[idx], met[idx]
            e1                   = e1[idx]

            n = idx.sum()

            #Redo the last few steps, but now calculating actual 
            #proper lifetimes for the stars
            T_maxbig = get_t_BAGB(np.maximum(m,mc),met+np.log10(self.Zsun/0.02))
            #T_maxbig = get_t_MS(np.maximum(m,mc),met+np.log10(self.Zsun/0.02))
            T_maxbig[T_maxbig>self.T_MW] = self.T_MW

            T_maxHVS = T_maxbig.copy()
            T_maxHVS[m<mc] = get_t_BAGB(m[m<mc],met[m<mc]+np.log10(self.Zsun/0.02))        
            #T_maxHVS[m<mc] = get_t_MS(m[m<mc],met[m<mc]+np.log10(self.Zsun/0.02))        
            T_maxHVS[T_maxHVS>self.T_MW] = self.T_MW

        else:
            e1       = np.random.random(n)
            #Redo the last few steps, but now calculating actual 
            #proper lifetimes for the stars
            #T_maxbig = get_t_BAGB(np.maximum(m,mc),met+np.log10(self.Zsun/0.02))
            T_maxbig = get_t_MS(np.maximum(m,mc),met+np.log10(self.Zsun/0.02))
            T_maxbig[T_maxbig>self.T_MW] = self.T_MW

            T_maxHVS = T_maxbig.copy()
            #T_maxHVS[m<mc] = get_t_BAGB(m[m<mc],met[m<mc]+np.log10(self.Zsun/0.02))        
            T_maxHVS[m<mc] = get_t_MS(m[m<mc],met[m<mc]+np.log10(self.Zsun/0.02))        
            T_maxHVS[T_maxHVS>self.T_MW] = self.T_MW           

        t_before = T_maxbig*e1
        t_rest = T_maxbig*(1-e1)

        tage = np.zeros(len(tflight))*u.Myr
        tage = t_before + tflight
        tage[tage>=T_maxHVS] = np.nan        

        idx = (~np.isnan(tage))

        m, mc                      = m[idx], mc[idx], 
        tage, t_before, tflight    = tage[idx], t_before[idx], tflight[idx]
        q, mem, met                = q[idx], mem[idx], met[idx]
        n                          = idx.sum()

        #Need the larger star's size to draw an orbital period

        #Get the mass of the larger star in the binary
        mbig = m.copy()
        mbig[mc>m] = mc[mc>m]

        Rmax, stagebefore, R, Lum, stage, T_eff = (np.empty(n) \
                                                for i in range(6))

        if self.amuseflag:

            #Stars with different metallicities can't be evolved at the same time 
            #with AMUSE (or so it seems to me), so need to do them in batches
            #Loop over unique metallicities
            #Should update, could be sped up
            for z in tqdm(np.unique(met),desc='Evolving HVSs'):   

                #indices with ith metallicity, subset
                idx = np.where(met==z)[0]

                #Evolve larger star in each of the binaries
                star = self._evo_star(mbig[idx].value,t_before[idx],met=z)
            
                #Extract bigger star radii and stages
                Rmax[idx] = star.radius.as_astropy_quantity().to('Rsun').value
                stagebefore[idx] = star.stellar_type.as_astropy_quantity()

                #Evolve each HVS
                star = self._evo_star(m[idx].value,tage[idx],met=z)

                #Extract HVS effective temperature, radius, 
                #luminosity, evolutionary stage
                T_eff[idx] = star.temperature.as_astropy_quantity().to('K').value
                R[idx] = star.radius.as_astropy_quantity().to('Rsun').value 
                Lum[idx] = star.luminosity.as_astropy_quantity().to('Lsun').value 
                stage[idx] = star.stellar_type.as_astropy_quantity()
        else:

                #indices with ith metallicity, subset
                #idx = np.where(met==z)[0]
            
                #Rmaxtmp = np.empty(len(mbig[idx]))
                #Rtmp = np.empty(len(mbig[idx]))
                #T_efftmp = np.empty(len(mbig[idx]))
                #Lumtmp = np.empty(len(mbig[idx]))

                #Extract bigger star radii
                for i in tqdm(range(len(mbig)), desc='Evolving HVSs'):   
                    #for z in tqdm(np.unique(met),desc='Evolving HVSs'):   
                    #Get radius of bigger star
                    Rmax[i] = Radius(mbig[i].value,met[i], 
                                        t_before[i].value)

                    #Rmax[idx] = Rmaxtmp
                    stagebefore[i] = 1

                    #Get HVS effective temperature, radius, luminosity
                    #for i in range(len(m[idx])):   
                    T_eff[i], R[i], Lum[i] = get_TempRadL(m[i].value,
                                                       met[i], tage[i].value)
                    stage[i] = 1
            
                #R[idx] = Rtmp
                #T_eff[idx] = T_efftmp
                #Lum[idx] = Lumtmp            
                #stage[idx] = 1
    
                
        Rmax = Rmax*u.Rsun
        R = R*u.Rsun
        Lum = Lum*u.Lsun
        T_eff = T_eff*u.K

        #Inverse sample an orbital period for the binary
        uniform_for_P = np.random.uniform(0, 1, n)
        P = self._inverse_cumulative_logP(uniform_for_P,Rmax,m+mc)
        P = P.to('day')

        #Convert this to an orbital separation, Kepler's 3rd law
        sep = np.cbrt((const.G*(m+mc)*P**2)/(4*np.pi**2)).to('Rsun')

        #Calculate velocity, see Sari+2010, Kobayashi+2012, Rossi+2014
        V_HVS = (np.sqrt( 2.*const.G.cgs*mc / sep )  \
                * ( self.M_BH/(mc+m) )**(1./6.)).to('km/s')

        #Cut out stars that are outside the imposed velocity range
        idx = (V_HVS > self.v_range[0]) & (V_HVS < self.v_range[1])

        m, mc                   = m[idx], mc[idx] 
        tage, t_before, tflight = tage[idx], t_before[idx], tflight[idx] 
        q, mem, met             = q[idx], mem[idx], met[idx]
        v, P, sep               = V_HVS[idx], P[idx], sep[idx] 
        stagebefore, stage      = stagebefore[idx], stage[idx]
        T_eff, Lum, R           = T_eff[idx], Lum[idx], R[idx]

        n = idx.sum()

        #You can use this line below to cut stars that aren't relevant for you, 
        #e.g. stars that are too small, wrong evolutionary stage, etc.
        #Can save a lot of computation time
        idx = (v>0.*u.km/u.s)
        
        # Distance from GC at 3 pc
        r0 = np.ones(n)*self.centralr

        # Isotropic position unit vector in spherical coordinates
        phi0 = np.random.uniform(0,2*np.pi, n)*u.rad
        theta0 = np.arccos( np.random.uniform(-1,1, n))*u.rad

        # The velocity vector points radially.
        phiv0 = phi0.copy()
        thetav0 = theta0.copy()
    
        returns = {'r0':r0, 'phi0':phi0, 'theta0':theta0, 'v0':v, 'phiv0':phiv0, 'thetav0':thetav0, 'm':m, 'tage':tage, 'tflight':tflight, 'sep':sep, 'P':P, 'q':q, 'mem':mem, 'met':met, 'stage':stage, 'stagebefore':stagebefore, 'Rad':R, 'T_eff':T_eff, 'Lum':Lum, 'size':len(r0), 'alpha': self.alpha, 'gamma': self.gamma, 'kappa': self.kappa, 'M_BH': self.M_BH, 'tflightmax': self.tflightmax, 'rate': self.rate, 'v_range': self.v_range, 'm_range': self.m_range, 'Zsun': self.Zsun}

        return returns
                    
#####################################################

class BMBH(EjectionModel):

    '''
    Isotropic ejection from a 3 pc sphere centred on a 
    massive black hole binary (MBHB).
    '''

    T_MW = 13.8e3*u.Myr # MW maximum lifetime from Planck2015

    centralr = 3*u.pc # GC radius from which stars are ejected

    def __init__(self, name_modifier = None, Zsun = 0.0142, current_a=1e-3*u.pc,
                 tlb = np.nan, MC=4000*u.Msun, M_BH=4e6*u.Msun, a0 = 0.02*u.pc, 
                rho=7e4*u.Msun/(u.pc**3), sigma=100*u.km/u.s, 
                tflightmax = 50.*u.Myr, m_range = [0.5, 100]*u.Msun, 
                cut_factor=1, saveorbit=False, Met=np.random.uniform, 
                metargs=(-0.25,+0.25), vmin=550*u.km/u.s, imftype = 'Salpeter', 
                kappa = 2.35, orbitpath='./'):
        '''

        Parameters
        ----------
        name_modifier : str
            Add a modifier to the class name

        Zsun : float 
            Total solar metallicity,

        current_a : astropy quantity (distance)
            Present-day separation of MBHB

        tlb : astropy quantity (time)
            Time since the MBHB merged.
            * current_a OR tlb but NOT BOTH must be defined

        MC : astropy quantity (mass)
            Mass of companion BH

        M_BH : astropy quantity (mass)
            Mass of primary BH

        a0: astropy quantity (distance)
            Initial separation of MBHB

        rho: astropy quantity ( mass density)
            Stellar density of nuclear star cluster

        sigma: astropy quantity (velocity)
            Velocity dispersion of nuclear star cluster

        tflightmax: astropy quantity (time)
            Maximum flight time. Stars with a flight time longer than 
            this are automatically cut out

        m_range: length-2 list of astropy quantities (mass)
            Minimum and maximum possible masses of drawn stars

        cut_factor: integer
            Final sample will be cut down by a factor of cut_factor, 
            i.e. size = size / cut_factor. 
            E.g. set cut_factor to 10 means the returned sample will be 
            10 times smaller than it should be. Useful for situations where 
            you know you'll get way too many stars.

        save_orbit: Boolean
            If true, saves the orbital evolution of the MBHB and the mass 
            ejection history to file. Useful for diagnosis and/or debugging

        Met: Float or callable function
            If float, all ejected stars have total metallicity 
            xi = log10(Z/Zsun) = Met.
            If callable, mtallicities are assigned by calling the function.
            Default is uniform distribution.
            Arguments to this callable can be provided in order as metargs.
            E.g. if you want metallicites to be a Gaussian distribution centred
            on 0.0 with a standard deviation of 1.0, call ejection model with 
            (Met = np.random.normal, metargs=[0.0, 1.0])
        metargs: listlike 
            Arguments, in order, to the supplied metallicity distribution, 
            see above.

        vmin: astropy quantity (velocity)
            Minimum ejection velocity which counts as an 'ejected star'
            *Risky to change, hardening rate and mass ejection rate 
            fits change depending on this

        imftype: 'Salpeter' or 'Kroupa' or callable.
            Functional form of initial mass function. If 'Salpeter', IMF is 
            single power law. If 'Kroupa', IMF is broken power law. Power law 
            slopes are NOT necessarily assumed to follow Salpter/Kroupa slopes.
            Slope is set by kappa argument (see below). 
            Alternatively, a complete IMF defined by the imf package can 
            be passed explicitly and stellar masses are drawn from it.

        kappa: float
            log-slope of the IMF (if imftype='Salpeter') or log-slope of the 
            m>1.3*u.Msun regime of the IMF (if imftype=='Kroupa'). Irrelevant
            if imftype is a callable,

        

        '''

        from scipy.integrate import quad
        from scipy.interpolate import interp1d

        #Set some attributes based on arguments
        self.Zsun       = Zsun
        self.Met        = Met

        if current_a is not np.nan:
            self.current_a  = current_a
            self.tlb = np.nan
        if tlb is not np.nan:
            self.current_a = np.nan
            self.tlb = tlb
        if current_a is not np.nan and tlb is not np.nan:
            raise ValueError(BMBHMergerStateError)
        if current_a is np.nan and tlb is np.nan:
            raise ValueError(BMBHMergerStateError)
        
        #self.tlb        = tlb
        #self.current_a  = current_a
        self.cut_factor = cut_factor
        self.Mc         = MC
        self.kappa      = kappa
        self.M_BH       = M_BH
        self.q          = self.Mc / self.M_BH # BHB mass ratio
        self.a0         = a0
        self.rho        = rho
        self.sigma      = sigma 
        self.m_range    = m_range
        self.tflightmax = tflightmax
        self.saveorbit  = saveorbit
        self.orbitpath  = orbitpath
        self.vmin       = vmin
        self.imftype    = imftype
        self.metargs    = metargs
    
        #Load in hardening rate (H) and mass ejection rate (J) fits
        #approxs = Table.read('analytical_approximations.txt',format='ascii')

        approxs_file = importlib.resources.files(speedystar2).joinpath(
                            'utils/analytical_approximations.txt')
        approxs = Table.read(approxs_file,format='ascii')
        approxind = np.argmin(np.abs(approxs['q'] - self.q))


        self.fit_H = {'A':approxs['A_H'][approxind], 
                    'a0': approxs['a_0,H'][approxind], 
                    'gamma': approxs['gamma_H'][approxind] }
        self.fit_J = {'A':approxs['A_J'][approxind], 
                    'a0': approxs['a_0,J'][approxind], 
                    'gamma': approxs['gamma_J'][approxind], 
                    'alpha': approxs['alpha_J'][approxind], 
                    'beta': approxs['beta_J'][approxind] }
    
        if(name_modifier is not None):
            self._name = 'BMBH - ' + name_modifier
        else:
            self._name = 'BMBH'

    def calculate_vc(self, a):

        '''
        Calculated the binary orbital velocity

        Parameters
        -----------
        a: astropy quantity (length, can be array)
            Current separation of the binary

        Returns
        -----------
        vc: astropy quantity (velocity, can be array)
            Binary orbital velocity
        '''

        return (np.sqrt(const.G * self.M_BH * (1 + self.q) / a)).to("km/s")

    def v_CDF_glob(self):

        '''
        Generate a cumulative distribution function of HVS ejection velocities.
        See Eq. 11. This can be expensive to do over and over, so is saved 
        as a global variable in case ejection is being called many times 
        sequentially, e.g. in a loop

        Returns
        -----------
        x_range: numpy array
            range of v/v_c
        CDF: numpy array
            velocity CDF evaluated at each x
        '''

        makeCDF = True

        #If CDF already exists in globals, just grab it
        if 'x_rangeCDF' in globals() and 'vCDF' in globals():

            #Check to make sure companion mass and primary mass have 
            #remained unchanged
            if self.Mc == Mctmp and self.M_BH == MBHtmp:
                x_range = x_rangeCDF
                CDF = vCDF
                makeCDF = False

        #Make CDF if it doesn't already exist
        if makeCDF:

            #Define break velocity
            h = np.sqrt(2*self.q) / (1+self.q)
            wbreak = 1.2*h

            #Define velocity distribution, Eq. 11
            #Not distribution for v but distribution of w= v/v_c
            function = lambda x: 1 if x<wbreak else (x/wbreak)**-3.1
            xfunction = lambda x: function(x)/x
            
            #w = v/v_c. Set minimum and maximum w over which to integrate
            wmin, wmax = 1e-5, 10

            #define w array over which to integrate CDF
            x_range = np.logspace(np.log10(wmin), np.log10(wmax), 1000)

            # Calculate CDF normalization
            quad_min_max = quad(xfunction, wmin, wmax)[0]

            # Construct the CDF
            CDF = []
            lastval, lastx = 0, wmin
            for xr in x_range:
                Lval = quad(xfunction, lastx, xr)[0]/quad_min_max
                lastval, lastx = lastval+Lval, xr
                CDF.append(lastval)
            CDF = np.asarray(CDF)

        #Return CDF and w values
        return CDF, x_range

    def _inverse_cumulative_v(self, x, v_range, vc, CDF, x_range):

        '''
        Inverse sample from a provided velocity distribution.

        Parameters
        -----------
        x: numpy array
            The CDF values being inverse sampled at
        
        v_range: numpy array (length 2)
            Minimum and maximum velocity over which the CDF is determined

        vc: astropy quantity (velocity)
            MBHB circular velocity

        x_range, CDF: numpy arrays
            The list of v/vc values and the CDF evaluated at these points


        Returns
        -----------
        v: numpy array
            sampled velocities
        '''

        #Compute wmin and wmax
        wmin = float(v_range[0] / vc)
        wmax = float(v_range[1] / vc)

        #renormalize histogram
        f = interp1d(x_range,CDF)
        CDFlow = f(wmin)
        CDFhigh = f(wmax)
        newCDF = (CDF[(x_range>=wmin) & (x_range<=wmax)] - CDFlow) \
                / (CDFhigh-CDFlow)

        if(len(newCDF)<3):
            return np.zeros(len(x))

        newCDF[0] = 0
        newCDF[1] = 1

        #Sample from renormalized histogram
        random_number = np.atleast_1d(x)
        result = np.zeros(x.shape, dtype=np.float64)

        f = interp1d(newCDF, x_range[(x_range>=wmin) & (x_range<=wmax)])
        result = f(random_number)

        #Return sampled velocities
        return result

    def hardening_rate_HVSs(self, a):
        '''
        Binary hardening rate from star ejections from the MBHB slingshot 
        mechanism. Eq. 1

        Parameters
        ------------
        a:  astropy quantity (length, can be array)
            Current separation of the binary

        Returns
        ----------
        da/dt: astropy quantity (length/time)
            Binary hardening rate due to HVS ejections
        '''

        #Binary hardening rate from HVS ejection. See Eq. 9
        H = self.fit_H['A']*(1.+a/(self.ah*self.fit_H['a0'])) \
            **self.fit_H['gamma']
        
        return (-const.G * self.rho * H * a**2 / self.sigma).to(u.pc/u.s)
        
    def hardening_rate_GW(self, a):
        '''
        Hardening rate due to gravitational wave emission. Eq. 4

        Parameters
        -------------
        a: astropy quantity (length, can be array)
            current separation of the binary

        Returns
        ------------
        da/dt: astropy quantity (length/time)
            Hardening rate due to GW emission
        '''

        return ((-64./5.) * const.G**3 * const.c**(-5) * self.M_BH * self.Mc 
                    * (self.M_BH + self.Mc) / a**3).to(u.pc/u.s)

    def time_to_a(self, a):
        '''
        Time required from a binary to shrink from a0 to a via 
        HVS ejections and GW emission

        Parameters
        -----------
        a: astropy quantity (length, can be array)
           separations of the binary 

        Returns
        ----------------
        t: astropy quantity
            Time required for binary to shrink from a0 to each a
        '''

        #Defining dt/da
        dt = lambda a: (- 1 / ( self.hardening_rate_HVSs(a*u.pc) 
                    + self.hardening_rate_GW(a*u.pc) ) ).to('yr/pc').value

        t_range = []
        lastval, lastx = 0., self.a0.to(u.pc).value

        #Integrate to get total elapsed time
        for ai in tqdm(np.flip(a), desc='Integrating binary (HVS + GW)',total=len(a)):
                Lval = quad(dt, ai.to(u.pc).value, lastx)[0]
                lastval, lastx = lastval+Lval, ai.to(u.pc).value
                t_range.append(lastval)

        return t_range

    def time_to_a_hvs(self,a):
        '''
        Time required from a binary to shrink from a0 to a via 
        HVS ejections only

        Parameters
        -----------
        a: astropy quantity (length, can be array)
           separations of the binary 

        Returns
        ----------------
        t: astropy quantity
            Time required for binary to shrink from a0 to each a 
            via HVS ejections
        '''

        #Defining dt/da
        dt = lambda a: (- 1 / \
                ( self.hardening_rate_HVSs(a*u.pc) ) ).to('yr/pc').value

        t_range = []
        lastval, lastx = 0., self.a0.to(u.pc).value

        #Integrate to get total elapsed time
        for ai in tqdm(np.flip(a), desc='Integrating binary (HVS ejection only)',total=len(a)):

                Lval = quad(dt, ai.to(u.pc).value, lastx)[0]
                lastval, lastx = lastval+Lval, ai.to(u.pc).value
                t_range.append(lastval)

        return t_range

    def ejected_mass(self, a, a_range):
        """
        Mass ejected in stars while binary is in a particular range 
        of orbital separations

        Parameters
        ------------
        a: astropy quantity (length; can be array)
            Midpoint(s) of the orbital separation bin(s)
        a_range: astropy quantity (length, must be length of a+1)
            Edge(s) of the orbital separation bin(s)

        Returns
        ----------
        Delta M_ej: astropy quantity (mass)
            Total mass ejected while BH binary is in given separation bin(s)
        """

        #Restate separation bins as an inverse log fractions of ah
        log_aah_grid = np.log(self.ah/a_range)

        #spacings between ah/a bins
        delta_log_aah = np.array([abs(log_aah_grid[i].value \
            - log_aah_grid[i+1].value) for i in range(len(log_aah_grid)-1)])

        #Calculate dimensionless mass ejection rate J from fit parameters
        J = self.fit_J['A'] \
            * (a/(self.ah*self.fit_J['a0']))**self.fit_J['alpha'] \
            * (1. + (a/(self.ah*self.fit_J['a0']))**self.fit_J['beta']) \
            **self.fit_J['gamma']

        #Return total ejected mass
        return((J * (self.M_BH + self.Mc) * delta_log_aah).to('Msun'))

    def number_of_ejected_HVS(self, IMF, dMej):
        """
        Number of ejected HVSs for a given total ejection mass. Eq. 12

        Parameters
        ------------
        IMF: imf.imf instance
            Assumed initial mass function (see self.set_imf())
        dMej: astropy quantity (mass)
            Total ejected mass

        Returns
        ----------
        N_HVS: float
            Number of ejected HVSs
        """

        mbarimf = IMF.m_integrate(IMF.mmin,IMF.mmax)[0]*u.Msun

        return(dMej / mbarimf)

    def sampler(self):
        '''
        Samples from the ejection distribution to generate an ejection sample.
        The distribution mass and velocity is generated using a Montecarlo 
        approach (Rossi 2014). The functions _inverse_cumulative_mp, 
        _inverse_cumulative_a, _inverse_cumulative_q dictate the parameters of 
        the progenitor binary population. They are the inverse cumulative 
        distributions of the mass ratio, semi-major axis and primary mass 
        respectively. The velocity vector is assumed radial.

        Boundaries are imposed by default on the binary parameters:
            ::    0.1/mp<q<1, Rsun*(mp/Msun)<a<2000*Rsun, 0.1<mp<100

        Returns
        -------
            r0, phi0, theta0, v, phiv0, thetav0 : Quantity
                Initial phase space position in spherical coordinates
            m, tage, tflight
                Stellar mass of the HVS, age at observation and tflight 
                between ejection and observation

            n : int
                Size of the output ejection sample

        '''

        from .utils.mainsequence import t_MS
        from .utils.hurley_stellar_evolution import get_t_BAGB, get_t_MS
        from astropy.table import Table
        global vCDF, x_rangeCDF, time_range_glob, dMejs_glob, a_glob
        global Mctmp, MBHtmp, a0tmp, rhotmp, sigmatmp

        tzero = time.time()

        PI = np.pi

        #Hardening radius
        self.ah = (const.G * self.M_BH * self.q / (4*self.sigma**2)).to(u.pc) 

        #Orbital separations of BHB are separated into 1000 bins
        #Equally spaced in log space
        #Note, separations go from SMALL to LARGE, i.e. BACKWRDS in time
        # Some negative signs must therefore be used in some time calcs
        #If BMBH orbit has already been integrated and exists in globals, 
            #just grab it

        make_orbit = True
        if 'time_range_glob' in globals() and 'dMejs_glob' in globals() and 'a_glob' in globals():

            #Check to make sure relevant parameters of the MBHB 
            #and the GC haven't changed
            if self.Mc == Mctmp and self.M_BH == MBHtmp:
                if self.a0 == a0tmp and self.rho == rhotmp \
                    and self.sigma == sigmatmp:

                    time_range = time_range_glob.copy()
                    dMejs = dMejs_glob.copy()
                    a = a_glob.copy()
                    make_orbit = False

        #If MBHB orbit hasn't already been calculated, do that
        if make_orbit:

            Mctmp = self.Mc
            MBHtmp = self.M_BH
            a0tmp = self.a0
            rhotmp = self.rho
            sigmatmp = self.sigma

            a_range = np.logspace(np.log10(1e-7), np.log10(self.a0.value), 
                                1000) * u.pc

            #Binary separation where hardening rates due to HVS ejections and 
            #gravitational wave emission are equal       
            a_eq = a_range[np.argmin(np.abs(self.hardening_rate_HVSs(a_range) 
                            - self.hardening_rate_GW(a_range)))]
        
            #Midpoints of separation bins
            a = [(a_range[i].value + a_range[i-1].value)/2. \
                for i in range(1, len(a_range))] * u.pc


            #Get time required for binary to shrink from a0 to each a in grid
            time_range = self.time_to_a(a)
            time_range = np.flip(np.asarray(time_range)*1e-6*u.Myr)

            #Get time required for binary to shrink from a0 to each a in grid
            #if ONLY HVS ejections are hardening the binary.
            #Necessary to correct ejected mass
            time_rangehvs = self.time_to_a_hvs(a)
            time_rangehvs = np.flip(np.asarray(time_rangehvs)*1e-6*u.Myr)

            #Calculate mass ejection correction, Eq. 7
            corr2 = (-1 / np.diff(time_rangehvs,append=0.0*u.Myr)) \
                    / (-1 / np.diff(time_range,append=0.0*u.Myr))

            #Total stellar mass ejected when binary is in each 
            #separation bin, Eq. 7
            dMejs = self.ejected_mass(a, a_range)
            dMejs = dMejs * corr2 

            #Save the orbital evolution and mass ejection to file, if desired.
            if self.saveorbit:

                datalist = [time_range, a, dMejs.to('Msun'), 
                            self.hardening_rate_HVSs(a).to('pc/yr'), 
                            self.hardening_rate_GW(a).to('pc/yr')]
                namelist = ['time_range', 'a', 'dMej', 'rate_HVS', 'rate_GW']
                data_table = Table(data=datalist, names=namelist)
                
                data_table.write(self.orbitpath+'IMBH_inspiral_Mc_' + \
                                str(self.Mc.value)+'.fits', overwrite=True)


        #Save MBHB properties as global variables
        time_range_glob = time_range.copy()
        a_glob = a.copy()
        dMejs_glob = dMejs.copy()

        #Time at which separation is the hardening separation
        th = time_range[np.argmin(np.abs(a-self.ah))]   

        #Get total time from a=a0 to binary coalescence
        t_merged = max(time_range)
        print("Time binary merged =", t_merged)

        if self.tlb is not np.nan:
            t_lb = t_merged + self.tlb    

        elif self.current_a is not np.nan:
            t_lb = time_range[np.argmin(np.abs(a-self.current_a))]   

        else:
            raise ValueError(__mergerstateError__)

        if t_lb < self.tflightmax:
            print(a0Warning)

        #Need to do a timeline shift
        # The timeline now:
        # 0 ----------------------- t_merged
        # ^ binary starts at a=a0     ^ binary merges

        # The timeline we want to switch to
        # 0 ---------t0-------------------t_MW 
        # ^Big Bang  ^ binary is at a=a0   ^binary has been merging for t_lb
        #            |~~~~~~~<t_lb>~~~~~~~~|

        # 0 ---------t0---*----*-----*----t_MW 
        #                 ^    ^     ^ stars are ejected at t=t_ej
        # <t_ej> --> |~~~~|
        # <t_ej> --> |~~~~~~~~~|
        # <t_ej> --> |~~~~~~~~~~~~~~~|
        # Stars were ejected at time t_ej after t0

        # <t_flight> -->  |~~~~~~~~~~~~~~~|
        #      <t_flight> -->  |~~~~~~~~~~|             
        #           <t_flight>  -->  |~~~~|
        # ^ stars have been flying for t_flight

        #Shift the timeline. t=0 becomes t=th
        time_range += (self.T_MW - t_lb)

        # binary stars at separation a=a0 at t=t0
        # set t0
        t0 = time_range[-1]#-th#-(t_merged-t_lb)
        
        #set t_ej
        tejs = time_range - t0

        #set t_flight
        tflights = t_lb - tejs

        # Remove all stars which have a negative flight time.
        # Stars with a negative flight time have not yet been ejected, 
        # binary has not shrunk enough yet

        print("largest tflight:", max(tflights))
        idx_tflights_cut = (tflights.value>=0) & (tflights < self.tflightmax)
        tflights = tflights[idx_tflights_cut]
        a = a[idx_tflights_cut]

        #Set the initial mass function
        IMF = self.set_imf()

        dMejs = dMejs[idx_tflights_cut]

        #From the IMF and total ejected mass, get number of ejected stars
        Nhvs = self.number_of_ejected_HVS(IMF=IMF, dMej=dMejs)

        #Cut down sample size if required
        Nhvs /= self.cut_factor
        if(self.cut_factor>1):
            print(cutWarning)

        mps, seps, vcs, vs, tflight, tbefore, tage, mets = ( np.empty(0)
                                                        for i in range(8) )

        #Calculate the cumulative distribution of ejection velocities
        vCDF, x_rangeCDF = self.v_CDF_glob()

        for iN, N in enumerate(tqdm(Nhvs, desc="Generating HVSs")):

            #Sample from the IMF to get HVS stellar masses
            uniform_for_mp = np.random.uniform(0, 1, int(N))
            mp = self._inverse_cumulative_mp(IMF=IMF, x=uniform_for_mp)*u.Msun
            
            #Sample metallicities
            met = self.sample_met(len(mp))

            #Cut out stars outside of prescribed mass range
            idx = (mp > self.m_range[0]) & (mp < self.m_range[1])
            mp, met = mp[idx], met[idx]

            #Maximum age for HVS is the minimum between 
            #its lifetime to the tip of the AGB branch and the #age of the MW
            t_max = np.minimum(get_t_BAGB(mp, \
                                    met+np.log10(self.Zsun/0.02)), self.T_MW)

            # each star was already alive for some fraction e1 of its 
            # life before getting ejected. Eq. 14
            e1 = np.random.random(int(len(mp)))
            t_age_ej = e1 * t_max

            # age today
            tages = t_age_ej + tflights[iN]

            #Stars for whom t_age>t_max are stellar remnants today. 
            #Remove them from the sample.
            tages[tages>=t_max] = np.nan
            idx = (~np.isnan(tages))

            if idx.size != 0:

                # calculate circular velocity of the BH binary
                vc = self.calculate_vc(a[iN])

                tages, t_age_ej = tages[idx], t_age_ej[idx]
                mp, met = mp[idx], met[idx]

                uniform_for_v = np.random.uniform(0, 1, int(len(mp)))

                #Set maximum and minimum ejection velocity
                v_range = [self.vmin, 3*vc / (1 + self.q)]

                if(v_range[0]>=v_range[1]):
                    idx = np.full(len(mp),False)
                    v = np.zeros(len(mp))
                else:

                    #Use constructed CDF to sample from ejection 
                    #velocity distribution
                    v = self._inverse_cumulative_v(uniform_for_v, v_range, vc,
                                                     vCDF, x_rangeCDF)

                #Cut out velocities outside the prescribed range
                idx = (v > v_range[0]/vc) & (v < v_range[1]/vc)

                v = v[idx]
                mp, met = mp[idx], met[idx]
                tages, t_age_ej = tages[idx], t_age_ej[idx]
    
                if idx.size != 0:
                    mps = np.append(mps,mp.value)
                    vcs = np.append(vcs,vc.value*np.ones(mp.size))
                    vs = np.append(vs,v)
                    tflight = np.append(tflight,tflights[iN].value \
                                            * np.ones(mp.size))
                    seps = np.append(seps,a[iN].value * np.ones(mp.size))
                    tage = np.append(tage,tages.value)
                    tbefore = np.append(tbefore,t_age_ej.value)
                    mets = np.append(mets,met)

        m = np.array(mps) * u.Msun
        met = np.array(mets)
        vc = np.array(vcs) * u.km/u.s

        #Multiply v by vc to get actual velocities
        v = np.array(vs) * vc
        tflight = np.array(tflight) * u.Myr
        tage = np.array(tage) * u.Myr
        tbefore = np.array(tbefore) * u.Myr
        seps = np.array(seps)*u.pc

        n = m.size

        # Distance from GC at 3 pc
        r0 = np.ones(n) * self.centralr

        # Isotropic position unit vector in spherical coordinates
        phi0 = np.random.uniform(0, 2 * PI, n) * u.rad
        theta0 = np.arccos(np.random.uniform(-1, 1, n)) * u.rad

        phase = np.abs(np.pi*u.rad - phi0)

        idx1 = (m > 0.64*u.Msun) & (v>750*u.km/u.s)
        idx2 = (m > 0.85*u.Msun) | (phase > 2.85*u.rad)
        idx3 = (m > 1.3*u.Msun) | (phase > 1.35*u.rad)
        idx4 = v.value*tflight.value/m.value < 2.2e4*np.sqrt(m.value)
        idx = idx1*idx2*idx3*idx4

        #idx = (m.value > 0)

        r0, phi0, theta0, v, m, tage, tflight, vc = r0[idx], phi0[idx], theta0[idx], v[idx], m[idx], tage[idx], tflight[idx], vc[idx]
        tbefore, met = tbefore[idx], met[idx]
        seps = seps[idx]
        n = m.size

        # The velocity vector points radially.
        phiv0 = phi0.copy()
        thetav0 = theta0.copy()

        #Initialize radius, luminosity, effective temperature 
        #and evolutionary stage (today at at ejection) of the HVSs
        R, Lum, stage, stagebefore, T_eff = (np.empty(n) for i in range(5))

        #For each unique metallicity
        for z in tqdm(np.unique(met),desc='Evolving HVSs'):   

            #indices with ith metallicity, subset
            idx = np.where(met==z)[0]

            #Get star stage at moment of ejection
            star = self._evo_star(m[idx].value,tbefore[idx],met=z)
            stagebefore[idx] = star.stellar_type.as_astropy_quantity()

            #Get star stage at time of observation
            star = self._evo_star(m[idx].value,tage[idx],met=z)

            #Extract HVS effective temperature, radius, 
            #luminosity, evolutionary stage
            T_eff[idx] = star.temperature.as_astropy_quantity().to('K').value
            R[idx] = star.radius.as_astropy_quantity().to('Rsun').value 
            Lum[idx] = star.luminosity.as_astropy_quantity().to('Lsun').value 
            stage[idx] = star.stellar_type.as_astropy_quantity()        

        R = R*u.Rsun
        Lum = Lum*u.Lsun
        T_eff = T_eff*u.K

        returns = {'r0':r0, 'phi0':phi0, 'theta0':theta0, 'a_BHB':seps, 
                    'aah_BHB': seps/self.ah, 'v0':v, 'vc':vc, 'phiv0':phiv0, 
                    'thetav0':thetav0, 'm':m, 'tage':tage, 'tflight':tflight, 
                    'met':met, 'stage':stage, 'stagebefore':stagebefore, 
                    'Rad':R, 'T_eff':T_eff, 'Lum':Lum, 'size':len(r0), 
                    'Zsun': self.Zsun, 'Met': self.Met, 
                    'current_a': self.current_a, 'tlb': self.tlb, 
                    'cut_factor': self.cut_factor, 'Mc': self.Mc, 
                    'kappa': self.kappa, 'M_BH': self.M_BH, 'q': self.q, 
                    'a0': self.a0, 'rho': self.rho, 'sigma': self.sigma, 
                    'm_range': self.m_range, 'tflightmax': self.tflightmax, 
                    'saveorbit': self.saveorbit, 'orbitpath': self.orbitpath, 
                    'vmin': self.vmin }

        return returns                
