
PackageError = "One or more required external packages have not been installed. Please see requirements.txt"

PotentialError = 'Error: Propagation should be performed assuming a potential for which physical output has been turned on in galpy. This can be done by initializing a potential with astropy units or by using galpy.potential.turn_physical_on(). See https://docs.galpy.org/en/v1.9.0/getting_started.html#physunits for more information.'

propagatedWarning = 'Warning: sample appears not to have been propagated. Proceeding with mock photometry regardless'

import os
#try:
from astropy import units as u
import numpy as np
from tqdm import tqdm

#@propagate
def propagate(self, potential, dt=0.1*u.Myr, 
                  solarmotion=[-11.1, 12.24, 7.25],
                    zo=0.0208*u.kpc, orbit_path=None):

        '''
        Propagates the sample in the Galaxy forwards in time.

        Requires
        ----------
        potential : galpy potential instance
            Potential instance of the galpy library used to integrate the orbits

        Optional
        ----------
        dt : Quantity (time)
            Integration timestep. Defaults to 0.01 Myr
        solarmotion : length-3 list of floats
            UVW Solar peculiar velocity in km/s. 
            Galpy likes the U to be sign-flipped. Default is Schonrich+2010
        zo : Float or astropy distance quantity
             Offset of the Sun above or below the Galactic plane.
             Default is 20.8 pc (Bennett+Bovy 2019)
             If float, units are assumed to be kpc
        orbit_path : None or string
            If supplied, full equatorial and Galactocentric Cartesian orbits 
            are saved to orbit_path. Useful for debugging        
        '''

        import signal

        #try:
            #from galpydev.galpy.orbit import Orbit
            #from galpydev.galpy.util.conversion import get_physical
        from galpy.orbit import Orbit
        from galpy.util.conversion import get_physical
        from astropy.table import Table


        def handler(signum, frame):
            print('OOPS! A star took to long to integrate.'\
                  'Problematic star at index')

        from galpy.potential import evaluatePotentials

        #if u.get_physical_type('specific energy') \
        #!= u.get_physical_type(evaluatePotentials(potential, 8*u.kpc,0*u.kpc)):
        #    print(evaluatePotentials(potential,8*u.kpc,0*u.kpc))
        #    raise ValueError(PotentialError)
        
        self.solarmotion=solarmotion
        self.zo = zo

        #Read vo and ro from potential
        self.potential = potential
        self.vo = get_physical(potential)['vo']*u.km/u.s
        self.ro = get_physical(potential)['ro']*u.kpc

        # Integration time step
        self.dt = dt
        
        # Number of integration steps
        nsteps = np.ceil((self.tflight/self.dt).to('1').value).astype(int)

        # Impose 100 timesteps at minimum
        nsteps[nsteps<100] = 100

        # Initialize position in cylindrical coords
        rho   = self.r0 * np.sin(self.theta0)
        z     = self.r0 * np.cos(self.theta0)
        phi   = self.phi0
        phiv0 = self.phiv0

        #... and velocity
        vx = self.v0 * np.sin(self.thetav0) * np.cos(phiv0)
        vy = self.v0 * np.sin(self.thetav0) * np.sin(phiv0)
        vz = self.v0 * np.cos(self.thetav0)

        vR = vx*np.sin(phi+0.5 *np.pi*u.rad) + vy*np.sin(phi)
        vT = vx*np.cos(phi+0.5 *np.pi*u.rad) + vy*np.cos(phi)

        #Initialize a lot of stuff
        self.vx, self.vy, self.vz = (np.zeros(self.size)*u.km/u.s
                                                for i in range(3))
        self.x, self.y, self.z    = (np.zeros(self.size)*u.kpc 
                                                for i in range(3))
        self.ra, self.dec         = (np.zeros(self.size)*u.deg 
                                                for i in range(2))
        self.dist                 = np.zeros(self.size)*u.kpc          
        self.par                             = np.zeros(self.size)*u.mas
        self.pmra, self.pmdec     = (np.zeros(self.size)*u.mas/u.yr 
                                                for i in range(2))
        self.vlos                 = np.zeros(self.size)*u.km/u.s

        self.orbits                          = [None] * self.size

        self.b, self.l            = (np.zeros(self.size)*u.deg 
                                                for i in range(2))
 
        self.pmb, self.pml        = (np.zeros(self.size)*u.mas/u.yr 
                                                for i in range(2))   

        #Integration loop for the self.size orbits
        for i in tqdm(range(self.size),desc='Propagating...'):

            # Galpy will hang on propagating an orbit on rare occasion. 
            # Is only a problem if the flight time is long and/or the star 
            # makes multiple close pericentres to the SMBH. 
            # The signal alarm prevents this
            # signal.signal(signal.SIGALRM, handler)
            # signal.alarm(5)

            #Get timesteps
            #ts = np.linspace(0, 1, nsteps[i])*self.tflight[i]
            ts = np.linspace(-1, 0, nsteps[i])*self.tflight[i]

            #Initialize orbit using galactocentric cylindrical 
            #phase space info of stars
            #print([rho[i], vR[i], vT[i], z[i], vz[i], phi[i]])
            self.orbits[i] = Orbit(vxvv = [rho[i], vR[i], vT[i], z[i], vz[i], \
                                     phi[i]], solarmotion=self.solarmotion, \
                                     zo=zo, **get_physical(potential))

            self.orbits[i].integrate(ts, potential, method='dopr54_c')

            # Export the final position
            self.ra[i] = self.orbits[i].ra(ts, quantity=True)[-1]
            self.dec[i] = self.orbits[i].dec(ts, quantity=True)[-1]
            self.pmra[i] = self.orbits[i].pmra(ts, quantity=True)[-1]
            self.pmdec[i] = self.orbits[i].pmdec(ts, quantity=True)[-1]
            self.dist[i] = self.orbits[i].dist(ts, quantity=True)[-1]
            self.par[i] = u.mas / self.dist[i].to('kpc').value
            self.vlos[i] = self.orbits[i].vlos(ts, quantity=True)[-1]

            self.b[i] = self.orbits[i].bb(ts, quantity=True)[-1]
            self.l[i] = self.orbits[i].ll(ts, quantity=True)[-1]
            self.pmb[i] = self.orbits[i].pmbb(ts, quantity=True)[-1]
            self.pml[i] = self.orbits[i].pmll(ts, quantity=True)[-1]

            self.x[i] = self.orbits[i].x(ts, quantity=True)[-1]
            self.y[i] = self.orbits[i].y(ts, quantity=True)[-1]
            self.z[i] = self.orbits[i].z(ts, quantity=True)[-1]
            self.vx[i] = self.orbits[i].vx(ts, quantity=True)[-1]
            self.vy[i] = self.orbits[i].vy(ts, quantity=True)[-1]
            self.vz[i] = self.orbits[i].vz(ts, quantity=True)[-1]

            #Save the orbits of each HVS, if orbit_path is supplied
            if orbit_path is not None:
                #Only saves the orbits of the first 5e5 HVSs to prevent bloat

                if not os.path.exists(orbit_path):
                    raise SystemExit('Path '+orbit_path+' does not exist')

                if i<50000:
                    flightra    = self.orbits[i].ra(ts, quantity=True)
                    flightdec   = self.orbits[i].dec(ts, quantity=True)
                    flightdist  = self.orbits[i].dist(ts, quantity=True)
                    flightpmra  = self.orbits[i].pmra(ts, quantity=True)
                    flightpmdec = self.orbits[i].pmdec(ts, quantity=True)
                    flightvlos  = self.orbits[i].vlos(ts, quantity=True)
             
                    flightx  = self.orbits[i].x(ts, quantity=True)
                    flighty  = self.orbits[i].y(ts,quantity=True)
                    flightz  = self.orbits[i].z(ts,quantity=True)
                    flightvx = self.orbits[i].vx(ts,quantity=True)
                    flightvy = self.orbits[i].vy(ts,quantity=True)
                    flightvz = self.orbits[i].vz(ts,quantity=True)
                    flightL  = self.orbits[i].L(ts,quantity=True)
             
                    #Table of equatorial coordinates for the star with time
                    #datalist = [ts, flightra, flightdec, flightdist, 
                    #            flightpmra, flightpmdec, flightvlos]
                    #namelist = ['t', 'ra', 'dec', 'dist', 
                    #            'pm_ra', 'pm_dec', 'vlos']
                    #data_table = Table(data=datalist, names=namelist)

                    #Writes equatorial orbits to file. Each star gets own file
                    #data_table.write(orbit_path+'flight'+str(i)+'_ICRS.fits', 
                    #                 overwrite=True)

                    #Writes cartesian orbits to file. Each star gets own file
                    datalist=[ts, flightx, flighty, flightz, flightvx, flightvy, 
                                flightvz, flightL, flightra, flightdec, 
                                flightdist, flightpmra, flightpmdec, flightvlos]
                    namelist = ['t', 'x', 'y', 'z', 'v_x', 'v_y', 'v_z', 'L', \
                                'ra', 'dec', 'dist', 'pm_ra', 'pm_dec', 'vlos']

                    data_table = Table(data=datalist, names=namelist)
                    data_table.write(orbit_path+'flight'+str(i)+'_trev.fits', 
                                     overwrite=True)

        #signal.alarm(0)

        self.propagated = True

        #Get Galactocentric distance and velocity and Galactic escape velocity
        #as well as final azimuthal and polar coordinates
        #in Galactocentric sherical coordinates
        if(self.size>0):
            self.get_vesc(potential=potential)
            self.GCdist = np.sqrt(self.x**2. + self.y**2. \
                                    + self.z**2.).to(u.kpc)
            self.GCv = np.sqrt(self.vx**2. + self.vy**2. + \
                                    self.vz**2.).to(u.km/u.s)
            self.thetaf = np.arccos(self.z/self.GCdist)
            self.phif = np.arctan2(self.y,self.x)

        else:
            self.GCv = []*u.km/u.s 
            self.GCdist = []*u.kpc
            self.Vesc = []*u.km/u.s 
            self.thetaf = []*u.rad
            self.phif = []*u.rad

#@get_vesc
def get_vesc(self, potential):

        '''
        Returns the escape speed of a given potential 
        at each star's position in a propagated sample

        Requires
        ----------
        potential : galpy potential instance
            Potential instance of the galpy library used to integrate the orbits
        '''

        #try:
        from galpy.potential import evaluatePotentials
            #from galpydev.galpy.potential import evaluatePotentials


        if (not hasattr(self,'propagated') or not self.propagated):
            print(propagatedWarning)

        self.Vesc = np.zeros(self.size)*u.km/u.s
        
        R = np.sqrt(self.x**2 + self.y**2)
        z = self.z
        phi = np.arctan2(self.y,self.x)

        #potential.turn_physical_on()
        for i in range(self.size):
            #self.Vesc[i] = 220*np.sqrt(2*(evaluatePotentials(potential,
            #                            1e6*u.kpc,0*u.kpc) \
            #                - evaluatePotentials(potential,R[i],z[i])))*u.km/u.s
            self.Vesc[i] = np.sqrt(2*(evaluatePotentials(potential,1e6*u.kpc,
                            0*u.kpc,quantity=True) \
                            - evaluatePotentials(potential,R[i],z[i], 
                                                    phi=phi[i],quantity=True)))


def backprop(self, potential, dt=0.1*u.Myr, tint_max = 100*u.Myr, \
                  solarmotion=[-11.1, 12.24, 7.25], zo=0.0208*u.kpc, \
                orbit_path='./', equatorialSave = True, cartesianSave = True):

    '''
    Propagates the sample in the Galaxy backwards in time.

    Requires
    ----------
    potential : galpy potential instance
            Potentialused to integrate the orbits

    Optional
    ----------
    dt : astropy quantity (time)
            Integration timestep. Defaults to 0.1 Myr

    tint_max : astropy quantity (time)
            Maximum backwards integration time.

    solarmotion : length-3 list of floats
            UVW Solar peculiar velocity in km/s. 
            Galpy likes the U to be sign-flipped. Default is Schonrich+2010

    zo : Float or astropy distance quantity
             Offset of the Sun above or below the Galactic plane.
             Default is 20.8 pc (Bennett+Bovy 2019)
             If float, units are assumed to be kpc

    orbit_path : None or string
            Orbits are saved to orbit_path.

    equatorialSave : Boolean
            If True, backwards trajectories in the equatorial frame (ra, dec, 
            distance, proper motion, radial velocity) are saved to file.

    cartesianSave : Boolean
            If True, backwards trajectories in the Galactocentric Cartesian
            frame (x, y, z, vx, vy, vz) are saved to file.


    '''


    #Propagate a sample backwards in time. Probably obsolete

    #try:
    from galpy.orbit import Orbit
    #from galpy.util.coords import pmllpmbb_to_pmrapmdec, lb_to_radec
    #from galpy.util.coords import vrpmllpmbb_to_vxvyvz, lbd_to_XYZ
    from galpy.util.conversion import get_physical
    from astropy.table import Table
    import astropy.coordinates as coord


    #Creates path if it doesn't already exist
    if not os.path.exists(orbit_path):
            os.mkdir(orbit_path)

    self.solarmotion=solarmotion

    # Integration time step
    self.dt = dt

    # Maximum integration time
    #tint_max = 100.*u.Myr

    # Number of integration steps
    nsteps = int(np.ceil((tint_max/self.dt).to('1').value))

    # Initialize
    self.orbits = [None] * self.size

    #Integration loop for the n=self.size orbits
    #for i in range(self.size):
    for i in tqdm(range(self.size),desc='Backpropagating...'):
        #for i in range(1000):
        #print(self.name+' star index'+str(i))
        ts = np.linspace(0, 1, nsteps)*tint_max
        #ts = np.linspace(0, 1, nsteps[i])*self.tflight[i]

        #Initialize orbit instance using astrometry and motion of the Sun,
        #.flip() method reverses the orbit so we integrate backwards in time
        self.orbits[i] = Orbit(vxvv = [self.ra[i], self.dec[i], self.dist[i], \
                                self.pmra[i], self.pmdec[i], self.vlos[i]], \
                                    solarmotion=self.solarmotion, radec=True, \
                                    zo=zo, **get_physical(potential)).flip()

        self.orbits[i].integrate(ts, potential, method='dopr54_c')

        # Uncomment these and comment the rest of the lines in the for loop to return only final positions
        #self.dist[i], self.ll[i], self.bb[i], 
        #self.pmll[i], self.pmbb[i], self.vlos[i] = \
        #           self.orbits[i].dist(self.tflight[i], use_physical=True), \
        #           self.orbits[i].ll(self.tflight[i], use_physical=True), \
        #           self.orbits[i].bb(self.tflight[i], use_physical=True), \
        #           self.orbits[i].pmll(self.tflight[i], use_physical=True) , \
        #           self.orbits[i].pmbb(self.tflight[i], use_physical=True)  , \
        #           self.orbits[i].vlos(self.tflight[i], use_physical=True)

        self.backra, self.backdec, self.backdist, \
            self.backpmra, self.backpmdec, self.backvlos = \
                    self.orbits[i].ra(ts, quantity=True), \
                    self.orbits[i].dec(ts, quantity=True), \
                    self.orbits[i].dist(ts, quantity=True), \
                    self.orbits[i].pmra(ts, quantity=True), \
                    self.orbits[i].pmdec(ts, quantity=True), \
                    self.orbits[i].vlos(ts, quantity=True) 

        self.backx, self.backy, self.backz, self.backvx, self.backvy, \
                self.backvz = \
                    self.orbits[i].x(ts, quantity=True), \
                    self.orbits[i].y(ts, quantity=True), \
                    self.orbits[i].z(ts, quantity=True), \
                    self.orbits[i].vx(ts, quantity=True), \
                    self.orbits[i].vy(ts, quantity=True), \
                    self.orbits[i].vz(ts, quantity=True)

        #Assembles table of the equatorial/Cartesian coordinates 
        #for the particular star in each timestep
        datalist=[ts]
        namelist = ['t']        

        if cartesianSave:
            datalist.extend([self.backx, self.backy, self.backz, self.backvx, \
                    self.backvy, self.backvz])
            namelist.extend(['x', 'y', 'z', 'v_x', 'v_y', 'v_z'])

        if equatorialSave:
            datalist.extend([self.backra, self.backdec, self.backdist, \
                    self.backpmra, self.backpmdec, self.backvlos])
            namelist.extend(['ra', 'dec', 'dist', 'pm_ra', 'pm_dec', 'vlos'])

        data_table = Table(data=datalist, names=namelist)
        data_table.write(orbit_path+'/flight'+str(i)+'_backprop.fits', \
                            overwrite=True)

