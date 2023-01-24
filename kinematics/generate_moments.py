import astropy.units as u
from astropy.utils import data
from spectral_cube import SpectralCube

#cube = SpectralCube.read('HCN32_high_res_sc3_.fits') #HCN 3-2
cube = SpectralCube.read('ngc253.B6b.cent.12m7m.230300.contsub.cv0_4.cube_sc5.0_.fits') #CO 2-1 @230.538GHz (230.338 in NGC253)

HCN32_cube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=230.338*u.GHz)  

HCN32_subcube =  HCN32_cube.spectral_slab(-150*u.km/u.s,150*u.km/u.s) 

moment_0 = HCN32_subcube.moment(order=0)
moment_1 = HCN32_subcube.moment(order=1)

moment_0.write('moment0_CO32.fits', overwrite=True) 
moment_1.write('moment1_CO32.fits', overwrite=True) 
