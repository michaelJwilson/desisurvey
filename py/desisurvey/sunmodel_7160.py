import  numpy                                     as      np
import  pylab                                     as      pl
import  matplotlib.pyplot                         as      plt

from    astropy.table                             import  Table
from    desisurvey.skyfactor.py.sun.indices       import  index_model


def Y0(airmass, params):
  ##  [6.9489374,  -3.96722344]                                                                                                                                                                                                       
  return  params[0] + params[1] * airmass

def Y1(airmass, params):
  ##  [4.07622358, -0.09630329]                                                                                                                                                                                                       
  return  params[0] + params[1] * airmass

def Z(sunalt, params):
  ##  Dependence on sunalt. for a given airmass.                                                                                                                                                                                     
  return  np.exp(5. + sunalt / params[1]) - params[0]

def sunmodel(sunsep, airmass, sunalt):
  #  Airmass dependence of sunalt shape;  See sunalts.py;                                                                                                                                                                            
  if sunalt < -18.:
    return  np.ones_like(sunsep)

  if sunalt > -14:
    ##  Invalid above -14 deg. 
    return  1.e99

  _Y0    =  Y0(airmass, np.array([6.9489374,  -3.96722344]))
  _Y1    =  Y1(airmass, np.array([4.0762236, -0.09630329]))
  
  # One parameter model for dependence on sun alt. for given airmass, Z.                                                                                                                                                             
  _Z     =  Z(sunalt, np.array([_Y0, _Y1]))

  index  =  index_model(airmass, np.array([2.00688895, -0.01108773, -0.00291262]))

  # Dependence on sunsep for fixed sunalt and airmass.                                                                                                                                                                               
  result =  np.ones_like(sunsep) + _Z + 3.79912194e-04 * np.abs(sunsep - 1.2356e+02) ** index
  result =  np.clip(result, 1., None) 

  return  result
