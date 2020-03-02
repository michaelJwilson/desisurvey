import  numpy              as      np
import  pylab              as      pl
import  matplotlib.pyplot  as      plt

from    astropy.table      import  Table


def X0(airmass):  
  ##  v1:         return  -0.01363498  + 0.15419213 * airmass ** 0.75 
  ##  v2 (7160):  
  return  -0.01474017 + 0.07491027 * airmass ** 0.75

def X1(airmass):
  ##  v1:         return  -7.68546710e-06 + 5.75539255e-05 * airmass ** 0.75
  ##  v2 (7160):  
  return  -6.84340500e-06 + 2.79721586e-05 * airmass ** 0.75

def Y0(moonfrac, airmass):
  ##  Dependence on moonfrac for a given airmass.                                                                                                                                                
  _X0 = X0(airmass)

  Y   = moonfrac - 0.1

  return  _X0 * np.exp(Y / 0.29) - _X0

def Y1(moonfrac, airmass):
  ##  Dependence on moonfrac for a given airmass.                                                                                                                                               
  _X1 = X1(airmass)
  Y   = moonfrac - 0.1

  return  _X1 * np.exp(Y / 0.29) - _X1

def Z0(moonalt, moonfrac, airmass):
  ##  Dependence on moonalt for fixed moonfrac and airmass.                                                                                                                                      
  ##  Note:  np.clip(Z0(_abs, res.x), 1., None)                                                                                                                                                  
  _Y0    = Y0(moonfrac, airmass)

  result =  1. + _Y0 * np.log(moonalt / 0.07)
  
  return  result
  
def Z1(moonalt, moonfrac, airmass):
  ##  Dependence on moonalt for fixed moonfrac and airmass.                                                                                                                                      
  ##  Note:  np.clip(Y1(_abs, res.x), 0.0, None)                                                                                                                                                 
  _Y1 = Y1(moonfrac, airmass)

  result = _Y1 * np.log(moonalt / 0.07)
  
  return  result
  
def moonmodel(moonsep, moonalt, moonfrac, airmass):
  ##  Dependence on separation for fixed moon alt., moon frac. and airmass.                                                                                                                      
  ##  Note:  np.clip(model(sample['MOONSEP'], res.x), 1., None)                                                                                                                                  
  if moonalt < 0.5:
    return  np.ones_like(moonsep)

  _Z0    = np.clip(Z0(moonalt, moonfrac, airmass), 1., None)
  _Z1    = np.clip(Z1(moonalt, moonfrac, airmass), 0., None)

  result = _Z0 + _Z1 * np.abs(moonsep - 96.) ** 1.86328519
  result = np.clip(result, 1., None)
  
  return  result
