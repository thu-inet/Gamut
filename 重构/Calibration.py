from abc import ABCMeta, abstractclassmethod
from numpy import array
from numpy import polyfit, poly1d, log, exp

class Calibration():
    
    __meta_class__ = ABCMeta
    
    def __init__(self, *args, **kwargs):
        self.__info = kwargs.get("info")
    
    @abstractclassmethod
    def __call__(self, inp):
        pass
        
    def _get_info(self):
        print(f"Type: {self.__info['type']} Method:{self.__method}")
        print(f"{self.__info['input']:<8s}|{self.__info['input_unit']:<4s}|" + "|".join(f"{i:.3e}" for i in self.__data[:, 0]))
        print(f"{self.__info['output']:<8s}|{self.__info['output_unit']:<4s}|" + "|".join(f"{i:.3e}" for i in self.__data[:, 1]))
        
    def _data_check(self, data):
        if data is not None:
            if isinstance(data, array):
                pass
            elif isinstance(data, list):
                data = array(data)
            else:
                raise TypeError("Calibration data should be list or array-like!")
            
            if any(data<0):
                raise TypeError("Calibration data should be positive!")
            else:
                self.__data = data
    
    def _method_check(self, method):
        if method is not None:
            if method not in self.info['methods']:
                raise Warning(f"{self.info['type']} accepts interpolation methods including" + " ".join(self.info['methods']))

    @abstractclassmethod
    def _update(self):
        pass

class EnergyCalibration(Calibration):
    
    def __init__(self, data, method='linear', ):
        super().__init__(info={'type':'EnergyCalibration', 'input':'channel', 'output':'energy', 'input_unit':'chn', 'output_unit':'keV', 'methods':['linear', 'quadratic']})
        self._update(data, method)
        
    def _update(self, data=None, method=None):
        self._data_check(data)
        self._method_check(method)
        if self.__method == 'linear':
            deg = 1
        elif self.__method == 'quadratic':
            deg = 2
        self.__coefs = polyfit(self.__data[:, 0], self.__data[:, 1], deg)
        self.__func = poly1d(self.__coefs)
    
    def __call__(self, inp):
        return exp( self.__func(log(inp)))

class EfficiencyCalibration(Calibration):
    
    def __init__(self, data, method='log-log-linear'):
        super().__init__(info={'type':'EfficiencyCalibration', 'input':'energy', 'output':'efficiency', 'input_unit':'keV', 'output_unit':'%', 'methods':['log-log-linear']})
        self._update(data, method)
        
    def _update(self, data=None, method=None):
        self._data_check(data)
        self._method_check(method)
        if self.__method == 'log-log-linear':
            self.__coefs = polyfit(log(self.__data[:, 0]), log(self.__data[:, 1]), deg=1)
            self.__func = poly1d(self.__coefs)

    def __call__(self, inp):
        return exp( self.__func(log(inp)))

if __name__ == "__main__":
    
    ergcal = EnergyCalibration(data=[[10, 100], [20, 199]], method='linear')
    ergcal._print()
        
            