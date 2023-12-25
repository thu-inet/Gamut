class Peak():
    def __init__(self,centroid,peak_left_boundary,peak_right_boundary,gamma,energy=0):
        self.centroid=centroid
        self.peak_left_boundary=peak_left_boundary
        self.peak_right_boundary=peak_right_boundary

class Gamma():

    def __init__(self,nuclide,energy,half_life,intensity,NFflags,PFflags):
        self.nuclide=nuclide
        self.energy=energy
        self.half_life=half_life
        self.intensity=intensity
        self.NFflags=NFflags
        self.PFflags=PFflags

class Nuclide():

    def __init__(self,nuclide,fyield,nucid,gammas=[]):
        self.nuclide=nuclide
        self.fyield=fyield
        self.nucid=nucid
        self.gammas=gammas