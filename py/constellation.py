#!/usr/bin/env python
# coding: utf-8

# In[291]:


import numpy as np
import json
from typing import NamedTuple


# In[292]:


class Parameters(object):
    pass

Const = Parameters()
Const.earthRadius = 6378135;           # Экваториальный радиус Земли [km]
Const.earthGM = 3.986004415e+14;       # Гравитационный параметр Земли [m3/s2]
Const.earthJ2 = 1.082626e-3;           # First zonal harmonic coefficient in the expansion of the Earth's gravity field

group = Parameters()


# In[293]:


class Walker(NamedTuple):
    inclination:  float
    satsPerPlane: int
    planeCount:   int 
    f:            int
    altitude:     float
    maxRaan:      float
    startRaan:    float
                
class WalkerGroup(Walker):

    def getTotalSatCount(self): 
        return self.satsPerPlane * self.planeCount;
    
    def getInitialElements(self):
        startRaan = np.deg2rad(self.startRaan)
        maxRaan = np.deg2rad(self.maxRaan)
        inclination = np.deg2rad(self.inclination)
        altitude = self.altitude * 1000
        
        raans = np.linspace(startRaan, startRaan + maxRaan, self.planeCount + 1)
        raans = raans[0:-1] % (2 * np.pi)

        elements = np.zeros((self.getTotalSatCount(), 6))
        idx = 0
        raanIdx = 0
        for raan in raans:
            for satIdx in range(self.satsPerPlane):
                sma = Const.earthRadius + altitude
                aol = 2 * np.pi / self.satsPerPlane * (satIdx + 1) + 2 * np.pi / self.getTotalSatCount() * self.f * raanIdx

                elements[idx, :] = [sma, 0, 0, raan, inclination, aol]
                idx += 1
                    
            raanIdx += 1
            
        return elements


# In[294]:


class Constellation:
    
    def __init__(self, nameCode):
        self.totalSatCount = 0
        self.groups = []
        self.elements = []
        self.loadFromConfig(nameCode)
        
    def loadFromConfig(self, nameCode):
        f = open('constellationsTest.json')
        jsonData = json.loads(f.read())
        
        for entryIdx in range(len(jsonData)):
            if (jsonData[entryIdx]['name']).lower() == nameCode.lower():
                print("Загружена группировка " + nameCode)
                constellationData = jsonData[entryIdx]
                
                for groupIdx in range(len(constellationData['Walkers'])):
                    self.groups.append(WalkerGroup(*constellationData['Walkers'][groupIdx]))
                    self.totalSatCount += self.groups[groupIdx].getTotalSatCount()
                
                f.close()
                return 
            
        f.close()
        raise Exception('Группировка не найдена в файле')
        
    def getInitialState(self):
        self.elements = np.zeros((self.totalSatCount, 6));
        shift = 0;

        for group in self.groups:
            ending = shift + group.getTotalSatCount()
            self.elements[shift:ending, :] = group.getInitialElements()
            shift = ending
            

    def propagateJ2(self, epochs):
        self.stateEci = np.zeros((self.totalSatCount, 3, len(epochs)))

        inclination = self.elements[:, 4]
        sma         = self.elements[:, 0]
        Omega0      = np.sqrt(Const.earthGM / sma**3)
        aol0        = self.elements[:, 5]

        raanPrecessionRate = -1.5 * (Const.earthJ2 * np.sqrt(Const.earthGM) * Const.earthRadius**2) / (sma**(7/2)) * np.cos(inclination)
        draconicOmega      = np.sqrt(Const.earthGM / sma**3) * (1 - 1.5 * Const.earthJ2 * (Const.earthRadius / sma)**2) * (1 - 4 * np.cos(inclination)**2)

        breakpoint
        
        for epoch in epochs:
            aol = aol0 + epoch * draconicOmega
            Omega = Omega0 + epoch * raanPrecessionRate

            epochState = sma * [(np.cos(aol) * np.cos(Omega) - np.sin(aol) * np.cos(inclination) * np.sin(Omega)), 
                                (np.cos(aol) * np.sin(Omega) + np.sin(aol) * np.cos(inclination) * np.cos(Omega)), 
                                (np.sin(aol) * np.sin(inclination))]
            self.stateEci[:, :, epochs.index(epoch)]  = np.array(epochState).T 


# In[ ]:




