# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 22:15:57 2023

@author: Leraloger
"""

import numpy as np

# Класс для генерации сетки на сфере, полученной из спирали Фибоначчи.
class Sphere():
    
    # Инициальзация начальных значений полей.
    # radius -- радиус сферы;
    # pointCount -- количетсво точек в спирали.
    def __init__(self, radius: float, pointCount: int):
        self.radius = radius
        self.pointCount = pointCount
        
    # Генерация сетки по заданным параметрам.
    # pointList -- выходной вектор, содержащий 3-хмерные координаты точек
    # на сфере, рассчитанные относительно ее центра. 
    def generateSpherePointList(self) -> np.array:
        pointList = np.empty([self.pointCount, 3], dtype = np.float64)        
        # Угол золотого сечения.
        phi = np.pi - (np.sqrt(5.0) - 1.0)
        # Список идексов точек.
        idxList = np.arange(self.pointCount)
        # Расчет текущего смещение по оси спирали. Координата y.
        pointList[:, 1] = (1 - (idxList / float(self.pointCount - 1)) * 2)
        # Угол врашения текущей точки вокруг оси спирали.
        theta = phi * idxList
        # Расчет текущего радиуса спирали.
        tempRadius = np.sqrt(1 - pointList[:, 1] ** 2)
        # Вычисление координат спирали (x, z) в плоскости,
        # перпендикулрной оси спирали.
        pointList[:, 0] = tempRadius * np.cos(theta)
        pointList[:, 2] = tempRadius * np.sin(theta)   
        # Полученне точки рассчитаны для единичной сферы, поэтому необходимо
        # домножить на радиус сферы.
        pointList = self.radius * pointList
        
        return pointList