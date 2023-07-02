# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 22:25:29 2023

@author: Leraloger
"""

import random
import matplotlib.pyplot as plt
import numpy as np

# Класс для визуализации точек на сфере и их связи со спутниками.
# Визуализация происхожит с помощью функции scatter.
# Цвет точки обозначает принадлежность к определенному спутнику. 
# Цвета могут повторяться, т.к. задаются случайным образом.
# В условиях задачи не просили отображения, но для субъективной неточной оценки 
# распределения точек сферы по КА мне было проще добавить описываемый функционал.
class ColoredSurfaceVisualizer():
    
    # Инициальзация начальных значений полей.
    # colorCount -- количетсво различных цветов.
    def __init__(self, colorCount: int):
        self.colorCount = colorCount
        self.colorList = []
        
    # Создание набора строк, содержащих 16-битную кодировку RGB-цвета.
    def generateColorList(self):
        rgbColorLength = 6
        return np.array(["#" + ''.join([random.choice('0123456789ABCDEF')
                               for colorPartIdx in range(rgbColorLength)])
                for colorIdx in range(self.colorCount)])
    
    # Отображение поверхности в виде набора разноцветных точек.
    def showColoredSurface(self, pointList: np.array, colorMapIdx: np.array):
        # Заполнение пустого списка цветов.
        if not self.colorList:
            self.colorList = self.generateColorList()
        # Инициализаця массива строк с цветами для каждой точки. Нет проверки
        # на выход за границы массива, если количество цветов меньше чем 
        # количество спутников.
        colorMap = self.colorList[colorMapIdx]
        # создание объектов типа figure, axis и отрисовка точек на поверхности.
        fig, ax = plt.subplots(subplot_kw = {"projection": "3d"})
        ax.scatter(pointList[:, 0], pointList[:, 1], pointList[:, 2], c = colorMap)