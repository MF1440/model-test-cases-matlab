# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 20:03:11 2023

@author: Leraloger
"""
from random import randint
import time
import numpy as np

from Constellation import Constellation
from Sphere import Sphere
from ColoredSurfaceVisualizer import ColoredSurfaceVisualizer
from CellConnector import CellConnector
from baseClasses import Const

# Класс для решения поставленной тестовой задачи
class Solver():
    
    # Инициализация значений радиуса сферы Земли и количетсва точек в сетке.
    # radius -- радиус Земли;
    # spherePointCount -- количество точек на сфере;
    def __init__(self, radius: float, spherePointCount: int):
        self.radius             = radius
        self.spherePointCount   = spherePointCount
        self.walkerPositionList = []
        self.constellation      = None
        
    # Инициализация расчетной модели спутниковой группировки, загрузка данных из 
    # файла и проведение моделирования.
    # constellationName -- имя группировки;
    # epochCount -- количество шагов по времени;
    def initializeConstellation(self, constellationName: str, epochCount: int):
        # создание объекта типа Constellation, инициализация параметрами
        # группировки Stalink из конфига
        self.constellation = Constellation(constellationName)
        # вычисление элементов орбиты для всех КА в начальный момент
        self.constellation.getInitialState()
        # определение точек на оси времени, в которые будут проихзводиться расчёты
        epochList = list(range(epochCount))
        # расчёт положений всех КА в заданные моменты времени
        self.constellation.propagateJ2(epochList)
        # Координаты всех КА (в инерциальных осях) в некоторую эпоху
        epochIdx = randint(0, len(epochList) - 1)
        # Копирование положения всех КА в эпоху epochIdx в поле класса.
        self.walkerPositionList = self.constellation.stateEciList[:, :, epochIdx]

    # Расчет точек на сфере, привязка их к спутникам, получение cellId для
    # каждой точки, визуализация результата и возвращение результата в виде np.array.
    # result -- выходной двумерный массив. По индексу точки получается вектор из двух
    # элементов. Первый элемент -- cellId точки, второй -- номер связанного КА.
    def getResultForTestCase(self) -> np.array:
        print("Начало генерации сетки")
        start = time.time()
        # Инициализация генератора точек и получение сетки.
        sphere = Sphere(radius = self.radius, pointCount = self.spherePointCount)
        spherePointList = sphere.generateSpherePointList()
        print("Начало расчета привязки")
        # Инициализация обхекта типа CellConnector.
        connector = CellConnector(spherePointCount = self.spherePointCount,
                                  walkerPointList = self.walkerPositionList,
                                  gConst = Const.fivegConst,
                                  earthRadius = Const.earthRadius)
        # Расчет привязки к КА и вычисление cellId.
        connector.connectCellWithWalker(spherePointList)
        result = connector.getResult()
        # Инициализция класса для рисования и отрисовка результата.
        visualizer = ColoredSurfaceVisualizer(self.constellation.totalSatCount)
        visualizer.showColoredSurface(spherePointList, colorMapIdx = result[:, 1])
        end = time.time()
        print("Расчетное время: {0}".format(end - start))
        return result

if __name__ == "__main__":
    # Значения параметров задачи.
    # Из-за выбранного алгоритма построения сетки задается точное количество узлов.
    spherePointCount = 64000
    # Количество эпох.
    epochCount = 1002
    # Название группровки КА.
    constellationName = "Starlink"
    # Создание обхекта-решателя Solver с передачей параметров радиуса
    # генерируемой сферы и количества точек на сфере.
    solver = Solver(Const.earthRadius, spherePointCount)
    # Запуск расчета.
    solver.initializeConstellation(constellationName, epochCount = epochCount)
    # Получение результата.
    result = solver.getResultForTestCase()
 
