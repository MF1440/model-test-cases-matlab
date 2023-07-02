# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 21:50:51 2023

@author: Leraloger
"""

import numpy as np

class CellConnector():
    
    # Инициальзация начальных значений полей.
    # gConst -- константа из стандарта 5G, ограничивающая 
    # максимальное количество идентификаторов;
    # walkerPointlist -- список всех КА в группировке;
    # connectionCountList -- список, содержащий данные о текущем количестве
    # соединений всех КА;
    # walkerIdList -- список, содержащий индексы соединенного КА для
    # каждой точки на Земле; 
    # connectionIdxList -- список, содержащий cellId для каждой точки.
    def __init__(self, spherePointCount: int, walkerPointList: np.array,
                 gConst: int, earthRadius: float):
        self.gConst              = gConst
        self.earthRadius         = earthRadius
        self.walkerPointList     = walkerPointList
        self.connectionCountList = np.zeros([np.size(walkerPointList, 0)])
        self.walkerIdList        = np.zeros(shape = spherePointCount, dtype = int)
        self.connectionIdxList   = np.zeros(shape = spherePointCount) - 1

    def checkQuadrant(self, walkerPoint: np.array, spherePoint: np.array) -> bool:
        return (bool)(np.prod(walkerPoint * spherePoint >= 0))
    
    # Функция проецирует позиции КА на заранее заданную сферу
    # радиусом self.earthRadius.
    def projectVectorList(self):
        vecLengthList = np.sqrt(np.sum(self.walkerPointList ** 2, axis = 1))
        tempVec = (self.earthRadius / vecLengthList)
        self.walkerPointList *= np.stack([tempVec, tempVec, tempVec], axis = 1)
  
    # Функция, в которой происходит построение связи между точками на Земле
    # и спутниками на основе расчета минимального растояния между точкой Земли
    # и подспутниковой точкой.
    def connectCellWithWalker(self, spherePointList: np.array):
        # Проекция положиний КА на сферу Земли.
        self.projectVectorList()
        # Цикл по всем точкам на Земле.
        for idx in range(np.size(spherePointList, 0)):
            # Здесь можно было бы добавить функционал уменьшения размрности
            # массива КА, по которому происходит расчет длин.
            # Один из таковых -- по коорлинатам или углам. Например, если угол
            # между КА и точкой больше заданного, то вектор не принимает
            # участия в дальнейшем.
            # Все реализации, пришедшие мне на ум, к сожалению имеют
            # сопоставимую с расчетом длины сложность, поэтому не дадут
            # выйгрыша в скорости. 
            # Разность векторов каждого КА и текущей точки.
            currentDistance = self.walkerPointList - spherePointList[idx]
            # Расчет длины каждого вектора.
            currentDistance = np.sum(currentDistance ** 2, axis = 1)
            # Нахождение минимальной длины.
            minDistance     = np.min(currentDistance)
            # Нахождение индекса КА с минимальным расстоянием.
            # Предпочтительнее сравнивать через двойное неравенство,
            # но я предположу, что Python скопировад данные без ошибки и 
            # подобная процедура зашита в операцию сравнения модуля NumPy.
            walkerIdx       = np.where(currentDistance == minDistance)[0][0]
            
            # Увеличение счетчика соединений для найденного КА.
            self.connectionCountList[walkerIdx] += 1
            # Присвоение cellId точке на Земле.
            self.connectionIdxList[idx]         = self.connectionCountList[walkerIdx]
            # Присвоение номера связанного КА.
            self.walkerIdList[idx]              = walkerIdx
        # Проверка на выход за пределы максимального количетсва соединений
        # для однго КА.
        badWalkerIdxList = np.where(self.connectionCountList > self.gConst)
        # Если есть хотя бы один КА.
        if len(badWalkerIdxList[0]):
            # Генерация исклчения о превышении количества свзей.
            exceptionString = "Превышено количество связей для спутников с номерами {0}"
            raise Exception(exceptionString.format(badWalkerIdxList[0]))
      
    # Возвращение привязки точек на поверхности к группировке КА.
    # result -- врзвращаемый массив со значениями;
    # Срез [:, 0] содержит cellId;
    # Срез [:, 1] содержит номер КА.
    def getResult(self) -> np.array:
        spherePointCount = np.size(self.walkerIdList, 0)
        result = np.zeros((spherePointCount, 2), dtype = int)
        result[:, 0] = self.connectionIdxList
        result[:, 1] = self.walkerIdList
        return result
        