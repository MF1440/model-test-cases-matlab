# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 20:03:11 2023

@author: Leraloger
"""

from Constellation import Constellation
from random import randint

# создание объекта типа Constellation, инициализация параметрами группировки Stalink из конфига
constellation = Constellation('Starlink')

# вычисление элементов орбиты для всех КА в начальный момент
constellation.getInitialState()

# определение точек на оси времени, в которые будут проихзводиться расчёты
epochList = list(range(1002))

# расчёт положений всех КА в заданные моменты времени
constellation.propagateJ2(epochList)

# Координаты случайного КА (в инерциальных осях) после этого можно прочитать из constellation.stateEci
satIdx   = randint(0, constellation.totalSatCount - 1)
epochIdx = randint(0, len(epochList) - 1)
print("Положение КА-" + str(satIdx) + " на эпоху " + str(epochList[epochIdx]) + ":")
print(constellation.stateEciList[satIdx, :, epochIdx])
print("Общее количество аппаратов: {0}".format(constellation.totalSatCount))