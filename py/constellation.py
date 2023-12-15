# 1. В сам код внес только одно изменение, имя файла с данными начинается с маленькой буквы,
#    из-за этого код не работает под Linux.
# 2. По оформлению кода замечания не высказываю, на мой взгляд все нормально, насколько он 
#    удовлетворяет PEP8 не знаю, я давно не встречал чтобы его соблюдали. StyleGuide стараюсь
#    соблюдать в коде на C/C++.
# 3. Основное замечание по коду - то, что циклы можно убрать и заменить на векторные операции Numpy.

import numpy as np
import json
from typing import NamedTuple


class Parameters(object):
    pass


Const = Parameters()
Const.earthRadius = 6378135      # Экваториальный радиус Земли [m]
Const.earthGM = 3.986004415e+14  # Гравитационный параметр Земли [m3/s2]
Const.earthJ2 = 1.082626e-3      # Вторая зональная гармоника геопотенциала

group = Parameters()


class Walker(NamedTuple):
    inclination: float           # наклонение орбиты
    satsPerPlane: int            # число КА в каждой орбитальной плоскости группы
    planeCount: int              # число орбитальных плоскостей в группе
    f: int                       # фазовый сдвиг по аргументу широты между КА в соседних плоскостях
    altitude: float              # высота орбиты
    maxRaan: float               # максимум прямого восхождения восходящего узла (при распределении орбитальных плоскостей)
    startRaan: float             # прямое восхождение восходящего узла для первой плоскости


class WalkerGroup(Walker):

    def getTotalSatCount(self):
        return self.satsPerPlane * self.planeCount

    def getInitialElements(self):
        startRaan   = np.deg2rad(self.startRaan)
        maxRaan     = np.deg2rad(self.maxRaan)
        inclination = np.deg2rad(self.inclination)
        altitude    = self.altitude * 1000
        satCount    = self.getTotalSatCount()

        raans = np.linspace(startRaan, startRaan + maxRaan, self.planeCount + 1)
        raans = raans[:-1] % (2 * np.pi)

        # По второму измерению размерность 6, используются не все компоненты, часть  компонент заполнена постоянными значениями
        elements = np.zeros((satCount, 6))
        idx = 0

        # Можно убрать циклы
        for raanIdx, raan in enumerate(raans):
            for satIdx in range(self.satsPerPlane):
                sma = Const.earthRadius + altitude
                aol = 2 * np.pi * (satIdx / self.satsPerPlane + self.f * raanIdx / satCount)

                # Нули не используются
                elements[idx, :] = [sma, 0, 0, raan, inclination, aol]
                idx += 1

        return elements


class Constellation:

    def __init__(self, nameCode):
        self.totalSatCount = 0
        self.groups   = []
        self.elements = []
        self.stateEci = []
        self.loadFromConfig(nameCode)

    def loadFromConfig(self, nameCode):
        # В проекте файл называется с маленькой буквы
        f = open('constellationsTest.json')
        # Удобнее было бы все сделать через библиотеку Pandas
        jsonData = json.loads(f.read())

        # Можно заменить на enumerate
        for entryIdx in range(len(jsonData)):
            # nameCode.lower() - можно один раз заранее посчитать
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
        self.elements = np.zeros((self.totalSatCount, 6))
        shift = 0

        for singleGroup in self.groups:
            ending = shift + singleGroup.getTotalSatCount()
            self.elements[shift:ending, :] = singleGroup.getInitialElements()
            shift = ending

    def propagateJ2(self, epochs):
        self.stateEci = np.zeros((self.totalSatCount, 3, len(epochs)))

        inclination = self.elements[:, 4]
        sma = self.elements[:, 0]
        raan0 = self.elements[:, 3]
        aol0 = self.elements[:, 5]

        raanPrecessionRate = -1.5 * (Const.earthJ2 * np.sqrt(Const.earthGM) * Const.earthRadius**2) \
                           / (sma**(7/2)) * np.cos(inclination)

        draconicOmega      = np.sqrt(Const.earthGM / sma**3) \
                           * (1 - 1.5 * Const.earthJ2 * (Const.earthRadius / sma)**2) \
                           * (1 - 4 * np.cos(inclination)**2)

        # Можно сделать без цикла
        for epoch in epochs:
            aol = aol0 + epoch * draconicOmega
            raanOmega = raan0 + epoch * raanPrecessionRate

            # Это бросается в глаза больше всего.
            # Тригонометрические функции от одних и тех же аргументов вычисляются по несколько раз.
            # Можно рассмотреть вариант вычисления sin и cos на каждом шаге, используя формулы:
            # cos((n+1)x) = cos(nx)*cos(x) - sin(nx)*sin(x), аналоличную для sin и делать линейный переход
            # от n к n+1.
            epochState = sma * [
                (np.cos(aol) * np.cos(raanOmega) - np.sin(aol) * np.cos(inclination) * np.sin(raanOmega)),
                (np.cos(aol) * np.sin(raanOmega) + np.sin(aol) * np.cos(inclination) * np.cos(raanOmega)),
                (np.sin(aol) * np.sin(inclination))]

            self.stateEci[:, :, epochs.index(epoch)] = np.array(epochState).T
