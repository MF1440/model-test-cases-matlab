import numpy as np
from typing import NamedTuple


# Класс основных констант
class Const:
    earthRadius = 6378135          # Экваториальный радиус Земли [m]
    earthGM     = 3.986004415e+14  # Гравитационный параметр Земли [m3/s2]
    earthJ2     = 1.082626e-3      # Вторая зональная гармоника геопотенциала[?]
    fivegConst  = 1008             # Константа из стандарта 5G, ограничивающая
                                   # максимальное количество идентификаторов
    
# Класс, описывающий состояние одной спутниковой группировки.
# Возможно, "Walker" не самое подходящее имя для спутника, но, быть может, 
# есть определенная идея в таком обозначении.
class WalkerGroup(NamedTuple):
    inclination: float           # наклонение орбиты
    satsPerPlane: int            # число КА в каждой орбитальной плоскости
                                 # группы
    planeCount: int              # число орбитальных плоскостей в группе
    phaseShift: int              # фазовый сдвиг по аргументу широты между КА
                                 # в соседних плоскостях [?]
    altitude: float              # высота орбиты [?]
    maxRaan: float               # максимум прямого восхождения восходящего узла
                                 # (при распределении орбитальных плоскостей) [?]
    startRaan: float             # прямое восхождение восходящего узла для
                                 # первой плоскости [?]
    

    # Возвращает общее количество спутников в группировке
    def getTotalSatCount(self) -> int :
        return self.satsPerPlane * self.planeCount

    # Возвращает начальные положения всех спутников в групировке. 
    # В процессе переводит величины из одной размерности в другую.
    def getInitialElements(self) -> np.array :
        startRaan   = np.deg2rad(self.startRaan)
        maxRaan     = np.deg2rad(self.maxRaan)
        inclination = np.deg2rad(self.inclination)
        altitude    = self.altitude * 1000
        satCount    = self.getTotalSatCount()

        # Создаем массив из равномерно распределеного набора восходящих узлов от 
        # минимаотного до максимального значения
        raanList = np.linspace(startRaan, startRaan + maxRaan, self.planeCount + 1)
        raanList = raanList[:-1] % (2 * np.pi)

        elementList = np.zeros((satCount, 6))
        walkerIdx = 0

        for raanIdx, raan in enumerate(raanList):
            for satIdx in range(self.satsPerPlane):
                earthCenterDistance = Const.earthRadius + altitude
                #не уверен, что здесь именно latitudeAngle, но явно какой-то угол
                latitudeAngle = 2 * np.pi * (satIdx / self.satsPerPlane +
                                   self.phaseShift * raanIdx / satCount)

                elementList[walkerIdx, :] = [earthCenterDistance, 0, 0, raan,
                                    inclination, latitudeAngle]
                walkerIdx += 1

        return elementList

