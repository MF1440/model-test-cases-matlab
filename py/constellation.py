import json
import numpy as np
from baseClasses import WalkerGroup, Const

# Класс, описывающий состояние и моделирующий движение набора 
# спутниковых группировок на орбите Земли.  
class Constellation:

    # Инициальзация начальных значений полей и загрузка данных о группировке
    # из файла
    # nameCode -- название групировки
    def __init__(self, nameCode):
        self.totalSatCount = 0
        self.groupList     = []
        self.elementList   = []
        self.stateEciList  = []
        self.loadFromJsonConfig(nameCode)

    # Загрузка данных о спутниковых группировках в текуший объект
    # Расположение файла жестко прописано в коде функции.
    # nameCode -- название группировки
    def loadFromJsonConfig(self, nameCode: str):
        jsonFile = open('../ConstellationsTest.json')
        jsonData = json.loads(jsonFile.read())

        # цикл по всем эдементам массива группировок
        # поиск элемента массива (группировки) с нужным именем
        for entryIdx in range(len(jsonData)):
            # нашли группировку
            if (jsonData[entryIdx]['name']).lower() == nameCode.lower():
                constellationData = jsonData[entryIdx]
                print("Загружена группировка " + nameCode)

                for groupIdx in range(len(constellationData['Walkers'])):
                    self.groupList.append(WalkerGroup(*constellationData['Walkers'][groupIdx]))
                    self.totalSatCount += self.groupList[groupIdx].getTotalSatCount()

                jsonFile.close()
                return

        jsonFile.close()
        raise Exception('Группировка не найдена в файле')

    # Возвращает начальные параметры всех спутников во всех группировках,
    # хранаящихся в объекте.
    def getInitialState(self):
        self.elementList = np.zeros((self.totalSatCount, 6))
        groupIdxStart = 0

        for singleGroup in self.groupList:
            groupIdxEnd = groupIdxStart + singleGroup.getTotalSatCount()
            self.elementList[groupIdxStart:groupIdxEnd, :] = singleGroup.getInitialElements()
            groupIdxStart = groupIdxEnd

    # Модеоирование динамики движения спутников по орбите. 
    # epochList -- количество шагов по времени
    
    def propagateJ2(self, epochList):
        self.stateEciList = np.zeros((self.totalSatCount, 3, len(epochList)))

        # Не очень понятны названия переменых и их участие в алгоритме расчета
        sma         = self.elementList[:, 0]
        raan0       = self.elementList[:, 3]
        inclination = self.elementList[:, 4]
        aol0        = self.elementList[:, 5]

        raanPrecessionRate = -1.5 * (Const.earthJ2 * np.sqrt(Const.earthGM) *
                                     Const.earthRadius**2) \
                           / (sma**(7/2)) * np.cos(inclination)

        draconicOmega      = np.sqrt(Const.earthGM / sma**3) \
                           * (1 - 1.5 * Const.earthJ2 * (Const.earthRadius / sma)**2) \
                           * (1 - 4 * np.cos(inclination)**2)

        for epoch in epochList:
            aol       = aol0  + epoch * draconicOmega
            raanOmega = raan0 + epoch * raanPrecessionRate

            epochState = sma * [
                (np.cos(aol) * np.cos(raanOmega) - 
                 np.sin(aol) * np.cos(inclination) * np.sin(raanOmega)),
                (np.cos(aol) * np.sin(raanOmega) +
                 np.sin(aol) * np.cos(inclination) * np.cos(raanOmega)),
                (np.sin(aol) * np.sin(inclination))]

            self.stateEciList[:, :, epochList.index(epoch)] = np.array(epochState).T
