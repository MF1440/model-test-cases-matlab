import numpy as np
import json
from typing import NamedTuple
from earthConstants import *

class WalkerGroup(NamedTuple):

    inclination: float           # Наклонение орбиты [град.]
    satsPerPlane: int            # Число КА в каждой орбитальной плоскости группы
    planeCount: int              # Число орбитальных плоскостей в группе
    f: int                       # Фазовый сдвиг по аргументу широты между КА в соседних плоскостях [рад.]
    altitude: float              # Высота орбиты [км]
    maxRaan: float               # Максимум прямого восхождения восходящего узла (при распределении орбитальных плоскостей) [град.]
    startRaan: float             # Прямое восхождение восходящего узла для первой плоскости [град.]

    def getTotalSatCount(self):
        """
        Получение общего числа КА в группе
        :return: число КА в группе
        :rtype: int
        """
        return self.satsPerPlane * self.planeCount

    def getInitialElements(self):
        """
        Инициализация Кеплеровых элементов орбит группы
        
        :rtype: list
        :return: массив Кеплеровых элементов
        
        """
        startRaan   = np.deg2rad(self.startRaan)
        maxRaan     = np.deg2rad(self.maxRaan)
        inclination = np.deg2rad(self.inclination)
        altitude    = self.altitude * 1000
        satCount    = self.getTotalSatCount()

        raans = np.linspace(startRaan, startRaan + maxRaan, self.planeCount + 1)
        raans = raans[:-1] % (2 * np.pi)

        keplerianElems = np.zeros((satCount, 6))
        idx = 0

        for raanIdx, raan in enumerate(raans):
            for satIdx in range(self.satsPerPlane):
                sma = earthRadius + altitude
                aol = 2 * np.pi * (satIdx / self.satsPerPlane + self.f * raanIdx / satCount)

                keplerianElems[idx, :] = [sma, 0, 0, raan, inclination, aol]
                idx += 1

        return keplerianElems


class Constellation:

    def __init__(self, nameCode):
        self.totalSatCount = 0
        self.groups   = []
        self.keplerianElems = []
        self.stateEci = []
        self.loadFromConfig(nameCode)

    def loadFromConfig(self, nameCode):
        """
        Загрузка параметров из конфигурационного файла
        
        :param str nameCode: наименование группировки
        
        """
        f = open('../ConstellationsTest.json')
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
        """
        Получение перовначального состояния для каждой группы
        """
        self.keplerianElems = np.zeros((self.totalSatCount, 6))
        shift = 0

        for singleGroup in self.groups:
            ending = shift + singleGroup.getTotalSatCount()
            self.keplerianElems[shift:ending, :] = singleGroup.getInitialElements()
            shift = ending

    def propagateJ2(self, epochArray):
        """
        Пропагатор J2 для расчета положений спутников
        
        :param list epochArray: массив временных отметок на которых проводится расчет
        
        """
        self.stateEci = np.zeros((self.totalSatCount, 3, len(epochArray)))

        inclination = self.keplerianElems[:, 4]
        sma = self.keplerianElems[:, 0]
        raan0 = self.keplerianElems[:, 3]
        aol0 = self.keplerianElems[:, 5]

        raanPrecessionRate = -1.5 * (earthJ2 * np.sqrt(earthGM) * earthRadius**2) \
                           / (sma**(7/2)) * np.cos(inclination)

        draconicOmega      = np.sqrt(earthGM / sma**3) \
                           * (1 - 1.5 * earthJ2 * (earthRadius / sma)**2) \
                           * (1 - 4 * np.cos(inclination)**2)

        for epoch in epochArray:
            aol = aol0 + epoch * draconicOmega
            raanOmega = raan0 + epoch * raanPrecessionRate

            epochState = sma * [
                (np.cos(aol) * np.cos(raanOmega) - np.sin(aol) * np.cos(inclination) * np.sin(raanOmega)),
                (np.cos(aol) * np.sin(raanOmega) + np.sin(aol) * np.cos(inclination) * np.cos(raanOmega)),
                (np.sin(aol) * np.sin(inclination))]

            self.stateEci[:, :, epochArray.index(epoch)] = np.array(epochState).T
    
    
    def getGroupConnectivityMatrix(self, epochArray, maxDistance, groupIdx, allowAdjacentPlanes):
        """
        
        Получаем матрицу связности внутри одной группы спутников.

        :param list epochArray: временной отрезок, который был подан пропагатору J2.
        :param int minDistance: максимальное расстояние в [км], которое должно быть между спутниками, чтобы была возможность провести сеанс связи
        :param int groupIdx: номер группы (в ноль-индексации), внутри которой мы считаем матрицу связности. 
        :param bool allowAdjacentPlanes: флаг, который отвечает за то, учитываем ли мы соседние плоскости
        :return: две матрицы связности - с временными флагами и без
        :rtype: (np.array, np.array)
        """    
        
        maxDistanceMeters = maxDistance * 1000
        
        # Необходимо проверить, есть ли у нас ECI координаты, прежде чем считать расстояния
        if len(self.stateEci) == 0:
            print("Cперва запустите пропагатор!")
            return
            
        
        groupShift = 0 # Индекс первой размерности массива координат ECI, с которого начнутся спутники рассматриваемой группировки
        for thisGroupIdx in list(range(len(self.groups))):
            if thisGroupIdx == groupIdx:
                break
            
            groupShift += self.groups[thisGroupIdx].getTotalSatCount()
            
        thisGroupSatCount          = self.groups[groupIdx].getTotalSatCount() 
        thisGroupSatsPerPlaneCount = self.groups[groupIdx].satsPerPlane
        
        CsMatrixEpoches = np.zeros((thisGroupSatCount, thisGroupSatCount, len(epochArray[::10]), 2))
        
        print("Вычисляем...")
        
        for epochIdx, epoch in enumerate(epochArray[::10]):

            for satIdx1 in range(groupShift, groupShift + thisGroupSatCount - 1):
            
                adjPlaneDistIdx = None
                adjPlaneCsIdx  = None
            
                for satIdx2 in range(satIdx1 + 1, groupShift + thisGroupSatCount - 1):
                    
                    # Проверка, что оба спутника принадлежат одной орбитальной плоскости
                    satIdx1Plane = (satIdx1 - groupShift) // thisGroupSatsPerPlaneCount
                    satIdx2Plane = (satIdx2 - groupShift) // thisGroupSatsPerPlaneCount
                        
                    inSamePlane = satIdx1Plane == satIdx2Plane
                    inAdjPlane  = abs(satIdx1Plane - satIdx2Plane) == 1
                     
                    if allowAdjacentPlanes:
                        if not inAdjPlane and not inSamePlane:
                            continue
                    else:
                        if not inSamePlane:    
                            continue
                        
                    
                    r1 = self.stateEci[satIdx1, :, epoch]
                    r2 = self.stateEci[satIdx2, :, epoch]
                    
                    distance = np.sqrt(np.sum(np.square(r2 - r1)))
                    
                    if distance <= maxDistanceMeters:
                    
                        if inSamePlane:
                            # Проверка что у нас имеется 2 записи о сеансах связи
                            ThisSatCsArr   = np.argwhere(CsMatrixEpoches[:][satIdx2 - groupShift][:][0] == 1).tolist()
                            samePlaneCsArr = [tup for tup in ThisSatCsArr if 
                            tup[0] // thisGroupSatsPerPlaneCount == tup[1] // thisGroupSatsPerPlaneCount]
                          
                          
                            # Eсли имеется 2 сеанса, то подкорректируем сеансы так, чтобы в матрице оказались ближайшие по расстоянию
                            # Иначе просто добавим
                            if len(samePlaneCsArr) == 2:
                                samePlaneCsArr = [[arr[0]] + [satIdx2 - groupShift] + [arr[1]] for arr in samePlaneCsArr]
                                samePlaneCsArr = [arr + [0] for arr in samePlaneCsArr]
                                samePlaneDists = [arr[:-1] + [1] for arr in samePlaneCsArr]
                                
                                samePlaneCsArr     = list(map(tuple, samePlaneCsArr))
                                samePlaneDists   = list(map(tuple, samePlaneDists))
                                
                                firstDistance  = min(CsMatrixEpoches[samePlaneCsArr[0]], 
                                                     CsMatrixEpoches[samePlaneCsArr[1]])
                                                     
                                secondDistance = max(CsMatrixEpoches[samePlaneCsArr[0]], 
                                                     CsMatrixEpoches[samePlaneCsArr[1]])
                                
                                if CsMatrixEpoches[samePlaneDists[0]] < CsMatrixEpoches[samePlaneDists[1]]:
                                    firstDistanceIdx  = samePlaneDists[0]
                                    secondDistanceIdx = samePlaneDists[1]
                                    firstCsIdx  = samePlaneCsArr[0]
                                    secondCsIdx = samePlaneCsArr[1]
                                    
                                else:
                                    firstDistanceIdx  = samePlaneDists[1]
                                    secondDistanceIdx = samePlaneDists[0]
                                    firstCsIdx  = samePlaneCsArr[1]
                                    secondCsIdx = samePlaneCsArr[0]
                                    
                                secondCsIdxInv       = (secondCsIdx[1], secondCsIdx[0],
                                                        secondCsIdx[2], secondCsIdx[3])
                                    
                                secondDistanceIdxInv = (secondDistanceIdx[1], secondDistanceIdx[0],
                                                        secondDistanceIdx[2], secondDistanceIdx[3])
                                
                                firstCsIdxInv        = (firstCsIdx[1], firstCsIdx[0],
                                                        firstCsIdx[2], firstCsIdx[3])
                                    
                                firstDistanceIdxInv  = (firstDistanceIdx[1], firstDistanceIdx[0],
                                                        firstDistanceIdx[2], firstDistanceIdx[3])                                
                                
                                    
                                if distance > secondDistance:
                                    continue
                                    
                                elif distance <= secondDistance and distance > firstDistance:
                                    CsMatrixEpoches[secondDistanceIdx]    = 0
                                    CsMatrixEpoches[secondDistanceIdxInv] = 0
                                    CsMatrixEpoches[secondCsIdx]          = 0
                                    CsMatrixEpoches[secondCsIdxInv]       = 0
                                    
                                elif distance <= firstDistance:
                                    CsMatrixEpoches[firstDistanceIdx]    = 0
                                    CsMatrixEpoches[firstDistanceIdxInv] = 0
                                    CsMatrixEpoches[firstCsIdx]          = 0
                                    CsMatrixEpoches[firstCsIdxInv]       = 0
                                    
                        # Если спутники в соседних орб. плоскости, то так же проверим не взяли ли мы спутник поближе           
                        elif inAdjPlane:
                            if adjPlaneCsIdx is not None:
                                if distance < CsMatrixEpoches[adjPlaneDistIdx]:
                                    adjPlaneCsIdxInv    = (adjPlaneCsIdx[1], adjPlaneCsIdx[0], 
                                                            adjPlaneCsIdx[2], adjPlaneCsIdx[3])
                                                            
                                    adjPlaneDistIdxInv  = (adjPlaneDistIdx[1], adjPlaneDistIdx[0], 
                                                            adjPlaneDistIdx[2], adjPlaneDistIdx[3])
                                    
                                    CsMatrixEpoches[adjPlaneCsIdx]      = 0
                                    CsMatrixEpoches[adjPlaneCsIdxInv]   = 0
                                    CsMatrixEpoches[adjPlaneDistIdx]    = 0
                                    CsMatrixEpoches[adjPlaneDistIdxInv] = 0          
                                    
                                else:
                                    continue
                            #else:
                                
                            adjPlaneCsIdx   = (satIdx1 - groupShift, satIdx2 - groupShift, epochIdx, 0)
                            adjPlaneDistIdx = (satIdx1 - groupShift, satIdx2 - groupShift, epochIdx, 1)

                           
                        #Если сеансов меньше 2, то мы просто добавляем новый
                        CsMatrixEpoches[satIdx1 - groupShift][satIdx2 - groupShift][epochIdx][0] = 1
                        CsMatrixEpoches[satIdx2 - groupShift][satIdx1 - groupShift][epochIdx][0] = 1
                        CsMatrixEpoches[satIdx1 - groupShift][satIdx2 - groupShift][epochIdx][1] = distance
                        CsMatrixEpoches[satIdx2 - groupShift][satIdx1 - groupShift][epochIdx][1] = distance
                        
                            
        # Мы выдаем 2 матрицы - одна с временем для последующей маршрутизации с ее учетом
        # Вторая - просто говорит о наличии связи на протяжении всего времени
        CsMatrix = np.zeros((thisGroupSatCount, thisGroupSatCount))
        
        for satIdx1 in range(thisGroupSatCount):
            for satIdx2 in range(satIdx1, thisGroupSatCount):
                if np.any(CsMatrixEpoches[satIdx1, satIdx2, :, :] == 1):
                    CsMatrix[satIdx1][satIdx2] = 1
                    CsMatrix[satIdx2][satIdx1] = 1
        print("Готово!")
        return CsMatrix, CsMatrixEpoches            
                    
        
    def getGroupRouteWithMinHops(self, groupIdx, satIdxStart, satIdxFinish, CsMatrixEpoches):
        """
        Построение маршрута между двумя КА с помощью модификации алогоритма Дейкстры.

        :param groupIdx: номер группы, в которой мы строим маршрут
        :param int satIdxStart: id спутника, с которого мы хотим отправить данные.
        :param int satIdxFinish: id спутника, на который мы хотим отправить данные.
        :param np.ndarray CsMatrixEpoches: динамическая матрица связности, полученная в методе getGroupConnectivityMatrix().  
        :return: список вершин маршрута в порядке обхода
        :rtype: list or None
        """    

        import heapq
        
        satCount, _, epochCount, _ = CsMatrixEpoches.shape
    
        distances              = [float('inf')] * satCount
        distances[satIdxStart] = 0   
        prevNodes              = [None] * satCount
        priorityQueue          = [(0, satIdxStart)]
        
        while priorityQueue:
            curDistance, curNode = heapq.heappop(priorityQueue)
            
            if curNode == satIdxFinish:
                route = []
                
                while curNode is not None:
                    route.insert(0, curNode)
                    curNode = prevNodes[curNode]
                    
                return route
            
            for neighbor in range(satCount):
            
                if CsMatrixEpoches[curNode, neighbor, :, 0].any():
                
                    for epoch in range(epochCount):
                        if CsMatrixEpoches[curNode, neighbor, epoch, 0]:
                            if curDistance + CsMatrixEpoches[curNode, neighbor, epoch, 1] < distances[neighbor]:
                                distances[neighbor] = curDistance + CsMatrixEpoches[curNode, neighbor, epoch, 1]
                                prevNodes[neighbor] = curNode
                                heapq.heappush(priorityQueue, (distances[neighbor], neighbor))
        
        return None

