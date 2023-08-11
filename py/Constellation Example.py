from constellation import *
from random import randint
import numpy as np

# создание объекта типа Constellation, инициализация параметрами группировки Stalink из конфига
constellation = Constellation('Starlink')

# вычисление элементов орбиты для всех КА в начальный момент
constellation.getInitialState()

# определение точек на оси времени, в которые будут проихзводиться расчёты
epochs = list(range(1002))

# расчёт положений всех КА в заданные моменты времени
constellation.propagateJ2(epochs)

# Координаты случайного КА (в инерциальных осях) после этого можно прочитать из constellation.stateEci
satIdx = randint(0, constellation.totalSatCount - 1)
epochIdx = randint(0, len(epochs) - 1)
print("Положение КА-" + str(satIdx) + " на эпоху " + str(epochs[epochIdx]) + ":")
print(constellation.stateEci[satIdx, :, epochIdx])

csm, csme = constellation.getGroupConnectivityMatrix(epochs, 3000, 1, False)
print(csm)
print("Сеансов без учета соседних плоскостей:", np.count_nonzero(csm == 1)/2)

csm1, csme1 = constellation.getGroupConnectivityMatrix(epochs, 3000, 1, True)
print(csm1)
print("Сеансов с учетом соседних плоскостей:", np.count_nonzero(csm1 == 1)/2)

route = constellation.getGroupRouteWithMinHops(1, 0, 7, csme1)
print(route)