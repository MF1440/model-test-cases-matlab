clc
clear

tic
% Создание объекта типа Constellation, инициализация параметрами группировки Starlink из конфига
constellationConfigFileName = '../../data/constellationsTest.json';
constellationName = 'Starlink';
constellation = Constellation(constellationConfigFileName, constellationName);

% Вычисление элементов орбиты для всех КА в начальный момент
constellation.initInitialState();

% Определение точек на оси времени, в которые будут производиться расчёты
firstEpoch = 0;
lastEpoch = 6000;
epochStep = 1000;
epochs = (firstEpoch: epochStep: lastEpoch);

% Расчёт положений всех КА в заданные моменты времени
OrbitPropagators.propagateJ2(constellation, epochs)

% Инициализация из файла положений наземных станций и значений трафика на этих станциях
stationsConfigFileName = '../../data/testDataFile.mat';
newData = load('-mat', stationsConfigFileName);

% Генерация тестовой выборки станционных данных заданного размера
testStationsCount = 100000;
stationsCount = numel(newData.values);
idx = randperm(stationsCount);

stationsEcefCoords = newData.coordsEcef(idx(1:testStationsCount), :);
stationsTrafficArray = newData.values(idx(1:testStationsCount));

testEpochIdx = 1;
testEpoch = epochs(testEpochIdx);

% Вычисление сервисных зон
serviceZones  = calcServiceZones(stationsEcefCoords, stationsTrafficArray, ...
    constellation.state.eci(:,:,testEpochIdx), testEpoch);

toc

isVisualizationNeeded = true; 

if isVisualizationNeeded
    satVisIdx = 1;
    visualizeTestResults(stationsEcefCoords, constellation.state.eci(:,:,testEpochIdx), ...
        serviceZones.seriveSatIndexes, testEpoch, satVisIdx)
end


function visualizeTestResults(stationsEcefCoords, satsEciCoords, serviceSatIndexes, epoch, satVisIdx)
% Визуалзация подспутниковых точек и зоны обслуживания спутника с индексом satVisIdx

    figure(1)

    satsEcefCoords = Utils.eci2ecef(satsEciCoords, epoch);
    satsEcefCoords = Constants.AstroConstants.earthRadius * normalize(satsEcefCoords, 2, 'norm');
    
    scatter3(satsEcefCoords(:,1), satsEcefCoords(:,2), satsEcefCoords(:,3), 10*ones(length(satsEcefCoords(:,1)),1), '.')
    hold on
    
    idx = (serviceSatIndexes == satVisIdx);
    scatter3(satsEcefCoords(satVisIdx,1), satsEcefCoords(satVisIdx,2), satsEcefCoords(satVisIdx,3), 'r+')
    scatter3(stationsEcefCoords(idx,1), stationsEcefCoords(idx,2), stationsEcefCoords(idx,3), 50, 'k.')

    hold off

end


