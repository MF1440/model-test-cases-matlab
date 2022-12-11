clc
clear

tic
% Создание объекта типа Constellation, инициализация параметрами группировки Starlink из конфига
constellationConfigFileName = '../data/constellationsTest.json';
constellationName = 'Starlink';
constellation = Constellation(constellationConfigFileName, constellationName);

% Вычисление элементов орбиты для всех КА в начальный момент
constellation.initInitialState();

% Определение точек на оси времени, в которые будут проихзводиться расчёты
epochs = (0: 1000: 6000);

% Расчёт положений всех КА в заданные моменты времени
OrbitPropagators.propagateJ2(constellation, epochs)

% Координаты случайного КА (в инерциальных осях) после этого можно прочитать из constellation.state.eci
satIdx = ceil(constellation.totalSatCount * rand());
epochIdx = ceil(length(epochs) * rand());
disp(['Положение КА-' num2str(satIdx) ' на эпоху ' num2str(epochs(epochIdx)) ':']);
disp(constellation.state.eci(satIdx, :, epochIdx));

% Простой тест того, что рефакторинг кода не меняет результатов его работы
testEpochIdx = 7;
testSatIdx = 1630;
testSatCoords = 1.0e+06 * [-2.141680404794564   0.983991843448561   6.695407797436292];
testTreshold = 10^-15;
assert(norm(testSatCoords-constellation.state.eci(testSatIdx, :, testEpochIdx))/norm(testSatCoords) <= testTreshold);

toc
