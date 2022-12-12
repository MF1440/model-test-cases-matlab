%% Я считал, что координаты спутников уже посчитаны в классе Constellation.
%% Расчет на 7 точек времени у меня занял примерно 660 секунд

clc
clear

tic

% Создание объекта типа Constellation, инициализация параметрами группировки Starlink из конфига
constellation = Constellation('Starlink');

% Вычисление элементов орбиты для всех КА в начальный момент
constellation.calcInitialState();

% Определение точек на оси времени, в которые будут производиться расчёты
%epochs = (0: 1000: 2000);
epochs = 0;

% Расчёт положений всех КА в заданные моменты времени
constellation.propagateJ2(epochs);

% Создание объекта типа Traffic
traffic = Traffic(constellation.state.ecefPosition,'testDataFile.mat');

% Распределение пользователей по спутникам
traffic.calcSatTraffic();

toc