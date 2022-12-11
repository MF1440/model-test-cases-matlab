%% Вам бы прикладывать к заданию методичку по небесной механике,
%% а то довольно сложно расшифровывать сокращения (гуглится все не очевидно).
%% Меня очень смущает название draconicOmega, но я решил его не трогать, вдруг это что-то внутреннее)
%% Сами формулы я не проверял, надеюсь этого не требовалось (в задании вроде нет такого)

clc
clear

tic
% Создание объекта типа Constellation, инициализация параметрами группировки Starlink из конфига
constellation = Constellation('Starlink');

% Вычисление элементов орбиты для всех КА в начальный момент
constellation.calcInitialState();

% Определение точек на оси времени, в которые будут производиться расчёты
epochs = (0: 1000: 6000);

% Расчёт положений всех КА в заданные моменты времени
constellation.propagateJ2(epochs);

% Координаты случайного КА (в инерциальных осях) после этого можно прочитать из constellation.state.eci
satIdx = ceil(constellation.totalSatCount * rand());
epochIdx = ceil(length(epochs) * rand());
disp(['Положение КА-' num2str(satIdx) ' на эпоху ' num2str(epochs(epochIdx)) ':']);
disp(constellation.state.xyzCoordinate(satIdx, :, epochIdx));

toc