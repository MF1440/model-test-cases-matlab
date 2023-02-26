clc
clear

tic
%создание объекта типа Constellation, инициализация параметрами группировки Starlink из конфига
constellation = Constellation('Starlink');

%вычисление элементов орбиты для всех КА в начальный момент
constellation.getInitialState();

% определение точек на оси времени, в которые будут проихзводиться расчёты
epochList = (0: 1000: 6000);

% расчёт положений всех КА в заданные моменты времени
constellation.propagateJ2(epochList);

% Координаты случайного КА (в инерциальных осях) после этого можно прочитать из constellation.state.eci
satIdx = ceil(constellation.totalSatCount * rand());
epochIdx = ceil(length(epochList) * rand());
disp(['Положение КА-' num2str(satIdx) ' на эпоху ' num2str(epochList(epochIdx)) ':']);
disp(constellation.state.eci(satIdx, :, epochIdx));

toc
