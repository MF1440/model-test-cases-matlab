clc
clear

tic
%создание объекта типа Constellation, инициализация параметрами группировки Stalink из конфига
constellation = Constellation('Starlink');

%вычисление элементов орбиты для всех КА в начальный момент
constellation.getInitialState(0);

% определение точек на оси времени, в которые будут проихзводиться расчёты
epochs = (0: 1000: 6000);

% расчёт положений всех КА в заданные моменты времени
constellation.propagateJ2(epochs);

% Координаты случайного КА (в инерциальных осях) после этого можно прочитать из constellation.state.eci
satIdx = ceil(constellation.totalSatCount * rand());
epochIdx = ceil(length(epochs) * rand());
disp(['Положение КА-' num2str(satIdx) ' на эпоху ' num2str(epochs(epochIdx)) ':']);
disp(constellation.state.eci(satIdx, :, epochIdx));

toc
