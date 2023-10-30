clc
clear

tic
%создание объекта типа Constellation, инициализация параметрами группировки Stalink из конфига
constellation = Constellation('biWalkerGlobal');

%вычисление элементов орбиты для всех КА в начальный момент
constellation.getInitialState();

% определение точек на оси времени, в которые будут проихзводиться расчёты
epochs = (0: 1000: 6000);

% расчёт положений всех КА в заданные моменты времени
constellation.propagateJ2(epochs);

% Координаты случайного КА (в инерциальных осях) после этого можно прочитать из constellation.state.eci
satIdx = randi(constellation.totalSatCount);
epochIdx = ceil(length(epochs) * rand());
disp(['Положение КА-' num2str(satIdx) ' на эпоху ' num2str(epochs(epochIdx)) ' в ECI:']);
disp(constellation.state.eci(satIdx, :, epochIdx));

% Координаты случайного КА (в осях, связанных с вращающейся Землёй) можно прочитать из constellation.state.ecef
disp(['Положение КА-' num2str(satIdx) ' на эпоху ' num2str(epochs(epochIdx)) ' в ECEF:']);
disp(constellation.state.ecef(satIdx, :, epochIdx));
toc
