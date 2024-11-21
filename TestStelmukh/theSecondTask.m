clc
clear

%создание объекта типа Constellation, инициализация параметрами группировки Stalink из конфига
constellation = Constellation('biWalkerGlobal');

%вычисление элементов орбиты для всех КА в начальный момент
constellation.getInitialState();

% определение точек на оси времени, в которые будут производиться расчёты
epochs = 0:900:86400; 

% расчёт положений всех КА в заданные моменты времени
constellation.propagateJ2(epochs);

% граничные условия на угол раствора конуса
minAlpha = 0;
maxAlpha = 90;
accuracyAlpha = 3; 

% алгоритм бинарного поиска
while (maxAlpha - minAlpha) > accuracyAlpha
    targetAlpha = (minAlpha + maxAlpha) / 2;
    if isGlobalCoverage(constellation.state.ecef, targetAlpha)
        maxAlpha = targetAlpha; 
    else
        minAlpha = targetAlpha; 
    end
end

disp(targetAlpha);