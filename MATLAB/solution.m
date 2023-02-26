clc
clear

%% Исходные данные
% диапазон для угла альфа (град.)
alfaList = 2:2:70;

% данные об абанентах
load testDataFile.mat

%% Спутники

% создание объекта типа Constellation, инициализация параметрами группировки Starlink из конфига
constellation = Constellation('Starlink');
%вычисление элементов орбиты для всех КА в начальный момент
constellation.getInitialState();
% по условию момент времени произвольный
epoch = 1000;
% расчёт положений всех КА в заданный момент времени
constellation.propagateJ2(epoch);
% координаты спутников в инерциальной СК 
satListEci = constellation.state.eci;

%% преобразуем координаты спутников из eci в ecef

satListEcef = zeros(size(satListEci));

utc = [2023 2 25 12 0 0]; % произвольный момент времени

tic
parfor satIdx =  1:size(satListEci,1)
    satListEcef(satIdx,:) = eci2ecef(utc, satListEci(satIdx,:));
end
toc

%% график зависимости удовлетворённого спроса от угла alfa
tic
plotTotalTraffic(alfaList, satListEcef, coordsEcef, values, constellation);
toc

