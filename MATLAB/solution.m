clc
clear

%% Исходные данные
% диапазон для угла альфа (град.)
alfaListDegrees = 2:2:70;

% данные об абанентах
load testDataFile.mat

% Выделим абонентов с ненулевым запросом трафика

activeUserIdx = find(values);
activeUserList.Ecef = coordsEcef(activeUserIdx,:);      
activeUserList.Value = values(activeUserIdx);
activeUserList.originDist = sqrt( sum( (activeUserList.Ecef ).^2, 2 ));    % Удаленность абонентов от начала координат (можно было использовать R Земли, скорее всего)
activeUserCount = length(activeUserIdx);

clear activeUserIdx coordsEcef values

% Спутники
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

%% Ищем для каждого абонента ближайший спутник на каждой из орбит

tic
[minDistUserSat, satIdxList] = findNearestSat(activeUserList, satListEcef, constellation);
toc

%% Найдем соотношение углов

tic
satOriginDist = sqrt( sum( (satListEcef ).^2, 2 ));             % Удаленность спутников от начала координат
alfaMask = repmat(alfaListDegrees, activeUserCount,1);          % Альфа
logicalMask = zeros(activeUserCount, length(alfaListDegrees));  % Логическая маска

for groupIdx = 1:length(constellation.groupList)  
    
    a = satOriginDist(satIdxList(:,groupIdx));                  % Коротенькие переменные нужны были, чтобы не запутаться в формуле
    b = minDistUserSat(:,groupIdx);
    c = activeUserList.originDist;
    
    angleSatUser = acosd( (a.^2 + b.^2 - c.^2) ./ (2 .* a .* b) );
    
    temp = repmat(angleSatUser,1,length(alfaListDegrees));
    tempMask = alfaMask > temp;
    logicalMask = logicalMask | tempMask;
    
end
toc

%% Вычисляем трафик с учетом логической маски 

ValueRep= repmat(activeUserList.Value, 1, length(alfaListDegrees));
totalValue = sum(ValueRep .* logicalMask); 

%% Строим график
    
fig = figure;
fig.Position = [500,300,900,400];
plot(alfaListDegrees, totalValue, 'linewidth', 2); grid on
xlabel('\alpha\circ', 'fontsize', 14)
ylabel('TotalTraffic', 'fontsize', 14)
title('Зависимость удовлетворённого спроса от угла \alpha\circ', 'fontsize', 14, 'FontWeight', 'normal')

%% Функции 

function [minDistUserSat,satIdxList] = findNearestSat(activeUserList, satListEcef, constellation)
    
    % Расстояния от абонентов до ближайших спутников на каждой из орбит
    minDistUserSat = zeros(length(activeUserList.Value), length(constellation.groupList));
    
    % Индексы ближайших спутников на каждой из орбит
    satIdxList = zeros(length(activeUserList.Value),length(constellation.groupList));
    
    idx=1;
    
    for groupIdx = 1:length(constellation.groupList)
    
        totalSatCount = constellation.groupList{groupIdx}.totalSatCount;
        satEqualAlt = satListEcef(idx : idx + totalSatCount -1, :);
    
        [minDistUserSat(:,groupIdx),satIdxList(:,groupIdx)] = findNearestInGroup(activeUserList.Ecef, satEqualAlt);                                                                                                                                       
    
        idx = idx + totalSatCount;
    
    end
    
    
end

function[minDistUserSat,satInd] = findNearestInGroup(activeUserEcef, satListEcef)
    
    distUserSat = zeros(size(activeUserEcef,1),size(satListEcef,1));
    
    for sat = 1:size(satListEcef,1)   
       satPos = repmat(satListEcef(sat,:),size(activeUserEcef,1),1);
       distUserSat(:,sat) = sqrt( sum( (activeUserEcef - satPos).^2, 2 ));
    end
    
    % расстоние до ближайшего спутника и индекс этого спутника для каждого абонента    
    [minDistUserSat,satInd] = min(distUserSat,[],2);
    
end





