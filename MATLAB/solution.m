clc
clear

%% Исходные данные

alfaListDegrees = 2:2:70;        % диапазон для угла альфа (град.)

epoch = 1000;                    % произвольный момент времени
utc = [2023 2 25 12 0 0]; 

constellation_name = 'Starlink'; % Данные о спутниках

load testDataFile.mat;           % Данные об абонентах

%% Инициализация
% Выделим абонентов с ненулевым запросом трафика
users = getActiveUsers(coordsEcef, values);
clear coordsEcef values

% Выделим необходимые параметры группировки спутников
sats = getSatStruct(constellation_name,epoch,utc);

%%  Вычисляем логическую маску для абонентского спрса 
mask = zeros(users.count, length(alfaListDegrees)); 
alfaMatrixs = repmat(alfaListDegrees, users.count,1);

for groupIdx = 1:sats.groupCount        
    groupMask = findGroupMask(users, sats, groupIdx, alfaMatrixs);  
    mask = mask | groupMask;
end

%% Вычисляем трафик с учетом логической маски 
valuesMatrix= repmat(users.values, 1, length(alfaListDegrees));
totalValue = sum(valuesMatrix .* mask); 

%% Строим график
plotTotalValue(alfaListDegrees, totalValue)

%% Функции

function [users] = getActiveUsers(coordsEcef, values)
    % Функция выделяет абонентов с ненулевым запросом трафика
    
    activeUserIdx = find(values);
    users.ecef = coordsEcef(activeUserIdx,:);            % Координаты активных абонентов    
    users.values = values(activeUserIdx);                % Запрашиваемый трафик
    users.module = sqrt( sum( (users.ecef ).^2, 2 ));    % Удаленность абонентов от начала координат (можно было использовать R Земли)
    users.count = length(activeUserIdx);                 % Число активных абонентов

end

function [sats] = getSatStruct(constellation_name,epoch,utc)
    % Функция инициалирует группировку спутников, пересчитывает координаты
    % и возвращает структуру с необходимыми для дальнейших вычислений
    % данными
    
    
    % Инициализация группировки   
    constellation = Constellation(constellation_name);
    constellation.getInitialState();
    constellation.propagateJ2(epoch);
    satListEci = constellation.state.eci;

    % преобразуем координаты спутников из eci в ecef
    satListEcef = zeros(size(satListEci));

    parfor satIdx =  1:size(satListEci,1)
        satListEcef(satIdx,:) = eci2ecef(utc, satListEci(satIdx,:));
    end

    % Выделим группы спутников на разных высотах
    idx = 1;
    equalAltIdx = zeros(length(constellation.groupList),2);
    for groupIdx = 1:length(constellation.groupList)
        totalSatCount = constellation.groupList{groupIdx}.totalSatCount;
        equalAltIdx(groupIdx, 1) = idx;
        equalAltIdx(groupIdx, 2) = idx + totalSatCount -1;    
        idx = idx + totalSatCount;
    end
    %
    sats.ecef = satListEcef;
    sats.module = sqrt( sum( (satListEcef ).^2, 2 )); % Удаленность спутников от начала координат
    sats.equalAltIdx = equalAltIdx;
    sats.groupCount = length(constellation.groupList);

end

function[distanceMin, nearestSat] = findNearestSat(users, sats, groupIdx)
    %   Функция находит для всех абонентов ближайшие спутники в заданной
    %   группе
    %   distanceMin - расстояние абонент-ближайший спутник
    %   nearestSat - индекс ближайшего спутника    
    % Выделим спутники в текущей группе локальные переменные

    firstIdx = sats.equalAltIdx(groupIdx,1);
    lastIdx = sats.equalAltIdx(groupIdx,2);
    
    groupIdxRange = firstIdx:lastIdx;
    satsCount = length(groupIdxRange);
    
    localSats.ecef = sats.ecef(groupIdxRange,:);
    localSats.module = sats.module(groupIdxRange);
    distance = zeros(users.count,satsCount);
    
    for sat = 1:satsCount   
       satPos = repmat(localSats.ecef(sat,:),users.count,1);
       distance(:,sat) = sqrt( sum( (users.ecef - satPos).^2, 2 ));
    end
    
    % расстоние до ближайшего спутника и индекс этого спутника для каждого абонента    
    [distanceMin, localIdx] = min(distance,[],2);
    nearestSat = groupIdxRange(localIdx);
    
end

function [groupMask] = findGroupMask(users, sats, groupIdx, alfaMatrixs)
    
    % Функция проверяет, находится ли абонент в зоне радиовидимости
    % ближайшего спутника из заданной группы.
    % Возвращает логическую маску
    
    % Находим расстояние до ближайшего спутника и его индекс
    [distance, nearestSat] = findNearestSat(users, sats, groupIdx);
    
    a = sats.module(nearestSat);% Коротенькие переменные нужны, чтобы не запутаться в формуле
    b = distance;
    c = users.module;
    % Вычисляем угол
    angleSatUser = acosd( (a.^2 + b.^2 - c.^2) ./ (2 .* a .* b) );
    temp = repmat(angleSatUser,1,size(alfaMatrixs,2));
    % Сравниваем вычисленный угол с углом alfa (по всей сетке)
    groupMask = alfaMatrixs > temp;
    
end

function [] = plotTotalValue(alfaListDegrees, totalValue)

    fig = figure;
    fig.Position = [500,300,900,400];
    plot(alfaListDegrees, totalValue, 'linewidth', 2); grid on
    xlabel('\alpha\circ', 'fontsize', 14)
    ylabel('TotalTraffic', 'fontsize', 14)
    title('Зависимость удовлетворённого спроса от угла \alpha\circ', 'fontsize', 14, 'FontWeight', 'normal')
    
end
