function serviceZones = calcServiceZones(stationsEcef, stationsTrafficArray, satsEci, epoch)
% Функция вычилсяет распределение наземных станций по зонам обслуживания спутников
% и производит вычисление суммарного трафика в каждой зоне обслуживания
% Положение спутников в момент времени epoch определяется массивом satsEci 
% Положение наземных станций определяется массивом stationsEcef, траффик 
% через данную станцию - массив stationsTrafficArray 
% На выходе функция возвращает структуру serviceZones
% serviceZones.serviceSatIndexes - номер спутника к обслуживающей зоне которого относится данная станция
% serviceZones.satsTrafficArray - суммарный трафиик в данной зоне

    serviceZones.seriveSatIndexes = zeros(size(stationsEcef, 1), 1);
    serviceZones.satsTrafficArray = zeros(size(stationsEcef, 1), 1);

    satsEcef = Utils.eci2ecef(satsEci, epoch);

    [satsLon,     satsLat,     ~] = cart2sph(satsEcef(:,1),     satsEcef(:,2),     satsEcef(:,3));    
    [stationsLon, stationsLat, ~] = cart2sph(stationsEcef(:,1), stationsEcef(:,2), stationsEcef(:,3));
    
    satsSinLat = sin(satsLat);
    satsCosLat = cos(satsLat);
    
    stationsCount = length(stationsTrafficArray);

    for stationIdx = 1:stationsCount
        % Вычисление центрального угла между даной станцией и спутником
        angle = (satsSinLat * sin(stationsLat(stationIdx)) + ...
                (satsCosLat * cos(stationsLat(stationIdx))) .* cos(satsLon - stationsLon(stationIdx)));
        [~, indMin] = max(angle);
        % Зона обслуживания станции определяется на основании минимума расстояния от станции до подспутниковой точки
        serviceZones.serviceSatIndexes(stationIdx) = indMin;
        serviceZones.satsTrafficArray(indMin) = serviceZones.satsTrafficArray(indMin) + stationsTrafficArray(stationIdx);
    end

end

