function coveredSats = findVisibleSats(constellation, fileName, epochList, epochIdx, elevAngleMinDeg)
    
    epoch = epochList(epochIdx);
    
    stationList = jsondecode(fileread(fileName));
    
    coveredSats = {};
    
    satEciList = constellation.state.eci(:, :, epochIdx);
    
    stationEciList = zeros(length(stationList), 3);
    
    % Вычисление координат станций в системе ECI
    for stationIdx = 1: length(stationList)
        stationEciList(stationIdx,:) = calcEci(stationList(stationIdx).lat, stationList(stationIdx).lon, stationList(stationIdx).altitude, epoch);
    end
    
    for stationIdx = 1: length(stationEciList)

        visibleSatsPerStation = [];

        for satelliteIdx = 1: length(satEciList)

            stationSatVec = satEciList(satelliteIdx,:) - stationEciList(stationIdx,:);
            elevAngleRad = pi/2 - abs(calcAngle(stationSatVec, stationEciList(stationIdx,:)));
            
            % Проверка условия нахождения КА в зоне видимости шлюзовой станции
            if  rad2deg(elevAngleRad) >= elevAngleMinDeg
                visibleSatsPerStation(end+1) = satelliteIdx;
            end

        end

        coveredSats{end+1} = visibleSatsPerStation;
    end

