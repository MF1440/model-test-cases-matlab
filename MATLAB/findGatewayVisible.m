function visiblesArray = findGatewayVisible(constellation, epochArray, epochIdx)
% Функция определяет видимые КА для шлюзовых станций (ШС)  в заданный момент
% времени. Функция возвращает массив данных с указанием номера КА, номера
% ШС, угол места и азимут антены, а так же доплеровский сдвиг частоты
% радиосигнала.

 % Константы
    earthRadius = 6378135;           % Экваториальный радиус Земли [м]
    elevationMinAngle = 25;          % Минимальный угол места [град]
    pointSpeedLight = 299792458;     % Скорость света [м/с]
    carrierFrequency= 433e6;         % Несущая частота [Гц]
        
    % Получение географических координат шлюзовых станций (ШС) в начальный момент времени [град]
    gatewaysCoordinatesGeoArray = jsondecode(fileread('gatewaysTest.json'));
    
    satellitesCount = constellation.totalSatCount;          % Число КА
    gatewaysCount = length(gatewaysCoordinatesGeoArray);    % Число ШС
    
    % Получение инерциальных координат ШС на текущую эпоху  [м]
    for gatewayIdx = 1: gatewaysCount
        gatewaysCoordinatesIcs(gatewayIdx, 1:3) = transitionGeoToIcs([gatewaysCoordinatesGeoArray(gatewayIdx).lat, ...
                                                                      gatewaysCoordinatesGeoArray(gatewayIdx).lon],...
                                                                      epochArray(epochIdx));
    end
        
    % Получение инерциальных координат КА на момент epochIdx [м]
    satellitesCoordinatesIcsArray = constellation.state.eci(:, :, epochIdx);
    
    % Получение географических координат КА на момент epochIdx [м]
    for satellitesIdx = 1:satellitesCount
        satellitesCoordinatesGeoArray(satellitesIdx, 1:2) = transitionIcsToGeo(satellitesCoordinatesIcsArray(satellitesIdx,:)); % Географические координаты КА 
    end
    
    maxDistanceGroupList = zeros(2,1); % Инициализация списка максимальной дистанции видимости для группировок
    
    % Нахождение максимальной дистанции между КА и ШС в зоне видимости [м]
    maxDistanceGroupList(1) = calcMaxDistance(earthRadius, constellation.groups{1}.altitude*1000, elevationMinAngle);
    maxDistanceGroupList(2) = calcMaxDistance(earthRadius, constellation.groups{2}.altitude*1000, elevationMinAngle);
    firstGroupCount  = constellation.groups{1}.satsPerPlane * constellation.groups{1}.planeCount; % Число КА в первой группе
    secondGroupCount = constellation.groups{2}.satsPerPlane * constellation.groups{2}.planeCount; % Число КА во второй группе
    
    pairSatGatArray     = []; % Инициализацияя массива номеров пар КА и ШС 
    distancesSatGatList = []; % Инициализацияя дистанций между КА и ШС
    distanceSatGatPoint = 0;  % Инициализацияя буферной переменной дистанции между текущими КА и ШС  
    maxDistancePoint    = 0;  % Инициализацияя буферной переменной максимальной дистанции между текущими КА и ШС  
    
    % Нахождение дистанции между КА и ШС [м]
    for satellitesIdx = 1: satellitesCount
        for gatewayIdx = 1:  gatewaysCount           
            distanceSatGatPoint = sqrt((satellitesCoordinatesIcsArray(satellitesIdx, 1) - gatewaysCoordinatesIcs(gatewayIdx, 1)) ^ 2 ...
                                     + (satellitesCoordinatesIcsArray(satellitesIdx, 2) - gatewaysCoordinatesIcs(gatewayIdx, 2)) ^ 2 ...
                                     + (satellitesCoordinatesIcsArray(satellitesIdx, 3) - gatewaysCoordinatesIcs(gatewayIdx, 3)) ^ 2 );
                                   
            if (satellitesIdx <=  firstGroupCount)    % Выбор максимальной дистанции от группировки
                maxDistancePoint =  maxDistanceGroupList(1);
            else
                maxDistancePoint =  maxDistanceGroupList(2);
            end
        
            if (distanceSatGatPoint <= maxDistancePoint) % Условие видимости
                pairSatGatArray     = [pairSatGatArray; satellitesIdx, gatewayIdx];
                distancesSatGatList = [distancesSatGatList; distanceSatGatPoint];
            end           
        end % Окончание цикла по всем спутникам
    end % Окончание цикла по всем ШС
    
    pairSatGatCount         = length(pairSatGatArray);   % Число пар КА и ШС
    antennaOrientationArray = zeros(pairSatGatCount,2);  % Инициализацияя массива ориентации антены для пар КА и ШС (угол места и азимут антены)
    auxiliaryVector         = zeros(1,3);                % Инициализацияя буферной переменной разности векторов
     
    % Нахождение ориентации антены для пар КА и ШС
    for pairSatGatIdx = 1:pairSatGatCount
        auxiliaryVector = gatewaysCoordinatesIcs   (pairSatGatArray(pairSatGatIdx,2),:) ...
                        - satellitesCoordinatesIcsArray (pairSatGatArray(pairSatGatIdx,1),:); % Вспомогательный вектор
         % Вычислене угла места [град]
         antennaOrientationArray(pairSatGatIdx, 2) =  acosd ( satellitesCoordinatesIcsArray(pairSatGatArray(pairSatGatIdx,1),:) * auxiliaryVector' / ...
                                                            ( norm(auxiliaryVector) * norm(satellitesCoordinatesIcsArray(pairSatGatArray(pairSatGatIdx,1),:)))) ...
                                                             - 90;
         % Определение азимута антены [град]   
         antennaOrientationArray(pairSatGatIdx, 1) = 180 + atand(tand ( gatewaysCoordinatesGeoArray(pairSatGatArray(pairSatGatIdx,2)).lon ...
                                                                      - satellitesCoordinatesGeoArray(pairSatGatArray(pairSatGatIdx,1))  ...
                                                                      / sind(gatewaysCoordinatesGeoArray(pairSatGatArray(pairSatGatIdx,2)).lat)));                                                
    end % Окончание цикла по парам КА и ШС
    
    % Определение Доплеровского смещения ( Абсолютной величины. 
    % Для определения знака требуется знать производную угла места, которую
    % данная функция не вычисляет)
    
    satelliteSpeedPoint = 0;                        % Инициализацияя переменной скорости текущего КА 
    orbitRadiusPoint    = 0;                        % Инициализация переменной радиуса орбиты для текущего КА 
    orbitHeghtPoint     = 0;                        % Инициализация переменной радиуса орбиты для текущего КА
    frequencyOffsetList = zeros(pairSatGatCount,1); % Инициализацияя массива сдвига частоты
    satelliteViewAngle  = 0;                        % Угол обзора текущего КА   
    
    for pairSatGatIdx = 1:pairSatGatCount
        
      if (pairSatGatArray(pairSatGatIdx,1) <=  firstGroupCount)    % Подбор высоты и радиуса орбиты [м]
          orbitRadiusPoint = earthRadius + constellation.groups{1}.altitude * 1000;
          orbitHeghtPoint  = constellation.groups{1}.altitude * 1000;
      else
          orbitRadiusPoint = earthRadius + constellation.groups{2}.altitude * 1000;
          orbitHeghtPoint  = constellation.groups{2}.altitude*1000;
      end 
      
      satelliteSpeedPoint = sqrt(constellation.earthGM / orbitHeghtPoint); % Скорость спутника [м/с]
      satelliteViewAngle = asind( (earthRadius * sind (90 + antennaOrientationArray(pairSatGatIdx, 2)))/( earthRadius + orbitHeghtPoint) ); % Угол обзора текущего КА
      % Смещение частоты
      frequencyOffsetList(pairSatGatIdx,1) = 2 * carrierFrequency * satelliteSpeedPoint * cosd(90 - satelliteViewAngle) / pointSpeedLight;
    end % Окончание цикла по парам КА и ШС

   % Результируюшая структура данных 
    visiblesArray = struct('satelliteIdx',    pairSatGatArray(:,1),          'gatewayIdx',     pairSatGatArray(:,2), ...
                           'elevationAngle',  antennaOrientationArray(:, 2), 'antennaAzimut',  antennaOrientationArray(:, 1), ...
                           'distancesSatGat', distancesSatGatList(:,1),      'frequencyOffset', frequencyOffsetList(:,1)); 
  
end