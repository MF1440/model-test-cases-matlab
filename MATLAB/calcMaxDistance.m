function maxDistance = calcMaxDistance(earthRadius, orbitHeight, elevationMinAngle)
% Функция находит максимальную дистанцию видимости(наклонную дальность) между КА и ШС. 
% earthRadius - радиус Земли,
% orbitHeight - высота орбиты КА,
% elevationMinAngle - минимальный угол места ШС.
% Maxdistance - максимальная дистанция видимости (соответствует
% минимальному углу места)
 
    satelliteViewAngle = asind( (earthRadius * sind (90 + elevationMinAngle))/( earthRadius + orbitHeight) ); % Угол обзора КА
    centralEarthAngle  = 180 - (90 + elevationMinAngle) - satelliteViewAngle; % Центральный угол из центра Земли
    maxDistance        = earthRadius * sind(centralEarthAngle) / sind(satelliteViewAngle); % Максимальная наклонная дальность    
end