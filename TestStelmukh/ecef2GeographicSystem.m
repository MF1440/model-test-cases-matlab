function [latitude, longitude, altitude] = ecef2GeographicSystem(x, y, z)
    % в качестве модели Земли здесь используется идеальная сфера
    earthRadius = 6378135; 

    % Вычисление долготы
    longitude = atan2(y, x);
    longitude = rad2deg(longitude);

    % Вычисление широты
    latitude = atan2(z, sqrt(x^2 + y^2));
    latitude = rad2deg(latitude);

    % вычисление высоты над идеальной сферой
    distanceFromCenter = sqrt(x^2 + y^2 + z^2);
    altitude = distanceFromCenter - earthRadius;


end