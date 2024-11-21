function distance = haversineDistance(lat1, lon1, lat2, lon2, radius)
    % вычисление расстояния между двумя точками на сфере по формуле гаверсинуса
    dLat = deg2rad(lat2 - lat1);
    dLon = deg2rad(lon2 - lon1);

    a = sin(dLat/2)^2 + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * sin(dLon/2)^2;
    c = 2 * atan2(sqrt(a), sqrt(1-a));

    distance = radius * c;
end