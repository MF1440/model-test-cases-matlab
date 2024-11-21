function hasCoverage = isGlobalCoverage(satellitePositions, alpha)
    earthRadius = 6378135; 

    % создаем сетку точек на поверхности Земли 
    latitudes = -90:10:90; 
    longitudes = -180:10:180;

    % метка покрытия для каждой точки на поверхности Земли
    covered = false(length(latitudes), length(longitudes));

    for epochIdx = 1:size(satellitePositions, 3)
        for satIdx = 1:size(satellitePositions, 1)
            % получаем положение спутника в данный момент времени
            satPos = satellitePositions(satIdx, :, epochIdx);

            % преобразуем координаты спутника в географические
            [satLat, satLon, satAlt] = ecef2GeographicSystem(satPos(1), satPos(2), satPos(3));

            % вычисляем длину дуги от точки проекции спутника на Землю и
            % точки пересечения луча спутника и Земли под углом раствора
            coverageDistance = earthRadius * asin((earthRadius / (earthRadius + satAlt)) * sin(deg2rad(alpha)));

            % проверяем видимость каждой точки на поверхности Земли
            for latIdx = 1:length(latitudes)
                for lonIdx = 1:length(longitudes)
                    % вычисляем расстояние между спутником и точкой на поверхности
                    pointLat = latitudes(latIdx);
                    pointLon = longitudes(lonIdx);

                    distance = haversineDistance(satLat, satLon, pointLat, pointLon, earthRadius);

                    % если точка находится в зоне видимости спутника, отмечаем ее как покрытую
                    if distance <= coverageDistance
                        covered(latIdx, lonIdx) = true;
                    end
                end
            end
        end
    end

    % проверяем, покрыты ли все точки на поверхности Земли
    hasCoverage = all(covered(:));
end
