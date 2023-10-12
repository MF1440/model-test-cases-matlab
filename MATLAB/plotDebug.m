%% Отрисовка результата 

h = patch('faces', segmentTrnglDotNmbrArr, 'vertices', segmentCoordArrM);
set(h, 'FaceColor', 'w', 'EdgeColor', 'k')
hold on
axis equal off vis3d
% Все геостационарные спутники
plot3(geoSatCoordArrM(:, 1), ...
      geoSatCoordArrM(:, 2), ...
      geoSatCoordArrM(:, 3), '*g', 'MarkerSize', 20)
% Все спутники орбитальной группировки
plot3(satCoordArrM(:, 1), satCoordArrM(:, 2), satCoordArrM(:, 3), '.')

% Пройдёмся по всем сегментам
for segmentNumber = 1:segmentCount
    
    % Вынимаем информацию о наличии удоввлетворяющих условию спутников
    acceptableSat       = segmentRelSatStructArr(segmentNumber).satName;
    acceptableGeoSat    = segmentRelSatStructArr(segmentNumber).corrGeoSatName;
    
    % Если таковых нет - пропускаем процесс.
    if acceptableSat == 0
        continue
    end

    % Формируем массив точек, чтобы нарисовать вектора из центра на спутник
    satVecArr = zeros(numel(acceptableSat).*2, 3);
    satVecArr(1:2:end, :) = [satCoordArrM(acceptableSat, 1), ...
                             satCoordArrM(acceptableSat, 2), ...
                             satCoordArrM(acceptableSat, 3)];
    satVecArr(2:2:end, :) = [segmentBarycenterCoordArrM(segmentNumber, 1), ...
                             segmentBarycenterCoordArrM(segmentNumber, 2), ...
                             segmentBarycenterCoordArrM(segmentNumber, 3)].*ones(size(satVecArr(2:2:end, :)));
    % % Формируем массив точек, чтобы нарисовать вектора из центра на
    % геостационар
    geoSatVecArr = zeros(numel(acceptableSat).*2, 3);
    geoSatVecArr(1:2:end, :) = [geoSatCoordArrM(acceptableGeoSat, 1), ... 
                                geoSatCoordArrM(acceptableGeoSat, 2), ...
                                geoSatCoordArrM(acceptableGeoSat, 3)];
    geoSatVecArr(2:2:end, :) = [segmentBarycenterCoordArrM(segmentNumber, 1), ...
                                segmentBarycenterCoordArrM(segmentNumber, 2), ...
                                segmentBarycenterCoordArrM(segmentNumber, 3)].*ones(size(geoSatVecArr(2:2:end, :)));
    
    % Рисуем вектора
    plot3(satVecArr(:, 1), satVecArr(:, 2), satVecArr(:, 3), 'r')
    plot3(geoSatVecArr(:, 1), geoSatVecArr(:, 2), geoSatVecArr(:, 3), 'b')
    
end
