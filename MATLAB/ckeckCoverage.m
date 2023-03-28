clear

% Переменные параметры задачи
timePoint = 1000; % Момент времени, для которого будет решаться задача
alpha = 40; % угол обзора спутников
eps = 25; % минимальный угол места для связи со шлюзовой станцией
gatewayFilename = '../gatewaysTest.json'; % имя файла с описанием шлюзовых станций


% опции запуска
DRAW_COVERAGE = false; % отрисовка зон покрытия и непокрытой точки
PRELIMINARY_TRIALS = true; % будем ли предварительно проверять покрытие на сетке
meshPointNumber = 1000; % приблизительное число точек сетки для проверки


% Создание объекта типа Constellation, инициализация параметрами группировки Stalink из конфига
constellation = Constellation('Starlink');

% Вычисление элементов орбиты для всех КА в начальный момент
constellation.getInitialState();

% расчёт положений всех КА выбранный момент времени
constellation.propagateJ2(timePoint);

% записываем координаты и высоты спутников в отдельный массив
satelliteCoordArray = (constellation.state.eci(:, :, 1));
satelliteAltitudeArray = (constellation.state.elements(:, 1)) - ...
    constellation.earthRadius;


% Для спутников определяем возможность связи со шлюзовой станцией
% Сюда уже входит добавление спутников из одной орбитальной плоскости со
% спутником в прямой видимости от шлюзовой станции
isVisible = checkGatewayAvailability(satelliteCoordArray, gatewayFilename, ...
    eps, [cell2mat(constellation.groups(1:2)).totalSatCount], ...
    [cell2mat(constellation.groups(1:2)).satsPerPlane], ...
    constellation.earthRadius);

% находим основания всех сферических сегментов, образуемых областями 
% видимости спутников. Т.к. далее будем использовать только плоскости 
% оснований, то радиусы оснований не считаем
[sphericalSegmentCenterCoordArray, rad, sphericalSegmentNormal] = ...
    getSphericalSegments(satelliteCoordArray(isVisible,:), ...
    satelliteAltitudeArray(isVisible,:), alpha, constellation.earthRadius);

% для ускорения сначала проверим покрытие 
if PRELIMINARY_TRIALS
    [meshCovered, ptsNotCovered] = ...
        checkCoverageOnRegularMesh(meshPointNumber, ...
        constellation.earthRadius, satelliteCoordArray(isVisible,:), alpha);
    if ~meshCovered
        disp('Есть необслуживаемые точки')
        if DRAW_COVERAGE
            drawCoverage(constellation.earthRadius, ... 
                sphericalSegmentCenterCoordArray, sphericalSegmentNormal);%#ok
            scatter3(ptsNotCovered(1),ptsNotCovered(2),ptsNotCovered(3),'f')
        end
        return
    end
end

% используя найденные сферические сегменты, найдем точки, покрытия
% спутниками которых достаточно для покрытия всей Земли
poi = getPoi(constellation.earthRadius, ...
    sphericalSegmentCenterCoordArray, sphericalSegmentNormal);

[allCovered, ptsNotCovered] = checkPtsCovered(poi, ...
    satelliteCoordArray(isVisible,:), alpha);
if allCovered
    disp('Все точки Земли обслуживаются')
else
    disp('Есть необслуживаемые точки')
end

% отрисовка зон покрытия, если покрытие не обеспечивается отрисовывается
% также пример точки с необслуживаемой окрестностью
if DRAW_COVERAGE
    drawCoverage(constellation.earthRadius, ... 
        sphericalSegmentCenterCoordArray, sphericalSegmentNormal);%#ok
    scatter3(ptsNotCovered(1),ptsNotCovered(2),ptsNotCovered(3),'f')
end


%% локальные функции

function isVisible = checkGatewayAvailability(satelliteCoordArray, gatewayFilename, eps, ...
    satGroupSizeArray, satsPerPlaneArray, radius)
    % Определение доступности связи со шлюзовой станцией для спутников
    % Входные данные: 
    %   satelliteCoordArray - массив с координатами спутников размера Nx3,
    %     где N - количество спутников
    %   gatewayFilename - имя файла с описанием шлюзовых станций
    %   eps - минимальный угол места для связи со шлюзовой станцией
    %   satGroupSizeArray - массив размеров групп спутников размера Mx1,
    %     где M - количество групп
    %   satsPerPlaneArray - количество спутников на орбиту в каждой группе,
    %     размер Mx1
    %   radius - радиус Земли, скаляр
    % Выходные данные:
    %   isVisible - логический массив Nx1, принимающий значение true для
    %     спутников, имеющих доступ к шлюзовой станции

    % чтение файла с описанием шлюзовых станций
    str = fileread(gatewayFilename);
    data = jsondecode(str);
    % перевод широты и долготы в декартовы координаты
    gatewayCoordArray  = [radius * cos([data.lat])' .* cos([data.lon])', ...
                          radius * cos([data.lat])' .* sin([data.lon])', ...
                          radius * sin([data.lat])'];

    % Находим спутники, имеющие связь со шлюзовой станцией
    isVisible = false(size(satelliteCoordArray,1),1);
    for satIdx = 1:size(satelliteCoordArray,1)
        for gatewayIdx = 1:size(gatewayCoordArray,1)
            % вычисляем угол места
            angle = getElevationAngle(satelliteCoordArray(satIdx,:),...
                gatewayCoordArray(gatewayIdx,:));
            if angle >= eps
                % добавляем КА, в т. ч. связанные с найденым КА с 
                % межспутниковой линией связи
                isVisible(getNeighbouringSatIdx(...
                    satIdx, satGroupSizeArray, satsPerPlaneArray)) = true;
            end
        end
    end

end
%--------------------------------------------------------------------------
function angle = getElevationAngle(satelliteCoords, viewerCoords)
    % Вычисление угла места
    % Входные данные:
    %   satelliteCoords - координаты спутника, вектор 1x3
    %   viewerCoords - координаты точки на Земле, для которой вычисляем 
    %     угол места, вектор 1x3
    % Выходные данные:
    %   angle - угол места в градусах

    viewerToSatellite = satelliteCoords - viewerCoords;
    viewerToEarthCenter = zeros(size(viewerCoords)) - viewerCoords;
    angle = atan2(norm(cross(viewerToSatellite,viewerToEarthCenter)), ...
        dot(viewerToSatellite,viewerToEarthCenter));
    angle = rad2deg(angle) - 90;
end
%--------------------------------------------------------------------------
function neighbouringSatIdxArray = getNeighbouringSatIdx(satIdx, ...
    satGroupSizeArray, satsPerPlaneArray)
    % Находим индексы всех спутников в той же орбитальной плоскости, что и 
    % спутник с индексом satIdx
    % Входные данные:
    %   satIdx - индекс спутника, скаляр
    %   satGroupSizeArray - массив количества спутников в каждой группе
    %     спутников, размер 1xN, N - количество групп
    %   satsPerPlaneArray - массив количества спутников на орбиту в каждой
    %     группе, размер 1xN
    % Выходные данные:
    %   neighbouringSatIdxArray - массив индексов спутников из той же 
    %     группы и той же орбитальной плоскости, что и satIdx, 
    %     размер 1xsatsPerPlaneArray

    groupSzCumsum = cumsum(satGroupSizeArray);
    groupIdx = find(satIdx <= groupSzCumsum, 1);
    groupSzCumsum = [0, groupSzCumsum];
    
    inPlaneIdx = mod(satIdx - groupSzCumsum(groupIdx) - 1, ...
        satsPerPlaneArray) + 1;
    
    leftNeighbour = satIdx - inPlaneIdx + 1;
    rightNeighbour = leftNeighbour + satsPerPlaneArray - 1;

    neighbouringSatIdxArray = leftNeighbour : rightNeighbour;
end
%--------------------------------------------------------------------------
function [centerArray, radiusArray, normalArray] = ...
    getSphericalSegments(satCoords, satAltitude, alpha, radius)
    % для набора спутников возвращает вычисляем отсекаемые от Земли их 
    % конусами сферические сегменты
    % Входные данные:
    %   satCoords - координаты спутников, размер Nx3, N - количество 
    %     спутников
    %   satAltitude - высоты спутников над землей, размер Nx1
    %   alpha - угол видимости
    %   radius - радиус сферы (Земли)
    % Выходные данные:
    %   centerArray - массив координат центров оснований сферических 
    %     сегментов, размер Nx3
    %   radiusArray - массив радиусов оснований, размер Nx1
    %   normalArray - массив нормальных векторов в плоскостям оснований, 
    %     размер Nx3

    % расстояние от центра сферы до центра основания вычисляется из 
    % геометрических соображений: если через x обозначить расстояние от
    % спутника до центра основания, то записывая теорему Пифагора для
    % треугольника, образованного центром сферы, центром основания и крайней
    % точкой видимости, получим квадратное уравнение для x с коэффициентами:
    a = ones(numel(satAltitude),1) * (1.0 + tand(alpha) * tand(alpha));
    b = - 2 * satAltitude - 2 * radius;
    c = satAltitude .^ 2 + 2 * radius * satAltitude;
    
    d = b .^ 2 - 4 * a .* c; % дискриминант
    
    % x1 = (- b + sqrt(d)) ./ (2 * a); это решение для другой стороны Земли
    x2 = (- b - sqrt(d)) ./ (2 * a);
    
    centerArray = satCoords ./ repmat(radius + satAltitude, 1, 3) .* ...
        repmat(radius + satAltitude - x2, 1, 3);
    radiusArray = tand(alpha) * x2;
    % нормальный вектор нормируем
    normalArray = satCoords ./ repmat(vecnorm(satCoords,2,2),1,3);
end
%--------------------------------------------------------------------------
function points = getPoi(r, centers, normalVectors)
    % Находим множество точек, проверки покрытия которых достаточно для
    % покрытия все Земли. Множество состоит из 1) точек пересечения 
    % оснований сферических сегментов, отсекаемых спутниками, 2) точек
    % пересечения оснований сферических сегментов с сечением Земли в
    % плоскости z=0 и 3) одной из точек этого сечения
    % Входные данные:
    %   r - радиус Земли, скаляр
    %   centers - массив центров оснований сферических сегментов, размер
    %     Nx3, где N - количество спутников (тех, что связаны со шлюзом)
    %   normalVectors - массив нормальных векторов оснований, размер Nx3

    eps = 1e-10; % малое возмущение

    % добавляем окружность, делящую сферу пополам - это граница подобластей
    centers = [centers; [0 0 0]];
    normalVectors = [normalVectors; [0 0 1]];
    
    % нормируем вычисления
    centers = centers / r;
    rScaled = 1;
    
    % ищем попарные пересечения оснований сферических сегментов
    numPoints = 0;
    points = zeros(size(centers,1) ^ 2, 3);
    for centerIdxI = 1:size(centers,1)
        for centerIdxJ = (centerIdxI + 1):size(centers,1)
            % ищем уравнение прямой проходящей через две плоскости через 
            % направляющий вектор :pt + t * direction

            % направляющий вектор - векторное произведение нормальных
            % векторов пересекающихся плоскостей
            direction = cross(normalVectors(centerIdxI,:),normalVectors(centerIdxJ,:));

            % убираем параллельные плоскости
            if (abs(sum(direction)) < eps)
                continue;
            end

            % ищем какую-нибудь точку пересечений плоскостей, зануляя одну
            % из координат
            rhs = [dot(normalVectors(centerIdxI,:), centers(centerIdxI,:)); ...
                dot(normalVectors(centerIdxJ,:), centers(centerIdxJ,:))];
            if abs(direction(1)) > eps
                linSys = [normalVectors(centerIdxI,2:3); normalVectors(centerIdxJ,2:3)];
                sol = linSys \ rhs;
                pt = [0; sol]; 
            elseif abs(direction(2)) > eps
                linSys = [normalVectors(centerIdxI,[1,3]); normalVectors(centerIdxJ,[1,3])];
                sol = linSys \ rhs;
                pt = [sol(1); 0; sol(2)];
            else
                linSys = [normalVectors(centerIdxI,[1,2]); normalVectors(centerIdxJ,[1,2])];
                sol = linSys \ rhs;
                pt = [sol; 0]; 
            end
            
            % считаем расстояние от центра сферы до прямой, если оно больше
            % радиуса Земли, то пересечения основания сферических
            % сегментов не пересекаются 
            distToLine = vecnorm(cross(pt,direction)) / vecnorm(direction);
            if distToLine > (rScaled + eps) 
                continue;
            end
    
            % ищем пересечения прямой со сферой Земли
            polynom = [sum(direction .^ 2), 2 * dot(direction, pt'), ...
                - rScaled ^ 2 + sum(pt .^ 2)];
            tSol = roots(polynom);

            % добавляем действительные различающиеся решения
            if isreal(tSol(1))
                numPoints = numPoints + 1;
                points(numPoints,:) = pt' + tSol(1) * direction;
            end
            solutionsAreDifferent = abs(tSol(2) - tSol(1)) > eps;
            if isreal(tSol(2)) && solutionsAreDifferent
                numPoints = numPoints + 1;
                points(numPoints,:) = pt' + tSol(2) * direction;
            end
        end
    end
    
    % обрезаем лишнюю выделенную память и добавляем фиксированную точку
    points = points(1:numPoints,:);
    points = [points; [1 0 0]];
    points = points * r;
end
%--------------------------------------------------------------------------
function [allCovered, ptsNotCovered] = checkPtsCovered(ptsArray, satCoord, alpha)
    % проверяем, все ли точки из массива ptsArray покрываются спутниками
    % Выходные данные:
    %   ptsArray - массив точек для проверки, Nx3
    %   satCoord - координаты спутников, Mx3
    %   alpha - угол видимости
    % Выходные данные:
    %   allCovered - логическая переменная, обозначающая "все точки из
    %     массива ptsArray покрываются спутниками"
    %   ptsNotCovered - пример точки, окрестность которой не покрыта
    %   спутниками

    eps = 1e-10; % малое возмущение

    for ptsIdx = 1:size(ptsArray,1)
        iCovered = false;
        for satIdx = 1:size(satCoord,1)
            % Вычисляем длины сторон треугольника: точка для проверки,
            % спутник, центр Земли
            centerToSat = vecnorm(satCoord(satIdx,:));
            satToPt = vecnorm(ptsArray(ptsIdx,:) - satCoord(satIdx,:));
            centerToPt = vecnorm(ptsArray(ptsIdx,:));
    
            % Вычисляем углы по теореме косинусов
            alphaAngle = acosd((centerToSat ^ 2 + satToPt ^ 2 - centerToPt ^ 2) / ...
                (2 * centerToSat * satToPt));
            epsAngle = acosd((centerToPt ^ 2 + satToPt ^ 2 - centerToSat ^ 2) / ...
                (2 * centerToPt * satToPt));

            % проверяем влазит ли точка в угол видимости и не лежит ли на
            % обратной стороне Земли (есть ли прямая видимость)
            if epsAngle > 90 && alphaAngle < (alpha - eps)
                iCovered = true;
                break
            end
        end

        % Если точка не покрыта ни одним из спутников, то прерываем
        % вычисления, ответ на задачу получен
        if ~iCovered
            ptsNotCovered = ptsArray(ptsIdx,:); % сохраняем пример точки
            allCovered = false;
            return
        end
    end

    % все точки закончились, т.е. все покрыты спутниками
    allCovered = true;
    ptsNotCovered = [nan, nan, nan];
end
%--------------------------------------------------------------------------
function drawCoverage(r, centers, normals)
    % Отрисовка области покрытия
    % Входные данные:
    %   r - радиус Земли, скаляр
    %   centers - массив центров оснований сферических сегментов,
    %     отсекаемых спутниками, Nx3
    %   normals - нормальные векторы к плоскостям оснований, Nx3

    eps = 0.001; % сдвиг масштаба для отрисовки зон покрытия

    % рисуем Землю
    resolution = 500;
    [x,y,z]=sphere(resolution);
    x = x * r;
    y = y * r;
    z = z * r;
    hold off
    surf(x, y, z,'FaceColor','red','EdgeColor','None');
    hold on
    coordArray = zeros(3,resolution + 1, resolution + 1);
    coordArray(1,:,:) = x;
    coordArray(2,:,:) = y;
    coordArray(3,:,:) = z;

    % для каждого спутника орисовываем зону покрытия
    for satIdx = 1:size(centers,1)
        xSat = x * (1 + eps);
        ySat = y * (1 + eps);
        zSat = z * (1 + eps);


        centerArray = zeros(size(coordArray));
        centerArray(1,:,:) = repmat(centers(satIdx,1), resolution + 1, ...
            resolution + 1);
        centerArray(2,:,:) = repmat(centers(satIdx,2), resolution + 1, ...
            resolution + 1);
        centerArray(3,:,:) = repmat(centers(satIdx,3), resolution + 1, ...
            resolution + 1);
        normalArray = zeros(size(coordArray));
        normalArray(1,:,:) = repmat(normals(satIdx,1), resolution + 1, ...
            resolution + 1);
        normalArray(2,:,:) = repmat(normals(satIdx,2), resolution + 1, ...
            resolution + 1);
        normalArray(3,:,:) = repmat(normals(satIdx,3), resolution + 1, ...
            resolution + 1);
        relativePos = squeeze(dot(normalArray, coordArray - centerArray,1));
        xSat(relativePos < 0) = nan;
        ySat(relativePos < 0) = nan;
        zSat(relativePos < 0) = nan;
        surf(xSat, ySat, zSat,'FaceColor','green','EdgeColor','None');
    end
end
%--------------------------------------------------------------------------
function [meshCovered, ptsNotCovered] = ...
    checkCoverageOnRegularMesh(ptNum, r, satCoord, alpha)
    % Проверка покрытия Земли спутниками на регулярной сетке
    % Входные данные:
    %   ptNum - приблизительное количество точек сетки, скаляр
    %   r - радиус Земли
    %   satCoord - массив координат спутников, размер Nx3, где N -
    %   количество спутников
    %   alpha - уголвидимости
    % Выходные данные:
    %   meshCovered - логическая переменная, означающая покрыты ли все
    %     точки сетки спутниками
    %   ptsNotCovered - пример непокрытой точки (если есть)

    resolution = ceil(sqrt(ptNum));
    [x, y, z] = sphere(resolution);
    x = x * r;
    y = y * r;
    z = z * r;
    [meshCovered, ptsNotCovered] = checkPtsCovered([x(:), y(:), z(:)], ...
        satCoord, alpha);
end