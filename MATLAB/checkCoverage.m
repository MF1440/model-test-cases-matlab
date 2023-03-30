clear

% Переменные параметры задачи
timePoint = 1000; % Момент времени, для которого будет решаться задача
alpha = 40; % угол обзора спутников
eps = 25; % минимальный угол места для связи со шлюзовой станцией
gatewayFilename = '../gatewaysTest.json'; % имя файла с описанием шлюзовых станций


% опции запуска
DRAW_COVERAGE = false; % отрисовка зон покрытия и непокрытой точки
PRELIMINARY_TRIALS = true; % будем ли предварительно проверять покрытие на сетке
MEASURE_ELAPSED_TIME = false; % измеряем ли время работы
meshPointNumber = 1000; % приблизительное число точек сетки для проверки

if MEASURE_ELAPSED_TIME
    tic %#ok
end

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
        if MEASURE_ELAPSED_TIME
            toc %#ok
        end
        return
    end
end

% используя найденные сферические сегменты, найдем точки, покрытия
% спутниками которых достаточно для покрытия всей Земли
poi = getPoi(constellation.earthRadius, ...
    sphericalSegmentCenterCoordArray, sphericalSegmentNormal);
% poi = getPoi(constellation.earthRadius, ...
%     sphericalSegmentCenterCoordArray, sphericalSegmentNormal);

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
if MEASURE_ELAPSED_TIME
    toc %#ok
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

    eps = 1e-10; % малое возмущение

    % считаем предельный угол видимости, при котором лучи, опущенные под
    % данным углом еще пересекают землю
    maxAlpha = asind(radius ./ (radius + satAltitude)) - eps;
    alpha = alpha * ones(size(satAltitude));
    alpha(alpha > maxAlpha) = maxAlpha(alpha > maxAlpha);

    % коэффициенты квадратного уравнения
    a = ones(numel(satAltitude),1) .* (1.0 + tand(alpha) .* tand(alpha));
    b = - 2 * satAltitude - 2 * radius;
    c = satAltitude .^ 2 + 2 * radius * satAltitude;
    
    d = b .^ 2 - 4 * a .* c; % дискриминант
    
    % x1 = (- b + sqrt(d)) ./ (2 * a); это решение для другой стороны Земли
    x2 = (- b - sqrt(d)) ./ (2 * a);
    
    centerArray = satCoords ./ repmat(radius + satAltitude, 1, 3) .* ...
        repmat(radius + satAltitude - x2, 1, 3);
    radiusArray = tand(alpha) .* x2;
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

    % добавляем экватор к списку окружности - это граница двух подобластей
    centers = [centers; [0 0 0]];
    normalVectors = [normalVectors; [0 0 1]];

    % Нормируем вычисления
    centers = centers / r;
    rScaled = 1;

    % генерируем индексы пар плоскостей
    nPlanes = size(centers, 1);
    [planeIdxI, planeIdxJ] = meshgrid(1:nPlanes, 1:nPlanes);
    mask = planeIdxI < planeIdxJ;
    planeIdxI = planeIdxI(mask);
    planeIdxJ = planeIdxJ(mask);

    % Находим направляющие вектора прямых - пересечений пар плоскостей
    directions = cross(normalVectors(planeIdxI, :), ...
        normalVectors(planeIdxJ, :));
    
    % Отбраковываем непересекающиеся плоскости
    nonParallelPlanes = vecnorm(directions, 2, 2) > eps;
    planeIdxI = planeIdxI(nonParallelPlanes);
    planeIdxJ = planeIdxJ(nonParallelPlanes);
    directions = directions(nonParallelPlanes, :);
    nPlanes = sum(nonParallelPlanes);

    % Чтобы найти точку на прямой пересечения плоскостей нужно найти
    % решение системы из двух уравнений плоскостей: 
    rhs = [dot(normalVectors(planeIdxI, :), centers(planeIdxI, :), 2), ...
        dot(normalVectors(planeIdxJ, :), centers(planeIdxJ, :), 2)];
    linSys = [normalVectors(planeIdxI,:), normalVectors(planeIdxJ,:)];

    % Найдем решения системы зануляя поочередно каждую из переменных xyz
    % Занулять только x нельзя, т.к. система может не проходить через x=0
    % какие-то из 3 систем могут быть вырожденными, решение будет inf
    solutionsAll = get3Sol(linSys, rhs);

    % Найдем индекс решения которое нам подходит (1 - x=0, 2 - y=0, 3 - z=0)
    indBool = abs(directions) > eps;
    indBool(:,2) = indBool(:,2) & ~indBool(:,1);
    indBool(:,3) = indBool(:,3) & ~indBool(:,1) & ~indBool(:,2);
    [~,ind] = max(indBool,[],2);

    % достанем из вектора solutionsAll нужное нам решение по индексам из ind
    subscriptI = sub2ind(size(solutionsAll), ...
        (1:size(solutionsAll,1))', ind * 2 - 1);
    subscriptJ = sub2ind(size(solutionsAll), ...
        (1:size(solutionsAll,1))', ind * 2);
    pt = [solutionsAll(subscriptI), solutionsAll(subscriptJ)];
    
    % добавим 0, который подставляли вместо одной из переменных xyz в
    % соответствующий столбец. Для этого добавим 0 в третий столбец и потом
    % переставим в нужный
    pt = [pt, zeros(numel(subscriptI),1)];
    permIdx = zeros(nPlanes, 3);
    permIdx(ind == 1, :) = repmat([3 1 2], sum(ind == 1),1);
    permIdx(ind == 2, :) = repmat([1 3 2], sum(ind == 2),1);
    permIdx(ind == 3, :) = repmat([1 2 3], sum(ind == 3),1);
    permIdx = permIdx';
    rowInd = repelem(1:nPlanes,3)';
    columnInd = permIdx(:);
    linInd = sub2ind(size(pt), rowInd, columnInd);
    pt = pt(linInd);
    pt = reshape(pt,[3,numel(subscriptI)])';

    % Вычислим расстояние от центра Земли до прямую пересечения плоскостей
    distToLine = vecnorm(cross(pt, directions, 2), 2, 2) ./ vecnorm(directions, 2, 2);

    % Отбракуем прямые, не пересекающие сферу Земли
    mask = distToLine < 1 + eps;
    directions = directions(mask, :);
    pt = pt(mask,:);

    % Найдем точки пересечения прямой со сферой, для этого решим квадратные
    % уравнения относительно параметра - множителя направляющего вектора,
    % полученное из уравнения сферы
    polynom = [sum(directions .^ 2, 2), 2 * dot(directions, pt, 2), ...
                - rScaled ^ 2 + sum(pt .^ 2, 2)];
    tSol = arrayfun(@(rowIdx) roots(polynom(rowIdx,:)), 1:size(pt,1), 'UniformOutput', false);
    tSol = cell2mat(tSol)';

    % Отбракуем комплексные корни
    mask = imag(tSol(:,1)) == 0;
    tSol = tSol(mask,:);


    % найдем точки пересечения прямых со сферой
    intersrcPts1 = pt + directions .* tSol(:,1);
    intersrcPts2 = pt + directions .* tSol(:,2);
    points = [intersrcPts1; intersrcPts2];

    % Добавляем произвольную точку границы в множество рассматриваемых
    % точек
    points = [points; [1 0 0]];

    % Восстанавливаем масштаб
    points = points * r;
end
%--------------------------------------------------------------------------
function x = get3Sol(linSystems, rhs)
    % Векторизовано решаем группы из трех систем уравнений, если числа в
    % строке linSystems обозначить как a1, a2, ...,a6, то решаются
    % уравнения, а числа в строке rhs обозначить как b1 b2
    % [a2 a3; a5 a6] [x11 x12] = [b1 b2]'
    % [a1 a3; a4 a6] [x21 x22] = [b1 b2]'
    % [a1 a2; a4 a5] [x31 x32] = [b1 b2]'
    % Входные данные
    %   linSystems - матрицы групп систем уравнений, записанные по принципу
    %     выше, размер Nx6, N - количество троек систем уравнений
    %   rhs - правые части систем, размер Nx2
    % Выходные данные
    %   x - решения, записанные в строки, в каждой строке 
    %     [x11 x12 x21 x22 x31 x32], размер Nx6

    % Если зануляем первый аргумент
    det1 = linSystems(:,2) .* linSystems(:,6) - ...
        linSystems(:,3) .* linSystems(:,5);
    det1Hat1 = rhs(:,1) .* linSystems(:,6) - rhs(:,2) .* linSystems(:,3);
    det1Hat2 = linSystems(:,2) .* rhs(:,2) - rhs(:,1) .* linSystems(:,5);
    x1 = [det1Hat1 ./ det1, det1Hat2 ./ det1];

    det2 = linSystems(:,1) .* linSystems(:,6) - ...
        linSystems(:,3) .* linSystems(:,4);
    det2Hat1 = rhs(:,1) .* linSystems(:,6) - rhs(:,2) .* linSystems(:,3);
    det2Hat2 = rhs(:,2) .* linSystems(:,1) - rhs(:,1) .* linSystems(:,4);
    x2 = [det2Hat1 ./ det2, det2Hat2 ./ det2];

    det3 = linSystems(:,1) .* linSystems(:,5) - ...
        linSystems(:,2) .* linSystems(:,4);
    det3Hat1 = rhs(:,1) .* linSystems(:,5) - rhs(:,2) .* linSystems(:,2);
    det3Hat2 = rhs(:,2) .* linSystems(:,1) - rhs(:,1) .* linSystems(:,4);
    x3 = [det3Hat1 ./ det3, det3Hat2 ./ det3];

    x = [x1 x2 x3];
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

    centerToSat = vecnorm(satCoord,2,2);
    numSat = size(satCoord,1);
    for ptsIdx = 1:size(ptsArray,1)
        satToPt = vecnorm(satCoord - repmat(ptsArray(ptsIdx,:),numSat,1),2,2);
        centerToPt = repmat(vecnorm(ptsArray(ptsIdx,:)),numSat,1);

        alphaAngle = acosd(...
            (centerToSat .^ 2 + satToPt .^ 2 - centerToPt .^ 2) ./ ...
            (2 * centerToSat .* satToPt));
        epsAngle = acosd(...
            (centerToPt .^ 2 + satToPt .^ 2 - centerToSat .^ 2) ./ ...
            (2 * centerToPt .* satToPt));
        isCovered = any(epsAngle > 90 & alphaAngle < (alpha - eps));
        if ~isCovered
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
% function [allCovered, ptsNotCovered] = checkPtsCovered(ptsArray, satCoord, alpha)
%     % проверяем, все ли точки из массива ptsArray покрываются спутниками
%     % Выходные данные:
%     %   ptsArray - массив точек для проверки, Nx3
%     %   satCoord - координаты спутников, Mx3
%     %   alpha - угол видимости
%     % Выходные данные:
%     %   allCovered - логическая переменная, обозначающая "все точки из
%     %     массива ptsArray покрываются спутниками"
%     %   ptsNotCovered - пример точки, окрестность которой не покрыта
%     %   спутниками
% 
%     eps = 1e-10; % малое возмущение
% 
%     for ptsIdx = 1:size(ptsArray,1)
%         iCovered = false;
%         for satIdx = 1:size(satCoord,1)
%             % Вычисляем длины сторон треугольника: точка для проверки,
%             % спутник, центр Земли
%             centerToSat = vecnorm(satCoord(satIdx,:));
%             satToPt = vecnorm(ptsArray(ptsIdx,:) - satCoord(satIdx,:));
%             centerToPt = vecnorm(ptsArray(ptsIdx,:));
%     
%             % Вычисляем углы по теореме косинусов
%             alphaAngle = acosd((centerToSat ^ 2 + satToPt ^ 2 - centerToPt ^ 2) / ...
%                 (2 * centerToSat * satToPt));
%             epsAngle = acosd((centerToPt ^ 2 + satToPt ^ 2 - centerToSat ^ 2) / ...
%                 (2 * centerToPt * satToPt));
% 
%             % проверяем влазит ли точка в угол видимости и не лежит ли на
%             % обратной стороне Земли (есть ли прямая видимость)
%             if epsAngle > 90 && alphaAngle < (alpha - eps)
%                 iCovered = true;
%                 break
%             end
%         end
% 
%         % Если точка не покрыта ни одним из спутников, то прерываем
%         % вычисления, ответ на задачу получен
%         if ~iCovered
%             ptsNotCovered = ptsArray(ptsIdx,:); % сохраняем пример точки
%             allCovered = false;
%             return
%         end
%     end
% 
%     % все точки закончились, т.е. все покрыты спутниками
%     allCovered = true;
%     ptsNotCovered = [nan, nan, nan];
% end
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