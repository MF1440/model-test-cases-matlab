%% Начало работы
clc
clear
allTaskTic = tic;

%% Генерация файла с мировыми константами

% Текущая директория
currentDir = fileparts(which(mfilename));
% Полное имя файла для проверки
filePath = fullfile(currentDir, 'constants.mat');

if ~exist(filePath, 'file')
    saveConstantToFile('constants.mat');
end

%% Определение констант задачи

geoSatCount     = 4;                    % Число геостационарных спутников.
currentTime     = ceil(10000.*rand());  % Время, когда необходимо провести исследование. [сек]
segmentCount    = 64000;                % На какое количество сегментов будет поделена Земля.
criteriaAngle   = deg2rad(2);           % Какой угол является критерием. 

%% Загрузка дополнительных мировых констант 

load('constants.mat', 'earthRadius', 'earthOmega', 'earthGeostatRadius');

%% Генерация орбитальной группировки 

% Создание объекта класса Constellation, инициализация параметрами группировки Stalink
% из конфига.
constellation = Constellation('Starlink');

% Вычисление элементов орбиты для всех КА в начальный момент времени.
constellation.getInitialState();

% Расчёт положений всех КА в заданный момент времени
constellation.propagateJ2(currentTime);

% Получение координат всех КА в текущий момент времени. 
satCoordArrM = squeeze(constellation.state.eci(:, :, 1));

% Число спутников в орбитальной группировке
satCount = constellation.totalSatCount;

%% Генерация геостационарных спутников

% Сетка азимутов спутников в начальный момент времени. 
geoSatAzimuthArrRad = [0: 2.*pi./geoSatCount: 2.*pi]';
geoSatAzimuthArrRad(end) = [];

% Азимуты геостационарных спутников в текущий момент времени. 
geoSatAzimuthArrRad = geoSatAzimuthArrRad + earthOmega.*currentTime;

% Декартовы координаты геостационарных спутников в текущий момент времни 
[geoSatXCoordArrM, ...
 geoSatYCoordArrM, ...
 geoSatZCoordArrM] = sph2cart(geoSatAzimuthArrRad, zeros(geoSatCount, 1), earthGeostatRadius);
geoSatCoordArrM = [geoSatXCoordArrM, geoSatYCoordArrM, geoSatZCoordArrM];

% Очищаем память от лишних для дальнейшей работы переменных
clear geoSatAzimuthArrRad geoSatXCoordArrM geoSatYCoordArrM geoSatZCoordArrM

%% Генерация секторов сферы

% Определение углов в сферической СК точек, равномерно распределённых по
% сфере; а также индексов точек, соответствующих триангуляции сферы.
[segmentElevationArrRad, ...
 segmentAzimuthArrRad, ...
 segmentTrnglDotNmbrArr] = sphereTriangulation(segmentCount);

% После определения индексов вершин треугольников, можно отслеживать
% эволюцию их положения в инерциальной СК с течением времени. 
% Поэтому определим азимуты равномерных точек в текущий момент времени. 
segmentAzimuthArrRad = segmentAzimuthArrRad + earthOmega.*currentTime;

% Переведём полученные координаты вершин треугольников секторов из
% сферических в декартовы
[segmentXCoordArrayM, ...
 segmentYCoordArrayM, ...
 segmentZCoordArrayM] = sph2cart(segmentAzimuthArrRad, segmentElevationArrRad, earthRadius);
segmentCoordArrM = [segmentXCoordArrayM, segmentYCoordArrayM, segmentZCoordArrayM];

% Определение координат барицентров треугольников. 
segmentBarycenterCoordArrM = (segmentCoordArrM(segmentTrnglDotNmbrArr(:, 1), :) + ...
                              segmentCoordArrM(segmentTrnglDotNmbrArr(:, 2), :) + ...
                              segmentCoordArrM(segmentTrnglDotNmbrArr(:, 3), :)) ./ 3;

% Определение перпендикуляров к сегментам. 
tempVec1Array = segmentCoordArrM(segmentTrnglDotNmbrArr(:, 2), :) - ...
                segmentCoordArrM(segmentTrnglDotNmbrArr(:, 1), :);
tempVec2Array = segmentCoordArrM(segmentTrnglDotNmbrArr(:, 3), :) - ...
                segmentCoordArrM(segmentTrnglDotNmbrArr(:, 1), :);
segmentPerpVecCoordArr = cross(tempVec1Array, tempVec2Array, 2);
segmentPerpVecCoordArr = segmentPerpVecCoordArr./vecnorm(segmentPerpVecCoordArr, 2) .* 10.^3;
% Домножили на 10^3, так как это не меняет направления вектора, но упрощает
% его учёт. 

clear tempVec1Array tempVec2Array

% На тот случай, если перпендикуляр оказался направлен внутрь Земли
isReversed = dot(segmentPerpVecCoordArr, segmentBarycenterCoordArrM, 2) < 0;
segmentPerpVecCoordArr(isReversed) = segmentPerpVecCoordArr(isReversed).*(-1);

% Очищаем память от лишних для дальнейшей работы переменных 
clear segmentElevationArrRad segmentAzimuthArrRad

%% Определение подходящих геостационарных спутников

% Определение векторов на геостационарные спутники из барицентров треугольников. 
tempSegmMatrix   = repmat(segmentBarycenterCoordArrM, 1, 1, geoSatCount);
tempGeoSatMatrix = repmat(geoSatCoordArrM, 1, 1, segmentCount);
tempGeoSatMatrix = permute(tempGeoSatMatrix, [3, 2, 1]);
segment2GeoSatVecArr = tempGeoSatMatrix - tempSegmMatrix;

clear tempGeoSatMatrix tempSegmMatrix

% Видны такие геостационарные спутники, что угол между нормалью сегмента и
% вектора от барицентра до спутника меньше pi/2, то есть
% \vec{a}*\vec{n}*cos(\theta) > 0.
% Полученная таблица соответствия геостационарных спутников сотам
% записывается в isVisible.

tempPerpVecArr = repmat(segmentPerpVecCoordArr, 1, 1, geoSatCount);
dotPerp2Segment2GeoSatVecArr = squeeze(dot(segment2GeoSatVecArr, tempPerpVecArr, 2)); 
isVisible = dotPerp2Segment2GeoSatVecArr > 0;

% Очищаем память от лишних для дальнейшей работы переменных 
clear dotPerp2Segment2GeoSatVecArr tempPerpVecArr 

%% Определение подходящих КА из орбитальной группировки

% Так как число геостацинарных спутников может быть больше 4 на порядок, то
% далее, с точки зрения памяти, выгодно вести вычисления при помощи цикла. 
% Будем проходиться по сотам, чтобы сразу заполнять для них массив
% структур. 

segmentRelSatStructArr(1:segmentCount, 1) = struct('satName', 0, 'relGeoSatName', 0);

segForTic = tic;

for segIdx = 1:segmentCount
    % Набор векторов от текущего сектора до видимых геостационаров.
    currSeg2GeoSatVecArr = segment2GeoSatVecArr(segIdx, :, isVisible(segIdx, :));
    
    % Если видимых геостационаров нет, то и считать нет смысла. 
    if isempty(currSeg2GeoSatVecArr)
        continue
    end

    % Определим вектора от барицентра сектора на спутники орбитальной группировки.
    segment2SatVecArr = satCoordArrM - segmentBarycenterCoordArrM(segIdx, :);
    
    % Для уменьшения нагрузки на память, можно предварительно отсеять те
    % спутники из группировки, что лежат от нормали на 92 (90 + 2 из условия)
    % градуса и более. С другой стороны, лишние операции в цикле увеличат
    % время исполнения. Поэтому делать этого не будем и сразу создадим временные 
    % переменные большего размера.
    tempSeg2GeoSatVecMatrix = repmat(currSeg2GeoSatVecArr, satCount, 1, 1);
    tempSeg2SatVecMatrix    = repmat(segment2SatVecArr, 1, 1, sum(isVisible(segIdx, :)));
    
    % Скалярное произведение и угол из него.
    dotSeg2Sat2GeoSatMatrix = dot(tempSeg2GeoSatVecMatrix, tempSeg2SatVecMatrix, 2);
    cosSeg2Sat2GeoSatMatrix = dotSeg2Sat2GeoSatMatrix ./ (vecnorm(tempSeg2GeoSatVecMatrix, 2, 2) .* vecnorm(tempSeg2SatVecMatrix, 2, 2));
    
    % Определяем, какие полученные косинусы удовлетворяют условию.
    isAcceptable = (cosSeg2Sat2GeoSatMatrix >= cos(criteriaAngle));
    
    % Получаем два списка: подходящих КА и соответствующих им
    % геостационаров.
    [acceptableSat, tempAcceptableGeoSatNumber] = find(squeeze(isAcceptable));
    tempCurrGeoSat = find(isVisible(segIdx, :));
    acceptableGeoSat = tempCurrGeoSat(tempAcceptableGeoSatNumber);
    
    % Если были найдены спутники, удовлетворяющие условию, то производится
    % запись в соотвтетствующую структуру из массива.
    if ~isempty(acceptableSat)
        segmentRelSatStructArr(segIdx).satName        = acceptableSat;
        segmentRelSatStructArr(segIdx).corrGeoSatName = acceptableGeoSat;
    end
end

disp(['Прогонка по секторам завершена за ', num2str(toc(segForTic)), ' сек.'])
disp(['Выполнение задачи завершено. Затрачено: ', num2str(toc(allTaskTic)), ' сек.'])




