function [elevationArrayRadian, ...
          azimuthArrayRadian, ...
          triangleDotNumberArray] = sphereTriangulation(segmentCount)
% Данная функция равномерно распределяет по единичной сфере такое
% количство точек, что при разбиении сферы на треугольники с вершинами
% в этих точках, их получается заданное число segmentCount. В том
% случае, если segmantCount нечётное, тогда генерируется segmantCount+1
% треуголников. 
%
% На вход:
% 1. segmentCount - число треуголников, которыми в результате должна быть
% покрыта сфера. 
%
% На выход:
% 1. elevationArrayRadian - массив зенитных углов равномерно наброшенных на
% сферу точек.
% 2. azimuthArrayRadian - массив азимутальных углов равномерно наброшенных
% на сферу точек. 
% 3. triangleDotNumberArray - массив размером (segmentCount, 3), содержащий
% в себе номера точек треугольников. Каждая строка определяет три номера - 
% три вершины треугольника.
%
% Источник:
% 1. https://github.com/AntonSemechko/S2-Sampling-Toolbox

    arguments
        segmentCount {mustBePositive} = 64000;
    end
    
    funcTic = tic;

    % Число точек для размещения на сфере из характеристики Эйлера. 
    dotCount = 0.5.*segmentCount + 2;
    
    % Золотые сечение и угол. 
    constGoldenRatio = (1 + sqrt(5)) ./ 2;       
    constGoldenAngle = 2 .* pi .* (1 - 1 ./ constGoldenRatio);
    
    % Массивы зенитного и азимутального углов.
    % Сетка зенитных углов определяется так, чтобы она была равномерной по оси
    % Oz при переходе к декартовым координатам.
    % Сетка азимутальных углов равномерная.
    dotIndexesArray         = [0:(dotCount-1)]';
    elevationArrayRadian    = pi./2 - acos(1 - 2 .* dotIndexesArray ./ (dotCount-1));
    azimuthArrayRadian      =   dotIndexesArray .* constGoldenAngle;
    clear dotIndexesArray
    
    % Перевод из сферической СК в декартову, с последующим формированием
    % массива векторов. Это необходимо сделать заранее, так как функция
    % convhull работает только с точками в декартовой СК.
    [xArray, yArray, zArray] = sph2cart(azimuthArrayRadian, elevationArrayRadian, 1);
    dotCoordinateArray = [xArray, yArray, zArray];
    clear xArray yArray zArray
    
    % Разбиение сферы на треугольники с вершинами в заданных точках.
    triangleDotNumberArray = convhull(dotCoordinateArray);

    disp(['Триангуляция сферы выполнена за ', num2str(toc(funcTic)), ' сек.'])
end









