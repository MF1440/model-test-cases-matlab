clc
clear all

tic
minOrbitCount = 10;        % минимальное рассматриваемое количество орбит 
maxOrbitCount = 70;        % максимальное рассматриваемое количество орбит
minSatCountPerOrbit = 10;  % минимальное рассматриваемое количество спутников на орбите
maxSatCountPerOrbit = 70;  % максимальное рассматриваемое количество спутников на орбите
phi = randi(360);          % азимутальный угол в сферической системе координат для рандомной точки на поверхности земли
theta = randi(180) - 90;   % полярный угол в сферической системе координат для рандомной точки на поверхности земли
epoch = 0;                 % эпоха, в который считается вероятность покрытия точки (phi, theta)
q = 0.9;                   % вероятность безотказной работы каждого спутника в заданный момент времени
timeMoment = 0.5;          % момент времени, в который расчитывается вероятность, заданный как доля периода, изменяется от 0 до 1


bestConfiguration = findBestConfiguration(minOrbitCount, ...
                                        maxOrbitCount, ...
                                        minSatCountPerOrbit, ...
                                        maxSatCountPerOrbit, ...
                                        epoch);


if ~isempty(bestConfiguration)
    probPointCover = calcProbPointCover(bestConfiguration, ...
                                        deg2rad(phi), ...
                                        deg2rad(theta), ...
                                        epoch, ...
                                        q, ...
                                        2*pi*timeMoment);
    
    disp(['Максимальным покрытием является покрытие земли с 60-й параллели северного полушария до 60-й параллели южного полушария'])
    disp(['Минимальное количество спутников необходимое для максимального покрытия: ', ...
        int2str(bestConfiguration(1) * bestConfiguration(2))])
    disp(['Количество спутников на орбите составляет ', int2str(bestConfiguration(1)), ...
        ' при количестве орбит ', int2str(bestConfiguration(2))])
    disp(['Вероятнсть того, что в момент времени t = ', num2str(timeMoment), ...
        ' точка с координатами phi = ', num2str(phi), ', theta = ',  num2str(theta), ...
        ' не будет иметь покрытие: ', num2str(probPointCover)])
else
    disp(['Полное покрытие не найдено'])
end

toc
function bestConfiguration = findBestConfiguration(minOrbitCount, ...
                                                    maxOrbitCount, ...
                                                    minSatCountPerOrbit, ...
                                                    maxSatCountPerOrbit, ...
                                                    epoch)
    % Находит лучшую конфигурацию количества орбит и спутников необходимую для полного
    % покрытия земли с 60 параллели северного получария до 60 параллели
    % южного полушария
    
    earthCoverage = EarthCoverage();
    arraySatOrbitCount = sortSatOrbitCount(minOrbitCount, ...
                        maxOrbitCount, ...
                        minSatCountPerOrbit, ...
                        maxSatCountPerOrbit);
    minSatNeaded = ceil(2 * sin(deg2rad(AstroConstants.inclination)) / ...
                                        (1 - cos(earthCoverage.gammaAngle)));  % минимально возможное количество спутников для покрытия земли без перекрытий
    [~, minIdx] = max(arraySatOrbitCount(3,:)>=minSatNeaded);  

    bestConfiguration = [];
    
    for indexSatOrbitCount=minIdx:length(arraySatOrbitCount)  % цикл по парам (количество спутников, количество орбит) начинается с минимально возможного значения
        orbitCount = arraySatOrbitCount(1,indexSatOrbitCount);
        satPerOrbitCount = arraySatOrbitCount(2,indexSatOrbitCount);
        disp(['Индекс массива поиска ', int2str(indexSatOrbitCount), ': ',...
            int2str(satPerOrbitCount) ,'x', int2str(orbitCount)])
        for phase = 1:orbitCount
            flag = earthCoverage.checkEarthCoverage(satPerOrbitCount, orbitCount, phase, epoch);
            if(flag)
                bestConfiguration = [satPerOrbitCount, orbitCount, phase];
                break  
            end
        end  % конец цикла по фазовым сдвигам 

        if(flag)
            disp(['орбит ', int2str(orbitCount), ' спутников ', ...
                int2str(satPerOrbitCount), ' фазовый сдвиг ', int2str(phase), ...
            ' спутников (всего ', int2str(satPerOrbitCount * orbitCount),')'])
            break  
        end  
        disp(['Нет полного покрытия для ', int2str(satPerOrbitCount * orbitCount), ...
            ' спутников (орбит ', int2str(orbitCount), ' спутников на орбите ' , int2str(satPerOrbitCount), ...
            ')'])
    end  % конец цикла по парам (количество спутников, количество орбит)
end

function arraySatOrbitCount = sortSatOrbitCount(minOrbitCount, ...
                                                maxOrbitCount, ...
                                                minSatCountPerOrbit, ...
                                                maxSatCountPerOrbit)
    % Сортирует массив из количества спутников на орбите и количества орбит 
    % по возрастанию общего числа спутников

    orbitsCountArray = minOrbitCount:maxOrbitCount;
    satPerOrbitArray = minSatCountPerOrbit:maxSatCountPerOrbit;
    totalSatCountArray = orbitsCountArray' * satPerOrbitArray;
    orbitsArray = repmat(orbitsCountArray, 1,length(satPerOrbitArray))';
    satArray = repelem(satPerOrbitArray, length(orbitsCountArray))';
    arraySatOrbitCount = sortrows([orbitsArray, satArray, totalSatCountArray(:)],3)';
end 

function probPointCover = calcProbPointCover(bestConfiguration, phi, theta, epoch, q, dPhi)
    % Считает вероятность покрытия заданной точки хотя бы одним КА в момент времени epoch

    earthCoverage = EarthCoverage();
    earthCoverage.getSatSphereCoords(bestConfiguration(1), ...
                                bestConfiguration(2), ...
                                bestConfiguration(3), ...
                                epoch, dPhi);
    coverageDensity = earthCoverage.checkPointCover(phi, theta);
    probPointCover = (1-q)^coverageDensity;    
end
