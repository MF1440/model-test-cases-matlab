clc
clear 

tic

%создание объекта типа Constellation, инициализация параметрами группировки Stalink из конфига
constellation = Constellation('Starlink');

%вычисление элементов орбиты для всех КА в начальный момент
constellation.getInitialState();


for i = 1:length(constellation.groups)
%для удобства сохраняем необходимые параметры в отдельные переменные
altitudeGroup = constellation.groups{i}.altitude * 1000; % перевод высоты орбиты из км в м
earthRadius   = constellation.earthRadius;% Радиус Земли
satCount      = constellation.groups{i}.totalSatCount;% Количество спутников в группировке

% Определение минимального угла обзора, при котором выполняется глобальное покрытие
viewingAngle = findMinViewAngle(altitudeGroup,satCount,earthRadius);

disp(['Группировка спутников с высотой орбиты ' num2str(constellation.groups{i}.altitude)...
    ' км обеспечивает полное покрытие при угле обзора ' num2str(viewingAngle) ' градуса' ]);

end

toc

function viewingAngle = findMinViewAngle(altitude,numberOfSat,earthRadius)
% Функция выполняет поиск минимального угла обзора при котором площадь,
% покрываемая спутниковой группировкой больше или равна площади поверхности
% Земли.

% Площадь поверхности Земли
earthArea = 4*pi*earthRadius^2;

% Шаг угла
angleStep = 1.5*pi/180;
% Массив углов места
elevationAngleArray = pi/2:-angleStep:0;
%Создаем пустые массивы
centralAngle1  = zeros(length(elevationAngleArray)); % массив центральных углов
localArea11 = centralAngle1; % массив эквивалентных площадей локальной рабочей зоны для отдельного взятого ИСЗ
summuryArea1 = centralAngle1; % массив общей площади, покрываемой группировкой спутников

i=1;
    for elevationAngle = elevationAngleArray
     
         centralAngle1(i) = acos(cos(elevationAngle)/(1+altitude/earthRadius)) - (elevationAngle);
         localArea11(i) = 2*pi*earthRadius^2*(1-cos(centralAngle1(i)));
         summuryArea1(i) = localArea11(i)*numberOfSat;

         if summuryArea1(i)>= earthArea           
             gammaAngle = pi/2+elevationAngle;
             viewingAngle = asin((sin(gammaAngle)*earthRadius)/((altitude+earthRadius)));
             viewingAngle = viewingAngle*180/pi;
             break
         end
     
        i=i+1;
    end
end