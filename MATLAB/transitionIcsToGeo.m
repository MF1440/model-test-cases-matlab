function coordinatesGeo = transitionIcsToGeo(coordinatesICS)
% Функция перехода из ИСК в географическую систему
% координат 
% coordinatesICS = [X,Y,Z] - координаты неподвижной ИСК 
% coordinatesGeo = [lat,lon] - Географические координаты,

    % Константы 
    Radian = 57.2957795; % 1 радиан в градусах
    
    % Получение географических координат 
    coordinatesGeo(1,1) = atand ( coordinatesICS(3) / sqrt(coordinatesICS(1)^2 + coordinatesICS(2)^2));
    coordinatesGeo(1,2) = acosd ( coordinatesICS(1) / sqrt(coordinatesICS(1)^2 + coordinatesICS(2)^2));

    % Учет квадранта 
    if (coordinatesICS(2) <= 0) 
        coordinatesGeo(1,2) = -coordinatesGeo(1,2) + 180;
    end
    
    coordinatesGeo(1,1) = coordinatesGeo(1,1); % lat    
    coordinatesGeo(1,2) = coordinatesGeo(1,2); % lon  
end