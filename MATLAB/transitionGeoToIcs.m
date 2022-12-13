function coordinatesIcs = transitionGeoToIcs(coordinatesGeo,epoch)
% Функция перехода из географической в инерциальную систему
% координат 
% coordinatesGeo = [lat,lon] - Географические координаты,
% pointRadius - Радиус сферы. (можно брать радиус Земли или радиус орбиты)
% coordinatesICS = [X,Y,Z] - координаты неподвижной ИСК 

    % константы
    earthRadius = 6378135;           % Экваториальный радиус Земли [m]
    earthRotV = 7.292115e-5;         % Угловая скорость вращения Земли [рад/c]
 
    % Получение координат в неподвижной ИСК 
    coordinatesIcsn(1,1) = earthRadius ...
                           * cosd(coordinatesGeo(1)) ...
                           * cosd(coordinatesGeo(2)); % X
    coordinatesIcsn(1,2) = earthRadius ...
                           * cosd(coordinatesGeo(1)) ...
                           * sind(coordinatesGeo(2)); % Y
    coordinatesIcsn(1,3) = earthRadius ... 
                           * sind(coordinatesGeo(1)); % Z

    %  Угол поворота Земли на текущую эпоху
    angleRotationEarth = epoch * earthRotV;
  
    % Матрица перехода из неподвижной ИСК в подвижную ИСК
    transferMatrixToRotationICS = [cos(angleRotationEarth) -sin(angleRotationEarth)  0;...
                                   sin(angleRotationEarth)  cos(angleRotationEarth)  0;...
                                   0                        0                        1];
    % Получение координат в подвижной ИСК                              
    coordinatesIcs= coordinatesIcsn*transferMatrixToRotationICS';                         
end