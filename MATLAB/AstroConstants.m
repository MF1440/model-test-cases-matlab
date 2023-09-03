classdef AstroConstants

    properties(Constant)
        % earth constants
        earthRadius = 6378135;           % Экваториальный радиус Земли [m]
        earthGM = 3.986004415e+14;       % Гравитационный параметр Земли [m3/s2]
        earthJ2 = 1.082626e-3;           % Вторая зональная гармоника геопотенциала
        
        % group constants
        altitude = 600;
        inclination = 60; 
        startRaan = 0;
        maxRaan = 360;
        alphaAngle = 50;
    end
end