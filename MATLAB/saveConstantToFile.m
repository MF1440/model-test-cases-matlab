function saveConstantToFile(fileName)
% Данная функция сохраняет в файл fileName константы, которые перечислены
% в теле, создавая, тем самым, удобный для использования справочник.
%
% На вход:
% 1. fileName - имя файла.

    arguments
        fileName {mustBeNonzeroLengthText} = 'constants.mat';
    end

    % Константы, которые необходимо сохранить для дальнейшей работы.
    earthRadius = 6378135;          % Экваториальный радиус Земли [m]
    earthGM     = 3.986004415e+14;  % Гравитационный параметр Земли [m3/s2]
    earthJ2     = 1.082626e-3;      % Вторая зональная гармоника геопотенциала
    earthOmega  = 2.*pi ./ (86164); % угловая скорость земли [рад/сек]
    
    % Производные константы 
    earthGeostatRadius = round(earthGM./earthOmega.^2).^(1./3);
    
    % Получение списка всех переменных в рабочем пространстве
    allVariables = who;

    % Исключение из списка имени самого файла
    saveVariables = setdiff(allVariables, fileName);

    % Сохранение переменных в файл .mat
    save(fileName, saveVariables{:});
end