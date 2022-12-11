function createConstellationFromJson(constellation, jsonFileName, constellationName)
    str = fileread(jsonFileName);
    data = jsondecode(str);
    dataThis = [];
    
    dataThis = data(strcmpi(constellationName, {data(:).name}));
    
    if isempty(dataThis)
        disp('Группировка не найдена в файле');
        return
    end
    
    for groupIdx = size(dataThis.Walkers, 1):-1:1
        % Наклонение орбитальной плоскости
        thisGroup.inclination = deg2rad(dataThis.Walkers(groupIdx, 1));
    
        % Число КА в каждой орбитальной плоскости группы
        thisGroup.satsPerPlane = dataThis.Walkers(groupIdx, 2);
    
        % число орбитальных плоскостей в группе
        thisGroup.planeCount = dataThis.Walkers(groupIdx, 3);
    
        % фазовый сдвиг по аргументу широты между КА в соседних плоскостях
        thisGroup.f = dataThis.Walkers(groupIdx, 4);
    
        % высота орбиты
        thisGroup.altitudeKilometers = dataThis.Walkers(groupIdx, 5);
    
        % максимум прямого восхождения восходящего узла (при распределении орбитальных плоскостей)
        thisGroup.maxRaan = deg2rad(dataThis.Walkers(groupIdx, 6));
    
        % прямое восхождение восходящего узла для первой плоскости
        thisGroup.startRaan = deg2rad(dataThis.Walkers(groupIdx, 7));
    
        thisGroup.totalSatCount = thisGroup.satsPerPlane * thisGroup.planeCount;
    
        constellation.groups{groupIdx} = thisGroup;
        constellation.totalSatCount = constellation.totalSatCount + thisGroup.totalSatCount;
    end
end

