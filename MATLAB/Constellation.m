classdef Constellation < handle
    % класс Constellation содержит в себе информацию обо всех КА во всех заданных группировках

    properties
        totalSatCount = 0;
        groups = {};
        state;
    end

    methods

        function this = Constellation(varargin)
            if isempty(varargin)
                return
            end
            this.loadFromConfigFile(varargin{1});
        end

        function loadFromConfigFile(this, groupName)
            % Загружает параметры группировки из json файла.
            % Файл должен содержать в себе список группировок.
            % Каждая группировака описывается структурой, в которой указано название группировки и
            % кеплеровы элементы для каждого из эшелонов
            % Пример:
            %   [{"name": "Starlink","Walkers": [[53, 50, 32, 1, 500, 360, 0]]}]

            fileName = ['..' filesep 'constellationsTest.json'];
            str = fileread(fileName);
            data = jsondecode(str);
            dataThis = [];

            for i = 1:length(data)
                if strcmpi(data(i).name, groupName)
                    dataThis = data(i);
                    break
                end
            end

            if isempty(dataThis)
                disp(['Группировка ' groupName ' не найдена в файле']);
                return
            end

            for i = 1:size(dataThis.Walkers, 1)
                group.inclination = deg2rad(dataThis.Walkers(i, 1));        % наклонение орбитальной плоскости в градусах
                group.satsPerPlane = dataThis.Walkers(i, 2);				% число КА в каждой орбитальной плоскости группы
                group.planeCount = dataThis.Walkers(i, 3);					% число орбитальных плоскостей в группе
                group.phaseShift = dataThis.Walkers(i, 4);				    % фазовый сдвиг по аргументу широты между КА в соседних плоскостях. Варьируется от 0 до количества орбит в группе
                group.altitude = dataThis.Walkers(i, 5);					% высота орбиты в км
                group.maxRaan = deg2rad(dataThis.Walkers(i, 6));            % максимум прямого восхождения восходящего узла (при распределении орбитальных плоскостей) в градусах
                group.startRaan = deg2rad(dataThis.Walkers(i, 7));			% прямое восхождение восходящего узла для первой плоскости в градусах
                group.totalSatCount = group.satsPerPlane * group.planeCount;

                this.groups{length(this.groups) + 1} = group;                
                this.totalSatCount = this.totalSatCount + group.totalSatCount;
            end
        end

        function getInitialState(this)
            % Записывает массив кеплеровых элементов орбит для всех эшелонов

            this.state.elements = zeros(this.totalSatCount, 6);
            shift = 1;

            for group = this.groups
                for i = 1:length(group{1})
                    ending = shift + group{1}(i).totalSatCount - 1;
                    this.state.elements(shift:ending,:) = this.getInitialElements(group{1});
                    shift = ending + 1;
                end
            end
        end

        function elements = getInitialElements(this, group)
            % Расчитывает кеплеровы элементы орбит для каждого эшелона

            raanArray = linspace(group.startRaan, group.startRaan + group.maxRaan, group.planeCount + 1);
            raanArray = mod(raanArray(1:end-1), 2 * pi);

            elements = zeros(group.totalSatCount, 6);
            idx = 1;
            raanIDX = 0;
            for raan = raanArray
                for i = 0:group.satsPerPlane-1
                    sma = AstroConstants.earthRadius + group.altitude * 1000;
                    aol = 2 * pi / group.satsPerPlane * i + 2 * pi / group.totalSatCount * group.phaseShift * raanIDX;

                    elements(idx, :) = [sma, 0, 0, raan, group.inclination, aol];
                    idx = idx + 1;
                end
                raanIDX = raanIDX + 1;
            end
        end        

        function propagateJ2(this, epochArray)
            % Расчитывает eci координаты всех КА в указанные эпохи

            this.state.eci = zeros(this.totalSatCount, 3, length(epochArray));

            sma         = this.state.elements(:, 1);
            inclination = this.state.elements(:, 5);            
            raan0       = this.state.elements(:, 4);
            aol0        = this.state.elements(:, 6);

            raanPrecessionRate = -1.5 * (AstroConstants.earthJ2 * AstroConstants.earthGM^(1/2) * ...
                                    AstroConstants.earthRadius^2) ./ (sma.^(7/2)) .* cos(inclination);
            draconicOmega      = sqrt(AstroConstants.earthGM ./ sma.^3) .* ...
                                    (1 - 1.5 * AstroConstants.earthJ2 .* (AstroConstants.earthRadius ./ sma).^2) .* ...
                                    (1 - 4 .* cos(inclination).^2);

            for epochIdx = 1:length(epochArray)
                aol = aol0 + epochArray(epochIdx) * draconicOmega;
                raanOmega = raan0 + epochArray(epochIdx) * raanPrecessionRate;

                this.state.eci(:, :, epochIdx)  = [sma .* (cos(aol) .* cos(raanOmega) - sin(aol) .* cos(inclination) .* sin(raanOmega)), ...
                                                   sma .* (cos(aol) .* sin(raanOmega) + sin(aol) .* cos(inclination) .* cos(raanOmega)), ...
                                                   sma .* (sin(aol) .* sin(inclination))];
            end
        end
    end
end
