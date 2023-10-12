classdef Constellation < handle

    properties
        totalSatCount {mustBeNonnegative} = 0;
        groups = {};
        state = struct();
    end

    methods

        function this = Constellation(constellationName)
            arguments
                constellationName {mustBeNonzeroLengthText};
            end
            
            if ~exist('constants.mat', 'file')
                saveConstantToFile('constants.mat');
            end

            this.loadFromConfigFile(constellationName);
        end

        function loadFromConfigFile(this, constellationName)
            arguments
                this Constellation;
                constellationName {mustBeNonzeroLengthText};
            end
            
            databaseFileName = 'constellationsTest.json';
            databaseString   = fileread(databaseFileName);
            dataStructArray  = jsondecode(databaseString);
            thisData = [];

            for idx = 1:length(dataStructArray)
                if strcmpi(dataStructArray(idx).name, constellationName)
                    thisData = dataStructArray(idx);
                    break
                end
            end

            if isempty(thisData)
                error('Группировка не найдена в файле');
            end
            
            for idx = 1:size(thisData.Walkers, 1)
                group.inclination   = deg2rad(thisData.Walkers(idx, 1));    % Наклонение орбитальной плоскости [рад.]
                group.satsPerPlane  = thisData.Walkers(idx, 2);             % Число КА в каждой орбитальной плоскости группы.
                group.planeCount    = thisData.Walkers(idx, 3);             % Число орбитальных плоскостей в группе
                group.f             = thisData.Walkers(idx, 4);             % Фазовый сдвиг по аргументу широты между КА в соседних плоскостях.
                group.altitude      = thisData.Walkers(idx, 5);	            % Высота орбиты [км]
                group.maxRaan       = deg2rad(thisData.Walkers(idx, 6));    % Максимум прямого восхождения восходящего узла (при распределении орбитальных плоскостей) [рад.]
                group.startRaan     = deg2rad(thisData.Walkers(idx, 7));    % Прямое восхождение восходящего узла для первой плоскости [рад.]
                group.totalSatCount = group.satsPerPlane * group.planeCount;

                this.groups{length(this.groups) + 1} = group;                
                this.totalSatCount = this.totalSatCount + group.totalSatCount;
            end
        end

        function getInitialState(this)
            arguments
                this Constellation;
            end

            this.state(1).elements = zeros(this.totalSatCount, 6);
            shift = 1;

            for group = this.groups
                ending = shift + group{1}.totalSatCount - 1;
                this.state.elements(shift:ending,:) = this.getInitialElements(group{1});
                shift = ending + 1;
            end
        end

        function elementsArray = getInitialElements(this, group)
            arguments
                this Constellation;
                group {mustBeA(group, 'struct')};
            end
            
            load('constants.mat', 'earthRadius');

            raanArray = linspace(group.startRaan, group.startRaan + group.maxRaan, group.planeCount + 1);
            raanArray = mod(raanArray(1:end-1), 2 * pi);

            elementsArray = zeros(group.totalSatCount, 6);
            counter = 1;
            raanIdx = 0;
            for raan = raanArray
                for idx = 0:group.satsPerPlane-1
                    sma = earthRadius + group.altitude * 1000; % Большая полуось [м]
                    aol = 2 * pi / group.satsPerPlane * idx + 2 * pi / group.totalSatCount * group.f * raanIdx; % Аргумент широты

                    elementsArray(counter, :) = [sma, 0, 0, raan, group.inclination, aol];
                    counter = counter + 1;
                end
                raanIdx = raanIdx + 1;
            end
        end        

        function propagateJ2(this, epochArray)
            arguments
                this Constellation;
                epochArray {mustBeNonempty};
            end
            
            load('constants.mat', 'earthGM', 'earthRadius', 'earthJ2');

            this.state(1).eci = zeros(this.totalSatCount, 3, length(epochArray));

            sma         = this.state.elements(:, 1);
            inclination = this.state.elements(:, 5);            
            raan0       = this.state.elements(:, 4);
            aol0        = this.state.elements(:, 6);

            raanPrecessionRate = -1.5 .* (earthJ2 * earthGM^(1/2) * earthRadius^2) ./ (sma.^(7/2)) .* cos(inclination);
            draconicOmega      = sqrt(earthGM ./ sma.^3) .* (1 - 1.5 .* earthJ2 .* (earthRadius ./ sma).^2) .* (1 - 4 .* cos(inclination).^2);

            for epochIdx = 1:length(epochArray)
                aol = aol0 + epochArray(epochIdx) .* draconicOmega;
                raanOmega = raan0 + epochArray(epochIdx) .* raanPrecessionRate;
                
                % Координаты КА в инерциальной центрированной по Земле СК в
                % метрах.
                this.state.eci(:, :, epochIdx)  = [sma .* (cos(aol) .* cos(raanOmega) - sin(aol) .* cos(inclination) .* sin(raanOmega)), ...
                                                   sma .* (cos(aol) .* sin(raanOmega) + sin(aol) .* cos(inclination) .* cos(raanOmega)), ...
                                                   sma .* (sin(aol) .* sin(inclination))];
            end
        end
    end
end
