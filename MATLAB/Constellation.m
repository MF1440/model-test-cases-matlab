classdef Constellation < handle

    properties
        totalSatCount = 0;
        groupArray = {};
        state;

        % константы
        earthRadius = 6378135;           % Экваториальный радиус Земли [m]
        earthGM = 3.986004415e+14;       % Гравитационный параметр Земли [m3/s2]
        earthJ2 = 1.082626e-3;           % Вторая зональная гармоника геопотенциала
    end

    methods

        function objConstel = Constellation(varargin)
            if isempty(varargin)
                return
            end
            objConstel.loadFromConfigFile(varargin{1});
        end

        function loadFromConfigFile(this, code)
            fileName = 'constellationsTest.json';
            str = fileread(fileName);
            data = jsondecode(str);
            dataThis = [];

            for dataIdx = data
                if strcmpi(dataIdx.name, code)
                    dataThis = dataIdx;
                    break
                end
            end

            if isempty(dataThis)
                disp('Группировка не найдена в файле');
                return
            end

            for idx = 1:size(dataThis.Walkers, 1)
                this.groupArray{idx}.inclination = deg2rad(dataThis.Walkers(idx, 1));            % наклонение орбитальной плоскости
                this.groupArray{idx}.satPerPlaneCount = dataThis.Walkers(idx, 2);                % число КА в каждой орбитальной плоскости группы
                this.groupArray{idx}.planeCount = dataThis.Walkers(idx, 3);                      % число орбитальных плоскостей в группе
                this.groupArray{idx}.phaseShift = dataThis.Walkers(idx, 4);                      % фазовый сдвиг по аргументу широты между КА в соседних плоскостях
                this.groupArray{idx}.altitudeKilometers = dataThis.Walkers(idx, 5);              % высота орбиты над поверхностью земли [km]
                this.groupArray{idx}.maxRaan = deg2rad(dataThis.Walkers(idx, 6));                % максимум прямого восхождения восходящего узла (при распределении орбитальных плоскостей)
                this.groupArray{idx}.startRaan = deg2rad(dataThis.Walkers(idx, 7));              % прямое восхождение восходящего узла для первой плоскости
                this.groupArray{idx}.totalSatCount = ...
                    this.groupArray{idx}.satPerPlaneCount * this.groupArray{idx}.planeCount;     % общее число КА в группе
                
                this.totalSatCount = this.totalSatCount + this.groupArray{idx}.totalSatCount;
            end
        end

        function calcInitialState(this)
            this.state.elements = zeros(this.totalSatCount, 4);
            shift = 1;

            for groupIdx = this.groupArray
                ending = shift + groupIdx{1}.totalSatCount - 1;
                this.state.elements(shift:ending,:) = this.calcInitialElements(groupIdx{1});
                shift = ending + 1;
            end
        end

        function elements = calcInitialElements(this, group)
            raanArray = linspace(group.startRaan, group.startRaan + group.maxRaan, group.planeCount + 1);
            raanArray = mod(raanArray(1:end-1), 2 * pi);

            elements = zeros(group.totalSatCount, 4);
            satNumber = 1;
            for planeIdx = 1:group.planeCount
                for satIdx = 1:group.satPerPlaneCount
                    planeRadius = this.earthRadius + group.altitudeKilometers * 1000;
                    directionAngle = 2 * pi / group.satPerPlaneCount * (satIdx - 1) + 2 * pi / group.totalSatCount * group.phaseShift * (planeIdx - 1);

                    elements(satNumber, :) = [planeRadius, raanArray(planeIdx), group.inclination, directionAngle];
                    satNumber = satNumber + 1;
                end % конец цикла по КА
            end % конец цикла по орбитам
        end        

        function propagateJ2(this, epochArray)
            this.state.eci = zeros(this.totalSatCount, 3, length(epochArray));

            planeRadius              = this.state.elements(:, 1);
            inclination              = this.state.elements(:, 3);            
            raanZeroEpoch            = this.state.elements(:, 2);
            directionAngleZeroEpoch  = this.state.elements(:, 4);

            raanPrecessionRate = -1.5 * (this.earthJ2 * this.earthGM^(1/2) * this.earthRadius^2) ./ (planeRadius.^(7/2)) .* cos(inclination);
            draconicOmega      = sqrt(this.earthGM ./ planeRadius.^3) .* (1 - 1.5 * this.earthJ2 .* (this.earthRadius ./ planeRadius).^2) .* (1 - 4 .* cos(inclination).^2);

            for epochIdx = 1:length(epochArray)
                directionAngle = directionAngleZeroEpoch + epochArray(epochIdx) * draconicOmega;
                raanOmega = raanZeroEpoch + epochArray(epochIdx) * raanPrecessionRate;

                this.state.xyzCoordinate(:, :, epochIdx)  = [planeRadius .* (cos(directionAngle) .* cos(raanOmega) - sin(directionAngle) .* cos(inclination) .* sin(raanOmega)), ...
                                                             planeRadius .* (cos(directionAngle) .* sin(raanOmega) + sin(directionAngle) .* cos(inclination) .* cos(raanOmega)), ...
                                                             planeRadius .* (sin(directionAngle) .* sin(inclination))];
            end
        end
    end
end
