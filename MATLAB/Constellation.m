classdef Constellation < handle

    properties
        totalSatCount = 0;
        groupList = {}; 
        state = struct;
        
        % Константы
        earthRadius = 6378135;           % Экваториальный радиус Земли [m]
        earthGM = 3.986004415e+14;       % Гравитационный параметр Земли [m3/s2]
        earthJ2 = 1.082626e-3;           % Вторая зональная гармоника геопотенциала
    end

    methods

        function this = Constellation(varargin)
            % varargin - ячейка с названием группировки
            if isempty(varargin)
                return
            end
            this.loadFromConfigFile(varargin{1});
        end

        function loadFromConfigFile(this, code)
            fileName = 'constellationsTest.json';
            str = fileread(fileName);
            data = jsondecode(str);                      
            dataThis = {};

            for constellationIdx = 1:length(data)  % Ищем нужную группировку по имени
                if strcmpi(data(constellationIdx).name, code)
                    dataThis = data(constellationIdx);   
                    break
                end
            end

            if isempty(dataThis)
                disp('Группировка не найдена в файле');
                return
            end

            for groupIdx = 1:size(dataThis.Walkers, 1)
                group.inclination = deg2rad(dataThis.Walkers(groupIdx, 1));        % наклонение орбитальной плоскости
                group.satsPerPlane = dataThis.Walkers(groupIdx, 2);				   % число КА в каждой орбитальной плоскости группы
                group.planeCount = dataThis.Walkers(groupIdx, 3);				   % число орбитальных плоскостей в группе
                group.phaseShift = dataThis.Walkers(groupIdx, 4);				   % фазовый сдвиг по аргументу широты между КА в соседних плоскостях
                group.altitude = dataThis.Walkers(groupIdx, 5);					   % высота орбиты
                group.maxRaan = deg2rad(dataThis.Walkers(groupIdx, 6));            % максимум прямого восхождения восходящего узла (при распределении орбитальных плоскостей)
                group.startRaan = deg2rad(dataThis.Walkers(groupIdx, 7));		   % прямое восхождение восходящего узла для первой плоскости
                group.totalSatCount = group.satsPerPlane * group.planeCount;

                this.groupList{length(this.groupList) + 1} = group;                
                this.totalSatCount = this.totalSatCount + group.totalSatCount;
            end
        end

        function getInitialState(this)
            this.state.elements = zeros(this.totalSatCount, 6);
            shift = 1;

            for group = this.groupList
                for i = 1:length(group{1})
                    ending = shift + group{1}(i).totalSatCount - 1;
                    this.state.elements(shift:ending,:) = this.getInitialElements(group{1});
                    shift = ending + 1;
                end        % Конец цикла по спутникам внутри группы
            end            % Конец цикла по группам
        end                % Конец описания метода

        function elements = getInitialElements(this, group)
            raanList = linspace(group.startRaan, group.startRaan + group.maxRaan, group.planeCount + 1);
            raanList = mod(raanList(1:end-1), 2 * pi);

            elements = zeros(group.totalSatCount, 6);
            idx = 1;
            raanIdx = 0;
            sma = this.earthRadius + group.altitude * 1000;
                        
            for raan = raanList
                for idxSat = 0:group.satsPerPlane-1                    
                    aol = 2 * pi / group.satsPerPlane * idxSat + 2 * pi / group.totalSatCount * group.phaseShift * raanIdx;
                    elements(idx, :) = [sma, 0, 0, raan, group.inclination, aol];
                    idx = idx + 1;
                end                    % Конец цикла по КА в орбитальной плоскости группы 
                raanIdx = raanIdx + 1;
            end                        % Конец цикла по орбитальным плоскостям 
        end                            % Конец описания функции getInitialElements          

        function propagateJ2(this, epochList)
            this.state.eci = zeros(this.totalSatCount, 3, length(epochList));

            sma         = this.state.elements(:, 1);
            inclination = this.state.elements(:, 5);            
            raan0       = this.state.elements(:, 4);
            aol0        = this.state.elements(:, 6);

            raanPrecessionRate = -1.5 * (this.earthJ2 * this.earthGM^(1/2) * this.earthRadius^2) ./ (sma.^(7/2)) .* cos(inclination);
            draconicOmega      = sqrt(this.earthGM ./ sma.^3) .* (1 - 1.5 * this.earthJ2 .* (this.earthRadius ./ sma).^2) .* (1 - 4 .* cos(inclination).^2);
            
            aol = repmat(aol0, 1, length(epochList)) + repmat(epochList, length(draconicOmega), 1) .* repmat(draconicOmega, 1, length(epochList));
            raanOmega = repmat(raan0, 1, length(epochList)) + repmat(epochList, length(raanPrecessionRate), 1) .* repmat(raanPrecessionRate, 1, length(epochList));
            
            for epochIdx = 1:length(epochList)
                
                this.state.eci(:, :, epochIdx)  = [sma .* (cos(aol(:,epochIdx)) .* cos(raanOmega(:,epochIdx)) - sin(aol(:,epochIdx)) .* cos(inclination) .* sin(raanOmega(:,epochIdx))), ...
                                                   sma .* (cos(aol(:,epochIdx)) .* sin(raanOmega(:,epochIdx)) + sin(aol(:,epochIdx)) .* cos(inclination) .* cos(raanOmega(:,epochIdx))), ...
                                                   sma .* (sin(aol(:,epochIdx)) .* sin(inclination))];
            end    % Конец цикла по эпохам 
            
        end        % Конец функции propagateJ2
    end            % Конец секции с описанием методов класса
end
