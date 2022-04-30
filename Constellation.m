classdef Constellation < handle
% Класс Constellation позволяет создать объект, описывающий спутниковую
% группировку, на основе параметров, считанных из конфигурационного файла,
% определить положения КА в первоначальный момент времени и расчитать
% дальнейшие положения КА в заданные моменты времени.
% Время задается в эпохах.

    properties
        totalSatCount = 0; % Общее количество спутников текущей орбитальной конфигурации
        groups = {}; % Массиы ячеек для хранящая параметров спутниковых группировок 
        state; % Структура для хранения результатов расчета 

        % константы
        earthRadius = 6378135;           % Экваториальный радиус Земли [m]
        earthGM = 3.986004415e+14;       % Гравитационный параметр Земли [m3/s2]
        earthJ2 = 1.082626e-3;           % Коэффициент первой зональной гармоники в расширении гравитационного поля Земли
    end

    methods

        function this = Constellation(varargin)
            % Коснструктор класса.
            % varargin принимает единственный параметр - наименование
            % спутниковой группироки, по которому осуществляется поиск
            %в конфигурационном файле constellationsTest.json
            if isempty(varargin)
                msg = 'Не указан обязательный параметр - наименование группировки.';
                error(msg)
            end
            this.loadFromConfigFile(varargin{1});
        end

        function loadFromConfigFile(this, code)
            % Метод, используемый конструктором класса для инициализации
            % параметров спутниковой группировки из конфигурационного
            % файла.
            fileName = 'constellationsTest.json';
            str = fileread(fileName);
            data = jsondecode(str);
            dataThis = [];

            for i = 1:length(data)
                if strcmpi(data(i).name, code)
                    dataThis = data(i);
                    break
                end
            end

            if isempty(dataThis)
                disp('Группировка не найдена в файле');
                return
            end

            for i = 1:size(dataThis.Walkers, 1)
                group.inclination   = deg2rad(dataThis.Walkers(i, 1)); % Наклонение орбит
                group.satsPerPlane  = dataThis.Walkers(i, 2); % Количество спутников в плоскости
                group.planeCount    = dataThis.Walkers(i, 3); % Количество плоскостей в группировке
                group.f             = dataThis.Walkers(i, 4);
                group.altitude      = dataThis.Walkers(i, 5); % Высота орбиты
                group.maxRaan       = deg2rad(dataThis.Walkers(i, 6)); % Максимальное значение градуса дуги окружности,
                                                                       % вдоль которой распределены плоскости орбит 
                group.startRaan     = deg2rad(dataThis.Walkers(i, 7)); % Начальное значение 
                group.totalSatCount = group.satsPerPlane * group.planeCount; % Общее количество спутников в группировке

                this.groups{length(this.groups) + 1} = group;                
                this.totalSatCount = this.totalSatCount + group.totalSatCount;
            end
        end

        function getInitialState(this)
            % Функция, формирующая структуру элементов орбит КА
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
            % Функция вычисляет элементы орбиты для всех КА в начальный
            % момент вемени
            
            % Вычисление положения орбитальных плоскостей
            raanArray = linspace(group.startRaan, group.startRaan + group.maxRaan, group.planeCount + 1);
            raanArray = mod(raanArray(1:end-1), 2 * pi);

            elements = zeros(group.totalSatCount, 6);
            idx = 1;
            raanIdx = 0;
            
            sma = this.earthRadius + group.altitude * 1000; %  *1000 - перевод высоты орбиты из км в м
            
            for raan = raanArray
                for i = 0:group.satsPerPlane-1
                    
                    % Угол дуги , на котором располагается каждый КА в отдельной орбитальной плоскости
                    % с учетом смещения КА в следующей плоскости
                    % относительно КА в предыдущей
                    aol = 2 * pi / group.satsPerPlane * i + 2 * pi / group.totalSatCount * group.f * raanIdx;

                    elements(idx, :) = [sma, 0, 0, raan, group.inclination, aol];
                    idx = idx + 1;
                end
                raanIdx = raanIdx + 1;
            end
        end        

        function propagateJ2(this, epochs)
            % Функция, расчитывающая положение всех КА в заданные моменты
            % времени
            this.state.eci = zeros(this.totalSatCount, 3, length(epochs));

            inclination = this.state.elements(:, 5);
            sma         = this.state.elements(:, 1);
            Omega0      = sqrt(this.earthGM ./ sma.^3);
            aol0        = this.state.elements(:, 6);

            raanPrecessionRate = -1.5 * (this.earthJ2 * this.earthGM^(1/2) * this.earthRadius^2) ./...
                                 (sma.^(7/2)) .* cos(inclination);
            draconicOmega      = sqrt(this.earthGM ./ sma.^3) .* (1 - 1.5 * this.earthJ2 .* (this.earthRadius ./ sma).^2) .* ...
                                 (1 - 4 .* cos(inclination).^2);

            for epochIdx = 1:length(epochs)
                aol = aol0 + epochs(epochIdx) * draconicOmega;
                Omega = Omega0 + epochs(epochIdx) * raanPrecessionRate;

                this.state.eci(:, :, epochIdx)  = [sma .* (cos(aol) .* cos(Omega) - sin(aol) .* cos(inclination) .* sin(Omega)), ...
                                                   sma .* (cos(aol) .* sin(Omega) + sin(aol) .* cos(inclination) .* cos(Omega)), ...
                                                   sma .* (sin(aol) .* sin(inclination))];
            end
        end
    end
end