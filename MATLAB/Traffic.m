classdef Traffic < handle

    properties
        trafficDistribution = {};
        satState;
        userState;

        % константы
        earthRadius = 6378135;           % Экваториальный радиус Земли [m]
    end

    methods
        function objTraffic = Traffic(satData, userDataFileName)
            if isempty(satData)
                return
            end
            objTraffic.loadUserData(userDataFileName);
            objTraffic.satDataOnEarth(satData);
        end

        function loadUserData(this, fileName)
            userData = load(fileName);
            
            if isempty(userData)
                disp('Данные пользователей не найдены');
                return
            end

            this.userState.userLocation = userData.coordsEcef;
            this.userState.userTraffic = userData.values;
        end

        function satDataOnEarth(this, satData)

            if isempty(satData)
                disp('Данные КА не найдены');
                return
            end

            this.satState.satOnEarth = zeros(size(satData));

            tetaAngle = atan(sqrt(satData(:, 1, :) .^2 + satData(:, 2, :) .^2) ./ satData(:, 3, :));
            phiAngle = atan(satData(:, 2, :) ./ satData(:, 1, :));

            this.satState.satOnEarth = [this.earthRadius .* sin(tetaAngle) .* cos(phiAngle), ...
                                        this.earthRadius .* sin(tetaAngle) .* sin(phiAngle), ...
                                        this.earthRadius .* cos(tetaAngle)];

        end

        function calcSatTraffic(this)

            nearestSat = zeros(length(this.userState.userLocation(:, 1)), length(this.satState.satOnEarth(1, 1, :)));
            satTraffic = zeros(length(this.satState.satOnEarth(:, 1, 1)), length(this.satState.satOnEarth(1, 1, :)));

            for epochIdx = 1:length(this.satState.satOnEarth(1, 1, :))
                for userIdx = 1:length(this.userState.userLocation(:, 1))                   
                    distance = sqrt((this.userState.userLocation(userIdx, 1) - this.satState.satOnEarth(:, 1, epochIdx)) .^2 + ...
                                    (this.userState.userLocation(userIdx, 2) - this.satState.satOnEarth(:, 2, epochIdx)) .^2 + ...
                                    (this.userState.userLocation(userIdx, 3) - this.satState.satOnEarth(:, 3, epochIdx)) .^2);
                    [~, satNum] = min(distance);

                    nearestSat(userIdx, epochIdx) = satNum;
                    satTraffic(satNum, epochIdx) = satTraffic(satNum, epochIdx) + this.userState.userTraffic(userIdx);

                    if mod(userIdx, 500000) == 0
                        disp(['Посчитано ' num2str(userIdx) ' пользователей']);
                        toc
                    end
                end % Конец цикла по пользователям
                disp(['Эпоха-' num2str(epochIdx) ' посчитана']);
            end % Конец цикла по эпохам

            this.trafficDistribution.satDistribution = nearestSat;
            this.trafficDistribution.trafficForSatDistribution = satTraffic;
        end
    end
end
