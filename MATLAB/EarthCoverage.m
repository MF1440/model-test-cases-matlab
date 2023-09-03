classdef EarthCoverage < handle
    % класс для проверки покрытия поверхности земли спутниками 

    properties
        meshDensityArray = [20, 50, 100, 200, 500, 1000];  % плотности сетки на которых будет считаться покрытие
        alphaAngle;                                        % заданный угол обзора КА
        gammaAngle;                                        % угол покрытия части поверхности одним спутником 
        satCoords;
        mesh;
    end

    properties (Dependent)
        coverageMatrix
    end

    methods

        function this = EarthCoverage()
            this.alphaAngle = deg2rad(AstroConstants.alphaAngle);
            this.gammaAngle = asin((AstroConstants.earthRadius + AstroConstants.altitude * 1000) / ...
                AstroConstants.earthRadius * sin(this.alphaAngle)) - this.alphaAngle;  
        end

        function flag = checkEarthCoverage(this, satCountPerOrbit, orbitCount, phase, epoch)
            % Проверяет полное покрытие земли эшелоном КА
            %
            % Для проверки покрытия используются сетки разного размера.
            % Заданная группировка проверяется сначала на крупной сетке и в
            % случае полного покрытия проверяется на более мелких сетках.
            % 
            % Возвращает true в случае полного покрытия и false в ином случае  

            this.getSatSphereCoords(satCountPerOrbit, orbitCount, phase, epoch);
            for meshDensity=this.meshDensityArray 
                this.mesh.phi = linspace(0, pi, 2*meshDensity);  % узлы сетки по сферической координате phi (меняется от 0 до pi)
                this.mesh.theta = linspace(0, deg2rad(AstroConstants.inclination), meshDensity);  % узлы сетки по сферической координате theta (меняется от 0 до pi/3)
                if min(this.coverageMatrix, [], "all") > 0  % проверяет найдено ли полное покрытие 
                    disp(['На сетке ', num2str(2*meshDensity),'x', num2str(meshDensity)...
                        ' найдена конфигурация с полным покрытием: '])
                    disp(['орбит ', int2str(orbitCount), ' спутников ', ...
                        int2str(satCountPerOrbit), ' фазовый сдвиг ', int2str(phase), ...
                    ' спутников (всего ', int2str(satCountPerOrbit * orbitCount),')'])
                    flag = true;
                else  % если полного покрытия нет, то выходит из цикла без проверки на более мелких сетках
                    flag = false;
                    break
                end 
            end  % конец цикла по плотности сетки
        end

        function getSatSphereCoords(this, satCountPerOrbit, orbitCount, phase, epoch)
            % Получает координаты КА в сферической системе координат

            config.Walkers = [[AstroConstants.inclination, ...
                            satCountPerOrbit, ...
                            orbitCount, ...
                            phase, ...
                            AstroConstants.altitude, ...
                            AstroConstants.maxRaan, ...
                            AstroConstants.startRaan]];
            constellation = Constellation(config);
            constellation.getInitialState();
            constellation.propagateJ2([epoch]);
            satEciCoords = constellation.state.eci(:, :, [1]);
            [this.satCoords.phi, this.satCoords.theta, ~] = cart2sph(satEciCoords(:,1), ...
                                                                    satEciCoords(:,2), ...
                                                                    satEciCoords(:,3));
        end

        function matrix = get.coverageMatrix(this)
            % Считает матрицу покрытия, каждый элемент которой равен
            % количеству КА покрывающих узел сетки

            matrix = zeros(length(this.mesh.phi), length(this.mesh.theta), "int8");
            for thetaIdx=1:length(this.mesh.theta)
                for phiIdx=1:length(this.mesh.phi)
                    coverageDensity = this.checkPointCover(this.mesh.phi(phiIdx), ...
                                                        this.mesh.theta(thetaIdx));
                    matrix(phiIdx, thetaIdx) = coverageDensity;
                end  % конец цикла по phi
            end  % конец цикла по theta
        end

        function coverageDensity = checkPointCover(this, phi, theta)
            % Проверяет сколько КА покрывает текущую точку (phi, theta)
            % 
            % Считает количество КА, угол которых с точкой (phi, theta) 
            % и вершиной в центре земли меньше угла gamma

            angleBetweenSatAndPoint = this.calcAngleBetweenSatAndPoint(phi, theta);
            coverageDensity = sum((angleBetweenSatAndPoint <= this.gammaAngle));
        end

        function angleBetweenSatAndPoint = calcAngleBetweenSatAndPoint(this, phi, theta)
            % Считает угол между КА и заданной точкой на сфере (phi, theta) с вершиной в центре земли
            
            angleBetweenSatAndPoint = acos(cos(this.satCoords.theta) .* cos(this.satCoords.phi) .* ...
                                            cos(theta) .* cos(phi) + ... 
                                            + cos(this.satCoords.theta) .* sin(this.satCoords.phi) .* ...
                                            cos(theta) .* sin(phi) + ...
                                            sin(this.satCoords.theta) .* sin(theta));
        end

    end
end