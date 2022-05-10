
function visibles = findGatewayZrv(constellation, epochIdx)
    earthRadius = 6378135; % Экваториальный радиус Земли [m]

    gateways = jsondecode(fileread('gatewaysTest.json'));
    distances = zeros(constellation.totalSatCount, length(gateways));
    gatewayPoints = zeros(length(gateways), 3);

    ranges = zeros(length(constellation.groupArray), 1);
    satellites = constellation.state.eci(:, :, epochIdx);

    for groupIdx = 1: length(constellation.groupArray)
        ranges(groupIdx) = 1000 * constellation.groupArray{groupIdx}.altitude / cosd(25);
    end

    for gatewayIdx = 1: length(gateways)
        gatewayPoints(gatewayIdx, 1) = earthRadius ...
                                     * cosd(gateways(gatewayIdx).lat) ...
                                     * cosd(gateways(gatewayIdx).lon);
        gatewayPoints(gatewayIdx, 2) = earthRadius ...
                                     * cosd(gateways(gatewayIdx).lat) ...
                                     * sind(gateways(gatewayIdx).lon);
        gatewayPoints(gatewayIdx, 3) = earthRadius ...
                                     * sind(gateways(gatewayIdx).lat);
    end

    for satelliteIdx = 1: length(satellites)
        for gatewayIdx = 1: length(gatewayPoints)
            distances(satelliteIdx, gatewayIdx) = sqrt((satellites(satelliteIdx, 1) - gatewayPoints(gatewayIdx, 1)) ^ 2 ...
                                                 + (satellites(satelliteIdx, 2) - gatewayPoints(gatewayIdx, 2)) ^ 2 ...
                                                 + (satellites(satelliteIdx, 3) - gatewayPoints(gatewayIdx, 3)) ^ 2);
        end
    end

    shift = 1;
    visibles = [];

    for groupIdx = 1: length(constellation.groupArray)
        ending = shift + constellation.groupArray{groupIdx}.totalSatCount - 1;
        head = (shift - 1) + unique(ceil(find(distances(shift: ending, :) ...
            <= ranges(groupIdx)) / length(gateways)));

        visibles = cat(1, visibles, head);
        shift = ending + 1;
    end
end

