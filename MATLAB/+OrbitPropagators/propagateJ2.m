function propagateJ2(constellation, epochs)
    constellation.state.eci = zeros(constellation.totalSatCount, 3, length(epochs));
    
    sma         = constellation.state.elements(:, 1);
    inclination = constellation.state.elements(:, 5);
    raan0       = constellation.state.elements(:, 4);
    aol0        = constellation.state.elements(:, 6);
    
    earthJ2     = Constants.AstroConstants.earthJ2;
    earthGM     = Constants.AstroConstants.earthGM;
    earthRadius = Constants.AstroConstants.earthRadius;
    
    raanPrecessionRate = -1.5 * (earthJ2 * earthGM^(1/2) * earthRadius^2) ./ (sma.^(7/2)) .* cos(inclination);
    draconicOmega = sqrt(earthGM ./ sma.^3) .* (1 - 1.5 * earthJ2 .* (earthRadius ./ sma).^2) ...
                    .* (1 - 4 .* cos(inclination).^2);
    
    for epochIdx = 1:length(epochs)
        aol = aol0 + epochs(epochIdx) * draconicOmega;
        raanOmega = raan0 + epochs(epochIdx) * raanPrecessionRate;
    
        constellation.state.eci(:, :, epochIdx)  = ...
            [sma .* (cos(aol) .* cos(raanOmega) - sin(aol) .* cos(inclination) .* sin(raanOmega)), ...
             sma .* (cos(aol) .* sin(raanOmega) + sin(aol) .* cos(inclination) .* cos(raanOmega)), ...
             sma .* (sin(aol) .* sin(inclination))];
    end
end