% BME 4015 Fall 2026
% Kids OR + Public Invention WAG Project
% Sarah Krishnamoorthy

%% PVC Scavenging System Model
% Uses Charlotte Pipe Schedule 40 PVC data.
% Compares PVC pipe sizes and predicts the operating point of Delta axial 
% fan and centrifugal blower

clear;
clc;
close all;

%% Gas properties
% Air approximation at ~20 deg C and 1 atm

rho = 1.204; % kg/m^3
mu = 1.81e-5; % Pa*s

%% OR information

siteName = input('OR Site Name: ', 's');

qMin = input('Lowest required scavenging flow rate (L/min): ');
qMax = input('Highest required scavenging flow rate (L/min): ');

if qMin <= 0
    error('Value must be greater than zero.');
elseif qMax <= qMin
    error('Value must be greater than the lowest flow rate.');
end

%% PVC geometry

straightLength = input('Total straight PVC pipe length (m): ');

no45 = input('Number of 45 degree elbows: ');
no90 = input('Number of 90 degree elbows: ');

%% Charlotte Pipe Schedule 40 PVC data

% Convert from in to mm.

% pipeIdMm = actual inside diameter in mm
pipeIdMm = [15.8 20.9 26.6 35.1 40.9 52.5 ...
    62.7 77.9 102.3 154.1 202.7];

% Straight-pipe length of one 45 degree elbow in m
l45M = [0.244 0.335 0.427 0.549 0.640 0.792 ...
    0.945 1.219 1.554 2.438 3.231];

% Straight-pipe length of one 90 degree elbow in m
l90M = [0.457 0.610 0.762 1.158 1.219 1.737 ...
    2.103 2.408 3.658 5.486 6.706];

fprintf('\nAvailable PVC inside diameters (mm):\n');
disp(pipeIdMm);

candidateSizes = input( ...
    'Enter PVC inside diameters to compare, e.g. [26.6 35.1 40.9]: ');

candidateSizes = unique(candidateSizes);

%% Check PVC input

if any(~ismember(candidateSizes, pipeIdMm))
    error('One or more selected diameters are not available.');
end

%% Fan performance data

% Delta ASB0412VHA-AF00 axial fan
% Approximate points read from manufacturer 12 V P-Q curve, needed for K

% Flow rate = L/min
axialFlow = [0 50 100 150 200 225 245 269];

% Pressure = mmH2O
axialPressureMmH2O = [6.422 5.0 3.6 2.35 1.7 1.1 0.5 0];

% Delta BFB0512HHA-CZMQ centrifugal blower
% Approximate points read from manufacturer 12 V P-Q curve

blowerFlow = [0 20 40 60 80 90 100 113];

blowerPressureMmH2O = [20.09 18.0 15.5 12.3 9.0 6.2 3.0 0];

%% Convert fan pressure to Pa

mmH2OToPa = 9.80665;

axialPressure = axialPressureMmH2O .* mmH2OToPa;
blowerPressure = blowerPressureMmH2O .* mmH2OToPa;

%% Model flow range

modelMaxFlow = max([qMax axialFlow blowerFlow]);

qModel = linspace(0.1, modelMaxFlow, 5000); % L/min
qModelM3s = qModel ./ (1000 * 60); % m^3/s

%% Fan curves across model flow range

axialPressureModel = interp1(axialFlow, axialPressure, qModel, 'linear', NaN);

blowerPressureModel = interp1(blowerFlow, blowerPressure, qModel, 'linear', NaN);

%% Prepare result storage

numberOfSizes = length(candidateSizes);

pressureAtMax = zeros(numberOfSizes, 1);

axialPredictedFlow = zeros(numberOfSizes, 1);
axialOperatingPressure = zeros(numberOfSizes, 1);
axialStatus = strings(numberOfSizes, 1);

blowerPredictedFlow = zeros(numberOfSizes, 1);
blowerOperatingPressure = zeros(numberOfSizes, 1);
blowerStatus = strings(numberOfSizes, 1);

systemPressureMatrix = zeros(numberOfSizes, length(qModel));

%% Evaluate each PVC diameter

for sizeNumber = 1:numberOfSizes

    diameterMm = candidateSizes(sizeNumber);
    dataIndex = find(pipeIdMm == diameterMm);

    diameter = diameterMm / 1000; % mm to m
    pipeArea = pi * diameter^2 / 4;

    %% Equivalent bend length

    l45 = l45M(dataIndex);
    l90 = l90M(dataIndex);

    totalEquivalentLength = no45 * l45 + no90 * l90;

    %% Gas velocity

    velocity = qModelM3s ./ pipeArea;

    %% Reynolds number

    reynoldsNumber = rho .* velocity .* diameter ./ mu;

    %% Darcy friction factor

    frictionFactor = zeros(size(reynoldsNumber));

    for i = 1:length(reynoldsNumber)

        if reynoldsNumber(i) < 2300

            % Laminar flow

            frictionFactor(i) = 64 / reynoldsNumber(i);

        elseif reynoldsNumber(i) <= 4000

            % Transitional flow
            %
            % Interpolate between laminar friction at Re = 2300
            % and Blasius turbulent friction at Re = 4000.
            % This prevents an artificial jump in pressure loss.

            fLaminar = 64 / 2300;
            fTurbulent = 0.3164 / (4000^0.25);

            fraction = (reynoldsNumber(i) - 2300) / (4000 - 2300);

            frictionFactor(i) = fLaminar + fraction * (fTurbulent - fLaminar);

        else

            % Smooth-pipe turbulent approximation

            frictionFactor(i) = 0.3164 / (reynoldsNumber(i)^0.25);

        end

    end

    %% Straight PVC pressure loss

    pipePressure = frictionFactor .* (straightLength / diameter) .* ...
        (rho .* velocity.^2 ./ 2);

    %% Elbow pressure loss

    bendPressure = frictionFactor .* (totalEquivalentLength / diameter) .* ...
        (rho .* velocity.^2 ./ 2);

    %% Total system pressure loss

    systemPressure = pipePressure + bendPressure;

    systemPressureMatrix(sizeNumber, :) = systemPressure;

    pressureAtMax(sizeNumber) = interp1(qModel, systemPressure, qMax);

    %% Axial fan operating point

    axialDifference = abs(axialPressureModel - systemPressure);

    axialDifference(isnan(axialDifference)) = Inf;

    [~, axialIndex] = min(axialDifference);

    axialPredictedFlow(sizeNumber) = qModel(axialIndex);
    axialOperatingPressure(sizeNumber) = systemPressure(axialIndex);

    %% Axial fan status

    if axialPredictedFlow(sizeNumber) < qMin
        axialStatus(sizeNumber) = "FAIL";

    elseif axialPredictedFlow(sizeNumber) <= qMax
        axialStatus(sizeNumber) = "PASS";

    else
        axialStatus(sizeNumber) = "ABOVE TARGET";
    end

    %% Centrifugal blower operating point

    blowerDifference = ...
        abs(blowerPressureModel - systemPressure);

    blowerDifference(isnan(blowerDifference)) = Inf;

    [~, blowerIndex] = min(blowerDifference);

    blowerPredictedFlow(sizeNumber) = qModel(blowerIndex);
    blowerOperatingPressure(sizeNumber) = systemPressure(blowerIndex);

    %% Centrifugal blower status

    if blowerPredictedFlow(sizeNumber) < qMin
        blowerStatus(sizeNumber) = "FAIL";

    elseif blowerPredictedFlow(sizeNumber) <= qMax
        blowerStatus(sizeNumber) = "PASS";

    else
        blowerStatus(sizeNumber) = "ABOVE TARGET";
    end

end

%% Results Display

fprintf('\nRESULTS - %s\n', siteName);
fprintf('Required flow: %.1f - %.1f L/min\n', qMin, qMax);
fprintf('----------------------------------------\n');

for sizeNumber = 1:numberOfSizes

    fprintf('\nPVC inside diameter: %.1f mm\n', ...
        candidateSizes(sizeNumber));

    fprintf('Pressure loss at %.1f L/min: %.2f Pa\n', ...
        qMax, pressureAtMax(sizeNumber));

    fprintf('Axial fan: %.1f L/min - %s\n', ...
        axialPredictedFlow(sizeNumber), ...
        axialStatus(sizeNumber));

    fprintf('Centrifugal blower: %.1f L/min - %s\n', ...
        blowerPredictedFlow(sizeNumber), ...
        blowerStatus(sizeNumber));

end

%% Lowest hydraulic resistance

if numberOfSizes > 1

    [lowestPressure, lowestIndex] = min(pressureAtMax);

    fprintf('\nLowest resistance of tested sizes:\n');

    fprintf('PVC ID = %.1f mm at %.2f Pa\n', ...
        candidateSizes(lowestIndex), ...
        lowestPressure);

end

%% Plot PVC resistance over required flow range

figure;

targetRange = qModel >= qMin & qModel <= qMax;

hold on;

legendLabels = strings(numberOfSizes, 1);

for sizeNumber = 1:numberOfSizes

    plot(qModel(targetRange), systemPressureMatrix(sizeNumber, ...
        targetRange),':');

    legendLabels(sizeNumber) = sprintf('PVC resistance, ID = %.1f mm', ...
        candidateSizes(sizeNumber));

end

grid on;

xlabel('Flow Rate (L/min)');
ylabel('Pressure Loss (Pa)');

title('PVC Pressure Loss Comparison');

legend(legendLabels, 'Location', 'northwest');

hold off;

%% Plot fan curves and system resistance

figure;

hold on;

fanLegend = strings(numberOfSizes + 2, 1);

%% PVC system curves

for sizeNumber = 1:numberOfSizes

    plot(qModel, systemPressureMatrix(sizeNumber, :), ':');

    fanLegend(sizeNumber) = sprintf('PVC resistance, ID = %.1f mm', ...
        candidateSizes(sizeNumber));

end

%% Individual fan curves

plot(qModel, axialPressureModel);
plot(qModel, blowerPressureModel);

fanLegend(numberOfSizes + 1) = "Delta axial fan";
fanLegend(numberOfSizes + 2) = "Delta centrifugal blower";

%% Plot operating point markers

for sizeNumber = 1:numberOfSizes

    plot(axialPredictedFlow(sizeNumber), axialOperatingPressure(sizeNumber), ...
        'o');

    plot(blowerPredictedFlow(sizeNumber), blowerOperatingPressure(sizeNumber), ...
        's');

end

%% Required flow limits

xline(qMin, '--', 'Minimum Required Flow');
xline(qMax, '--', 'Maximum Required Flow');

grid on;

xlabel('Flow Rate (L/min)');
ylabel('Pressure (Pa)');

title('System Resistance and Fan Performance');

legend(fanLegend, 'Location', 'northeast');

% Turn on MATLAB data tips
datacursormode on;

hold off;