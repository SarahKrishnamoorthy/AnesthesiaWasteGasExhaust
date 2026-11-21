% BME 4015 Fall 2026
% Kids OR + Public Invention WAG Project
% Code written by Sarah Krishnamoorthy

%% TO COMPARE AXIAL VS CENTRIFUGAL FAN, HOW MUCH RESISTANCE IN PIPE

clear; clc; close all;

%% Initialize gas properties
% WAG is approximated to be air at ~20 deg C and 1 atm.

rho = 1.204; % kg/m^3, air density
mu = 1.81e-5; % Pa*s, dynamic viscosity of air

%% Piping properties and geometry

% Using WATO EX-65 operating range: high flow, low vacuum AGSS
Qmin = 50; % L/min
Qmax = 80; % L/min

%% User input

inputD = input('PVC pipe inner diameter (mm): ');
L = input('PVC pipe length (m): ');
%% Bend information

numBends = input('Number of bends in PVC pipe: ');

bendType = input(['Bend type:\n' ...
    '1 = 45 degree bend\n' ...
    '2 = Long-radius 90 degree bend\n' ...
    '3 = Standard 90 degree elbow\n' ...
    'Enter selection: ']);

if bendType == 1
    Kbend = 0.35;

elseif bendType == 2
    Kbend = 0.30;

elseif bendType == 3
    Kbend = 0.90;

else
    error('Invalid bend type selected.');
end
Kbend = input('Loss coefficient K for each bend: ');

D = inputD / 1000; % Convert mm to m

%% Manufacturer recommended flow range

Qrec = linspace(Qmin, Qmax, 500); % L/min

Q = Qrec ./ (1000 * 60); % Convert L/min to m^3/s

%% Find pipe cross-sectional area and velocity

A = pi * D^2 / 4; % m^2

V = Q ./ A; % m/s, average gas velocity

%% Find Reynolds number

R = rho .* V .* D ./ mu;

%% Find friction factor and flow regime

f = zeros(size(R)); % Store friction factor values
flow = strings(size(R));  % Store flow regime

for i = 1:length(R)

    if R(i) < 2300

        flow(i) = "Laminar";
        % Darcy friction factor for laminar flow
        f(i) = 64 / R(i);

    elseif R(i) <= 4000

        flow(i) = "Transitional";
        % Preliminary smooth-pipe approximation
        f(i) = 0.3164 / (R(i)^0.25);

    else

        flow(i) = "Turbulent";
        % Blasius approximation for smooth pipe
        f(i) = 0.3164 / (R(i)^0.25);

    end

end

%% Estimate pressure loss

% Pressure loss from straight PVC pipe
dPpipe = f .* (L / D) .* (rho .* V.^2 ./ 2);

% Pressure loss from bends
Ktotal = numBends * Kbend;

dPbends = Ktotal .* (rho .* V.^2 ./ 2);

% Total pressure loss
dP = dPpipe + dPbends;

%% Display results at 50, 65, and 80 L/min

flows = [50 65 80];

fprintf('\n');
fprintf('PVC FLOW RESULTS\n');
fprintf('\n');

for i = 1:length(flows)

    velocity = interp1(Qrec, V, flows(i));
    reynolds = interp1(Qrec, R, flows(i));
    friction = interp1(Qrec, f, flows(i));

    pipePressure = interp1(Qrec, dPpipe, flows(i));
    bendPressure = interp1(Qrec, dPbends, flows(i));
    totalPressure = interp1(Qrec, dP, flows(i));

    if reynolds < 2300

        regime = "Laminar";

    elseif reynolds <= 4000

        regime = "Transitional";

    else

        regime = "Turbulent";

    end

    fprintf('\nAt %.0f L/min:\n', flows(i));
    fprintf('Velocity = %.3f m/s\n', velocity);
    fprintf('Reynolds number = %.0f\n', reynolds);
    fprintf('Flow regime = %s\n', regime);
    fprintf('Darcy friction factor = %.4f\n', friction);
    fprintf('Straight PVC pressure loss = %.2f Pa\n', pipePressure);
    fprintf('Bend pressure loss = %.2f Pa\n', bendPressure);
    fprintf('TOTAL pressure loss = %.2f Pa\n', totalPressure);

end

%% Fan/blower design requirements

P50 = interp1(Qrec, dP, 50);
P65 = interp1(Qrec, dP, 65);
P80 = interp1(Qrec, dP, 80);

fprintf('\n');
fprintf('FAN / BLOWER DESIGN REQUIREMENTS\n');
fprintf('\n');

fprintf('PVC inner diameter: %.1f mm\n', inputD);
fprintf('PVC pipe length: %.1f m\n', L);
fprintf('Number of bends: %.0f\n', numBends);
fprintf('K value per bend: %.2f\n\n', Kbend);

fprintf('At 50 L/min, fan must provide at least %.2f Pa\n', P50);
fprintf('At 65 L/min, fan must provide at least %.2f Pa\n', P65);
fprintf('At 80 L/min, fan must provide at least %.2f Pa\n', P80);

fprintf('\nRequired pressure range over WATO flow range: %.2f to %.2f Pa\n', ...
    P50, P80);

%% Plot system resistance

figure;

plot(Qrec, dPpipe);
hold on;

plot(Qrec, dP);

grid on;

xlabel('Flow Rate (L/min)');
ylabel('Pressure Loss (Pa)');
title('PVC Scavenging Pathway Pressure Loss');

legend('Straight PVC Only', ...
       'PVC + Bends', ...
       'Location', 'northwest');

hold off;