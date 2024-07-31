clc; clear; close all; tic; addpath("Common");
% =================================================
%
%         Orbital Perturbation Comparison
%
%                   By Ian Byres
%
%  ================================================
%% Notes

% Numerically propagates orbit from initial velocity and position
% conditions. Computes unperturbed keplerian orbit as well as perturbed to
% compare results.

% Perturbations included:
%   - 3rd body gravitation (moon & sun)
%   - Earth Oblateness (using J2 effect, further terms ignored)
%   - Solar Radiation Pressure (cylindrical shadow model)
%   - Atmospheric drag (using USSA76 model)
%   - Relativistic correction to acceleration (types of frame dragging
%     ignored) {Source 3}

% All units converted first and computations done in [km,kg,s]

% By default includes ISS data for Jan 01, 2024 UTC inputted {Source 2}

% References:
%   1. Curtis, Howard D. Orbital Mechanics for Engineering Students. 
%      4th ed., Butterworth-Heinemann, 2021. 
%   2. https://spotthestation.nasa.gov/trajectory_data.cfm
%   3. McCarthy, Dennis D. “IERS TECHNICAL NOTE 21.” International Earth 
%      Rotation Service, July 1996, Accessed 8 Jan.2024 
%      https://ilrs.gsfc.nasa.gov/docs/1996/iers_1996_conventions.pdf

%% ===================== EDIT HERE ====================================

tStart = [2024, 1, 2, 12, 0, 0]; % [Year, month, day, hour, minute, second]
tInt = [0, 0, .5, 0]; % Integrating time interval in [days,hours,minutes,seconds]
tSpan = [1, 0, 0, 0]; % Integrating time span in [days,hours,minutes,seconds]

startPos = [-3389.35 ; -5868.95; 523.149]; % x;y;z [km]
startVel = [3.8816 ; -2.7829 ; -5.9884]; % x;y;z [km/s]

ScMass = 451378; % [kg]
ScCD = 1.8; % Drag coefficient of spacecraft
ScDragArea = 2040.50; % Area inducing drag [m^2]
solPArea = 1500; % Solar Pressure Area [m^2]
solPCoeff = 1.8; % Radiation pressure coefficient




%% ===================== DO NOT EDIT =================================
% Adjusting inputs
EarthRadAvg = 6378.1366; % [km]
muEarth = (3.986004418e14)/(1e9); % [km^3/s^2]
angSpeedEarth = 7.2921159e-5; % Rotation speed of Earth [rad/s]
J2Earth = 0.00108263; % J2 constant of Earth

muMoon = (4.9048695e12)/(1e9); % [km^3/s^2]
muSun = (1.32712440018e20)/(1e9); % [km^3/s^2]

cLight = 299792.458; % Speed of light in [km/s]
RadiationIntensity = 1367; % Radiation intensity [W/m^2] (m in W cancels with denominator)

solPArea = solPArea/(1e6);
ScDragArea = ScDragArea/(1e6); % Converting to [km^2]

tInt = (tInt(1)*86400)+(tInt(2)*3600)+(tInt(3)*60)+(tInt(4));
tSpan = (tSpan(1)*86400)+(tSpan(2)*3600)+(tSpan(3)*60)+(tSpan(4));

% Pulling Ephemeris Info & ODE Inputs
odeTime = (0:tInt:tSpan)'; % finding integrating time span for ode
ephTime = tStart(ones(length(odeTime),1), :);
for i = 1:length(ephTime)
    ephTime(i,6) = ephTime(i,6) + ((i-1)*tInt);
end
ephTime = juliandate(ephTime); % Creating array of individual ephemeris times 

sunPos = planetEphemeris(ephTime,'Earth','Sun'); % Position of sun relative to Earth at all time indices
moonPos = planetEphemeris(ephTime,'Earth','Moon'); % Position of moon relative to Earth at all time indices

clear i ephTime;

initCond = [startPos;startVel]; % Setting initial conditions
inputArray = [muEarth; ScMass; angSpeedEarth; ScCD; ScDragArea; ...
              EarthRadAvg; J2Earth; cLight; muMoon; muSun; ...
              RadiationIntensity; solPCoeff; solPArea];

% Integrating ODE
[~,QUPSC] = ode89(@(t,q) unPurtOrbit(t,q,muEarth), odeTime, initCond); % Unpurturbed SC
QUPSC(:,4:6) = []; % Clears Velocity array to save workspace capacity

[~,QSC] = ode89(@(t,q) purtOrbit(t,q,inputArray,odeTime,sunPos,moonPos), odeTime, initCond); % Purturbed SC
QSC(:,4:6) = []; % Clears Velocity array to save workspace capacity
                 % tUPSC and tSC removed for space b/c they match odeTime

orbitPositionError = vecnorm((QSC-QUPSC)');

% Plotting Sphere
f1 = figure(1);
f1.Name = 'Orbit Plots';
hold on;

[XEIMG, map] = rgb2ind(imread('Earth.jpg'),128);
[xEIMG,yEIMG,zEIMG] = sphere(50);
xEIMG = EarthRadAvg*xEIMG;
yEIMG = EarthRadAvg*yEIMG;
zEIMG = EarthRadAvg*zEIMG;
props.FaceColor = 'texture';
props.EdgeColor = 'none';
props.Cdata = flipud(XEIMG);
surface(xEIMG,yEIMG,zEIMG,props);
colormap(map);

% Plotting Orbits
plot3(QUPSC(1,1),QUPSC(1,2),QUPSC(1,3),'linewidth',4,'Marker','o','Color','w')
plot3(QSC(end,1),QSC(end,2),QSC(end,3),'linewidth',4,'Marker','o','Color','b')
plot3(QUPSC(:,1),QUPSC(:,2),QUPSC(:,3),'linewidth',2.25,'Color','#14c400');
plot3(QSC(:,1),QSC(:,2),QSC(:,3),'linewidth',0.5,'Color','#b51818');

% Formatting Plots
title('Trajectory Propagations');
xlabel('x [km]');
ylabel('y [km]');
zlabel('z [km]');
legend('','Starting Position','End Position','Unpurturbed Trajectory'...
    ,'Purturbed Trajectory');
set(gca,'Color','#0a0a0d');
axis equal;
grid off;
view([18 29]);
f1.Color = "#0a0a0d";
f1.WindowState = 'maximized';

% Plotting Perturbation Measurement
f2 = figure(2);
f2.Name = 'Position Error Buildup';
hold on;
plot((odeTime/3600),orbitPositionError,"LineWidth",1.5); % Plots in hours
title('Keplerian Orbit Position Error Magnitude vs. Time','Interpreter','latex');
xlabel('Time [hrs]');
ylabel('Difference in Position [km]');
legend('Position Difference in [km]');
grid on;
f2.WindowState = 'maximized';

% Displaying Final Position
fprintf('Ending position of perturbed S/C in [km]: \n')
disp(QSC(end,:));

toc;
%% ================= FUNCTIONS ====================
% Purturbed Orbit Function
function qdot = purtOrbit(t,q,inputArray,odeTime,sunPosArray,moonPosArray)
    % Pulling misc variables from inputArray
    muEarth = inputArray(1);
    ScMass = inputArray(2);
    angSpeedEarth = inputArray(3);
    ScCD = inputArray(4);
    ScDragArea = inputArray(5);      % Should hard code some values in here
    EarthRadAvg = inputArray(6);
    J2Earth = inputArray(7);
    cLight = inputArray(8);
    muMoon = inputArray(9);
    muSun = inputArray(10);
    RadiationIntensity = inputArray(11);
    solPCoeff = inputArray(12);
    solPArea = inputArray(13);

    % Pulling R and V vector from input q
    R = [q(1);q(2);q(3)]; % vector position
    x = R(1); y = R(2); z = R(3); % x,y,z components of position
    r = vecnorm(R); % magnitude of position
    V = [q(4);q(5);q(6)]; % vector of velocity
    v = vecnorm(V); % magnitude of velocity

    % Finding current index # for ephemeris data
    [~,ephIndex]=min(abs(odeTime-t));

    SunPos = sunPosArray(ephIndex,:)'; % Sun wrt Earth
    sunPos = vecnorm(SunPos);
    RSun = SunPos - R; % Sun wrt SC
    rSun = vecnorm(RSun);
    MoonPos = moonPosArray(ephIndex,:)'; % Moon wrt Earth
    moonPos = vecnorm(MoonPos);
    RMoon = MoonPos - R; % Moon wrt SC
    rMoon = vecnorm(RMoon);

    % ----------------------------------

    % All accelerations computed in [km/s^2]

    % Gravitational acceleration from Earth
    a = -(muEarth/(r^3))*R; % Keplerian acceleration from Earth

    % Gravitational acceleration from Moon
    MoonA = muMoon*((RMoon/(rMoon^3)) - (MoonPos/(moonPos^3)));
    a = a + MoonA; % [Comment if undesired] added 3rd body lunar acceleration

    % Gravitational acceleration from Sun
    SunA = muSun*((RSun/(rSun^3)) - (SunPos/(sunPos^3)));
    a = a + SunA; % [Comment if undesired] added 3rd body solar acceleration

    % Atmospheric Drag
    VrelAtm = V - cross(angSpeedEarth*[0;0;1],R); % Velocity relative to atmosphere's air (due to Earths rotation)
    rho = (atmosphere((vecnorm(R)-EarthRadAvg))) * 1000^3; % Air density at altitude (USSA76) in [kg/km^3]
    Drag = (-0.5 *rho*vecnorm(VrelAtm)*ScCD*ScDragArea)*VrelAtm;
    DragA = Drag/ScMass;
    a = a + DragA; % [Comment if undesired] added atmospheric acceleration

    % Earth Oblateness (J2)
    J2A = [ (x/r)*(5*((z^2)/(r^2))-1)  ; ...
            (y/r)*(5*((z^2)/(r^2))-1)  ; ...
            (z/r)*(5*((z^2)/(r^2))-3) ];
    J2A = (3*J2Earth*muEarth*(EarthRadAvg^2))/(2*(r^4)) * J2A;
    a = a + J2A; % [Comment if undesired] added oblateness acceleration

    % Solar Pressure
    solPTheta = acos(( (dot(RSun,R))/(rSun*r) ));
    solPTheta1 = acos(EarthRadAvg/r);
    solPTheta2 = acos(EarthRadAvg/rSun);
    if (solPTheta1 + solPTheta2) <= solPTheta
        shadow = 0;
    else
        shadow = 1;
    end
    SolarPressure = ((-shadow*RadiationIntensity*solPCoeff*solPArea)/(cLight*ScMass));
    SolPA = SolarPressure*(RSun/rSun);
    a = a + SolPA; % [Comment if undesired] added solar pressure acceleration

    % Relativistic Acceleration
    RelA = (muEarth/((cLight^2)*(r^3)))*(((4*muEarth/r)-v^2)*R+(4*dot(R,V)*V));
    a = a + RelA; % [Comment if undesired] added relativistic acceleration

    % ----------------------------------

    % Sending new velocity & acceleration data into loop
    qdot = [q(4);q(5);q(6);a(1);a(2);a(3)];
end




% Unpurturbed Orbit Function
function qdot = unPurtOrbit(~,q,muEarth)
    % Pulling R and V vector from input q
    R = [q(1);q(2);q(3)];
    % V = [q(4);q(5);q(6)];

% ----------------------------------
    % Gravitational force from Earth
    rEarth = vecnorm(R);
    a = -(muEarth/(rEarth^3))*R;
% ----------------------------------

    % Sending new velocity & acceleration data into loop
    qdot = [q(4);q(5);q(6);a(1);a(2);a(3)];
end




