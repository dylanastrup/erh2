% Jacob Wolf
% Heat Transfer for the copper pipe in ERH2's system
% Started 12/1/2022

clear all
clc

% Air Properties for Interpolating
temps = [20,25,30,35,40,45,50,60,70,80,90,100,120,140,160,180,200,250,300,350];  % Temperatures for interpolation, °C
specHeat = [1007,1007,1007,1007,1007,1007,1007,1007,1007,1008,1008,1009,1011,1013,1016,1019,1023,1033,1044,1056];  % Specific Heats of air at temps, J/kg*K
thermCond = [2.514,2.551,2.588,2.625,2.662,2.699,2.735,2.808,2.881,2.953,3.024,3.095,3.235,3.374,3.511,3.646,3.779,4.104,4.418,4.721]*10^-2;    % Thermal conductivities of air, W/m*K
dynVisc = [1.825,1.849,1.872,1.895,1.918,1.941,1.963,2.008,2.052,2.096,2.139,2.181,2.264,2.345,2.42,2.504,2.577,2.76,2.934,3.101]*10^-5;   % Dynamic Viscosity, kg/m*s
kinVisc = [1.516,1.562,1.608,1.655,1.702,1.750,1.798,1.896,1.995,2.097,2.201,2.306,2.522,2.745,2.975,3.212,3.455,4.091,4.765,5.475]*10^-5;    % Kinematic Viscosity, m^2/s
pran = [.7309,.7296,.7282,.7268,.7255,.7241,.7228,.7202,.7177,.7154,.7132,.7111,.7073,.7041,.7014,.6992,.6974,.6946,.6935,.6937];     % Prandtl Number

tempsExpans = [20,25,30,40,50,60,80,100,125,150,175,200,225,300,400];   % Temperatures for thermal expansion interpolation, °C
expans = [3.43,3.38,3.32,3.21,3.12,3.02,2.85,2.70,2.51,2.33,2.22,2.10,2.01,1.76,1.52]*10^-3;   % Thermal Expansion Coefficients, 1/K

% Bulk Air Properties
Tbulk = 20;     % Temp of bulk air, 20°C
grav = 9.81;    % Gravitational Acceleration, m/s^2


% Values to record from iterations
Texits = [];   % Exit Temperature of section, °C
Qdots = [];     % Heat transfer from secton, W


% Values to change/control
cylLength = 6*2.54/100;    % Length of Cylinder, m
inletTemp = 300;        % Inlet temp of air into cylinder, °C
exitTemp = 150;         % Exit temp of air from cylinder, °C
secLength = .1*2.54/100;    % Length of each section, m
diaOut = 0.5*2.54/100;       % Outer Diamater of pipe, m
diaIn = 0.45*2.54/100;       % Inner Diameter of pipe, m
massFlow = 0.03333;          % Mass Flow Rate in pipe, kg/s


% Process for each iteration: calculate internal & external convection
% coefficients, thermal resistance network, heat transfer through section,
% exit temperature. 

Texit = 325; % Temperature at exit from cylinder section, °C



% Interpolation: find where to interpolate to, find value. 
% For internal flow, at Texit exit temperature of previous section
tempInd = find(Texit<=temps,1);    % Finds first index in temps that the exit temp is under
tempDist = (Texit - temps(tempInd - 1))/(temps(tempInd) - temps(tempInd -1));    % Nondimensional distance from temp lower boundary to temp
specHeatIt = specHeat(tempInd-1)*(1-tempDist) + specHeat(tempInd)*tempDist;     % Specific Heat at that temperature
thermCondIt = thermCond(tempInd-1)*(1-tempDist) + thermCond(tempInd)*tempDist;  % Thermal Conductivity at that temp
dynViscIt = dynVisc(tempInd-1)*(1-tempDist) + dynVisc(tempInd)*tempDist;        % Dynamic Viscosity at that temp
kinViscIt = kinVisc(tempInd-1)*(1-tempDist) + kinVisc(tempInd)*tempDist;        % Kinematic Viscosity at that temp
pranIt = pran(tempInd-1)*(1-tempDist) + pran(tempInd)*tempDist;                 % Prandtl number at that temp

% Interpolation for External Convection Values, at Tfilm film temperature
Tsurf = 100;%Tsurf = (Texit + Tbulk)/2;       % Surface temp at section (mean of internal & external air)
Tfilm = (Tsurf + Tbulk)/2;          % Film Temperature, °C
tempFilmInd = find(Tfilm<=temps,1);
tempFilmDist = (Tfilm - temps(tempFilmInd-1))/(temps(tempFilmInd) - temps(tempFilmInd-1));
thermCondFilm = thermCond(tempFilmInd-1)*(1-tempFilmDist) + thermCond(tempFilmInd)*tempFilmDist;
kinViscFilm = kinVisc(tempFilmInd-1)*(1-tempFilmDist) + kinVisc(tempFilmInd)*tempFilmDist;
pranFilm = pran(tempFilmInd-1)*(1-tempFilmDist) + pran(tempFilmInd)*tempFilmDist;


% Interpolation for Thermal Expansion Coefficient, at Tfilm film
% temperature
tempExpansInd = find(Tfilm<=tempsExpans,1);
tempExpansDist = (Tfilm - tempsExpans(tempExpansInd - 1))/(tempsExpans(tempExpansInd) - tempsExpans(tempExpansInd - 1));
expansFilm = expans(tempExpansInd - 1)*(1 - tempExpansDist) + expans(tempExpansInd)*tempExpansDist;


% External Convection Values
grash = grav*expansFilm*(Tsurf - Tbulk)*cylLength^3/kinViscFilm^2;   % Grashof Number
ray = grash*pranFilm;    % Rayleigh Number
nussFlat = 0.68+0.67*ray^0.25/(1+(0.492/pranFilm)^(9/16))^(4/9); % Nusselt Number of flat plat
nussCyl = nussFlat*(1+0.3*(32^.5*grash^(-0.25)*cylLength/diaOut)^0.909);  % Nusselt Number of the cylinder
heatCoefExt = nussCyl*thermCondFilm/cylLength;       % Convective Heat Transfer Coefficient, External, W/m^2*K


