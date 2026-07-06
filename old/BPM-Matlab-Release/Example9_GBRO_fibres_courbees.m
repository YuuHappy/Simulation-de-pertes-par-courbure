clear all
close all
%clc

P = BPMmatlab.model;

% This example shows the propagation of an LP mode in a LMA fiber.
% The LMA fiber is divided into 7 segments, where some segments are
% straight, some are bent in the x direction and some are bent in
% the y direction. Plotting of the mode overlaps has been enabled by
% setting P.calcModeOverlaps = true.

% We use the getLabeledModeIdx() function to find out which index in the
% P.modes array correponds to a particular mode label (LP01, LP11e, etc.)

% ################################## ATTENTION ##################################
% Ce code ne peut considérer un x_trans comme à la Fig. 3 de Schermer2007.
% Le calcul de pertes par courbure est alors faussé pour les petits rayons
% de courbure sur un domaine de calcul très large. Il faut s'assurer de 
% réduire le domaine de calcul (significativement plus petit que le
% cladding) afin d'avoir des résultats valables
% ###############################################################################

%% General and solver-related settings
P.name = mfilename;
P.useAllCPUs = true; % If false, BPM-Matlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
P.useGPU = false; % !! (Default: false) Use CUDA acceleration for NVIDIA GPUs

%% Visualization parameters
P.calcModeOverlaps = true;  % Set it to true to calculate mode overlap integrals of propagating field with respect to different modes in the P.modes struct array
%updatestepsize = 1e-5;
updatestepsize = 0.1e-3;

%% Resolution-related parameters (check for convergence)
P.Lx_main = 30e-6;        % [m] x side length of main area
P.Ly_main = 30e-6;        % [m] y side length of main area
P.Nx_main = 150;          % x resolution of main area
P.Ny_main = 150;          % y resolution of main area
P.padfactor = 1.5;  % How much absorbing padding to add on the sides of the main area (1 means no padding, 2 means the absorbing padding on both sides is of thickness Lx_main/2)
P.dz_target = 0.5e-6; % [m] z step size to aim for
P.alpha = 3e14;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area

%% Problem definition
P.lambda = 1064e-9; % [m] Wavelength
P.n_background = 1.45; % [] (may be complex) Background refractive index, (in this case, the cladding)
a_core = 5E-6;
NA = 0.09;
P.n_0 = sqrt(P.n_background^2 + NA^2); % The reference refractive index, see README.md

P = initializeRIfromFunction(P,@calcRI,{P.n_0,a_core});

P = findModes(P,10); % Find up to 10 modes

%Mode injecté
modeIdx = getLabeledModeIndex(P,'LP01'); % LP01 mode
%modeIdx = getLabeledModeIndex(P,'LP11e'); % LP11e mode

P.E = P.modes(modeIdx);

%% Segment 1
P.Lz = 5e-3; % [m] z propagation distances for this segment
P.updates = ceil(P.Lz/updatestepsize);
P.bendDirection = 0; % [degrees] direction of the bending, in a polar coordinate system with 0° to the right (towards positive x) and increasing angles in counterclockwise direction
P.bendingRoC = Inf; % [m] radius of curvature of the bend

P = FD_BPM(P);

for jj=1:2
    %% Segment 2
    P.bendingRoC = 15e-3;
    P.Lz = P.bendingRoC*(pi/6); % [m] z propagation distances for this segment
    P.updates = ceil(P.Lz/updatestepsize);
    P.bendDirection = 0;

    P = FD_BPM(P);

    %% Segment 3
    P.bendingRoC = 15e-3;
    P.Lz = P.bendingRoC*(pi/6); % [m] z propagation distances for this segment
    P.updates = ceil(P.Lz/updatestepsize);
    P.bendDirection = 90;

    P = FD_BPM(P);

    %% Segment 4
    P.Lz = 15e-3*(pi/6); % [m] z propagation distances for this segment
    P.updates = ceil(P.Lz/updatestepsize);
    P.bendDirection = 0;
    P.bendingRoC = Inf;

    P = FD_BPM(P);
end

%% Segment 5
P.Lz = 20e-3; % [m] z propagation distances for this segment
P.updates = ceil(P.Lz/updatestepsize);
P.bendDirection = 0;
P.bendingRoC = Inf;

P = FD_BPM(P);

%% USER DEFINED RI FUNCTIONS
function n = calcRI(X,Y,n_background,nParameters)
disp(nParameters)
% n may be complex
n = n_background*ones(size(X)); % Start by setting all pixels to n_background
n(X.^2 + Y.^2 < nParameters{2}^2) = nParameters{1};
end