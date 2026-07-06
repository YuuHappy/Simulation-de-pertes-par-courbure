clear all
close all
%clc
tic;   % start timer
P = BPMmatlab.model;

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
P.calcModeOverlaps = true;  % !! Set it to true to calculate mode overlap integrals of propagating field with respect to different modes in the P.modes struct array
%updatestepsize = 1e-5;
updatestepsize = 0.1e-3;

Lval = 55e-6;

%% Resolution-related parameters (check for convergence)
P.Lx_main = Lval/sqrt(2);        % [m] x side length of main area
P.Ly_main = Lval/sqrt(2);        % [m] y side length of main area
P.Nx_main = round(200/50e-6*Lval);          % x resolution of main area
P.Ny_main = round(200/50e-6*Lval);          % y resolution of main area
P.padfactor = 1.5;  % How much absorbing padding to add on the sides of the main area (1 means no padding, 2 means the absorbing padding on both sides is of thickness Lx_main/2)
P.dz_target = 0.5e-6; % [m] z step size to aim for
P.alpha = 3e14;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area

%% Problem definition
P.lambda = 1065e-9; % [m] Wavelength

%Fiber
%NA = 0.07;
NA = 0.0779;
a_core = 0.5*(12E-6);
n_clad = 1.4515;  %https://refractiveindex.info/?shelf=main&book=SiO2&page=Malitson
a_clad = 0.5*(80E-6);
n_fa = 1.37;

n_core = sqrt(n_clad^2 + NA^2); % [] reference refractive index

P.n_0 = n_core; %The reference refractive index, see README.md
P.n_background = n_clad; % [] (may be complex) Background refractive index, (in this case, the cladding)

P = initializeRIfromFunction(P,@calcRI,{n_core,a_core,n_clad,a_clad});

%% Segment 1 (segment courbé en partant)
P.bendDirection = 0; % [degrees] direction of the bending, in a polar coordinate system with 0° to the right (towards positive x) and increasing angles in counterclockwise direction
P.bendingRoC = 27.5E-3; % [m] radius of curvature of the bend

P = findModes(P,8); % Find first modes


%Mode injecté
modeIdx = 1;
%modeIdx = getLabeledModeIndex(P,'LP01'); % LP01 mode
%modeIdx = getLabeledModeIndex(P,'LP11e'); % LP11e mode

P.E = P.modes(modeIdx);

if(true) 
    P.n_background = n_fa; % [] (may be complex) Background refractive index, (in this case, the double-cladding)
    P = initializeRIfromFunction(P,@calcRI,{n_core,a_core,n_clad,a_clad});
end

P.Lz = 1e-1; % [m] z propagation distances for this segment
%P.Lz = 10e-3; % [m] z propagation distances for this segment
P.updates = ceil(P.Lz/updatestepsize);

P = FD_BPM(P);
    
if(false)
    
    warning('Calcul peu fiable, voir commentaire sur x_trans ci-haut'); 
    Perte_dB_m = -10*log10(P.modeOverlaps(modeIdx,2:end)./P.modeOverlaps(modeIdx,1))./((1:P.updates)*P.Lz/P.updates);
    figure; plot(Perte_dB_m); title(sprintf('PERTE PEU FIABLE : %.3f dB/m',Perte_dB_m(end)));
end
elapsedTime = toc;   % stop timer and return time

fprintf('Execution time: %.4f seconds\n', elapsedTime);
%% USER DEFINED RI FUNCTIONS
function n = calcRI(X,Y,n_background,nParameters)
disp(nParameters)
% n may be complex
n = n_background*ones(size(X)); % Start by setting all pixels to n_background
n(X.^2 + Y.^2 < nParameters{4}^2) = nParameters{3};
n(X.^2 + Y.^2 < nParameters{2}^2) = nParameters{1};
end


