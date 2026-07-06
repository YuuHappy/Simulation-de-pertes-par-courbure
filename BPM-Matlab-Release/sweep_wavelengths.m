clear all
close all
clc

tic;

% --- Define wavelength sweep (meters) ---
lambda_values = linspace(900e-9, 1250e-9, 1);  % adjust range here
num_points = length(lambda_values);

% --- Storage ---
final_power = zeros(num_points,1);

for i = 1:num_points
    
    fprintf('Run %d/%d | lambda = %.1f nm\n', ...
        i, num_points, lambda_values(i)*1e9);

    %% --- Reinitialize model each iteration ---
    P = BPMmatlab.model;

    %% General settings
    P.name = 'FD_BPM_sweep';
    P.useAllCPUs = true;
    P.useGPU = false;

    %% Visualization OFF for speed
    P.calcModeOverlaps = false;   % disable extra overhead

    updatestepsize = 0.1e-3;
    Lval = 55e-6;

    %% Grid settings
    P.Lx_main = Lval/sqrt(2);
    P.Ly_main = Lval/sqrt(2);
    P.Nx_main = round(200/50e-6*Lval);
    P.Ny_main = round(200/50e-6*Lval);
    P.padfactor = 1.5;
    P.dz_target = 0.5e-6;
    P.alpha = 3e14;

    %% ✅ Set wavelength (swept parameter)
    P.lambda = lambda_values(i);

    %% Fiber parameters
    NA = 0.0779;
    a_core = 0.5*(12E-6);
    n_clad = 1.4515;
    a_clad = 0.5*(80E-6);
    n_fa = 1.37;

    n_core = sqrt(n_clad^2 + NA^2);

    P.n_0 = n_core;
    P.n_background = n_clad;

    P = initializeRIfromFunction(P,@calcRI,{n_core,a_core,n_clad,a_clad});

    %% Bend + modes
    P.bendDirection = 0;
    P.bendingRoC = 27.5E-3;

    P = findModes(P,8);

    % Inject LP01
    modeIdx = 1;
    P.E = P.modes(modeIdx);

    % Change background if needed
    P.n_background = n_fa;
    P = initializeRIfromFunction(P,@calcRI,{n_core,a_core,n_clad,a_clad});

    %% Propagation
    P.Lz = 1e-2;
    P.updates = ceil(P.Lz/updatestepsize);

    %% Run BPM
    P = FD_BPM(P);

    %% ✅ Store final relative power
    final_power(i) = P.powers(end);
end

%% --- Save results ---
results_table = table(lambda_values(:)*1e9, final_power, ...
    'VariableNames', {'Wavelength_nm','FinalRelativePower'});



% --- Convert parameters to strings ---
RoC_str = sprintf('RoC_%.0fmm', P.bendingRoC * 1e3);   % meters → mm

% Replace '.' with '-' in NA
NA_raw = sprintf('%.3f', NA);
NA_clean = strrep(NA_raw, '.', '-');
NA_str = ['NA_' NA_clean];


% --- Base filename ---
base_filename = sprintf('sweep_%s_%s.csv', RoC_str, NA_str);

filename = base_filename;
counter = 1;

% --- Increment if file exists ---
while isfile(filename)
    filename = sprintf('sweep_%s_%s_%d.csv', RoC_str, NA_str, counter);
    counter = counter + 1;
end

% --- Save ---
writetable(results_table, filename);

fprintf('\nSweep complete. Results saved to %s\n', filename);


% --- Save ---
writetable(results_table, filename);

fprintf('\nSweep complete. Results saved to %s\n', filename);

elapsedTime = toc;
fprintf('Total execution time: %.2f seconds\n', elapsedTime);


%% --- USER RI FUNCTION ---
function n = calcRI(X,Y,n_background,nParameters)
n = n_background*ones(size(X));
n(X.^2 + Y.^2 < nParameters{4}^2) = nParameters{3};
n(X.^2 + Y.^2 < nParameters{2}^2) = nParameters{1};
end