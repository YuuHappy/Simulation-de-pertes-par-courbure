clear all
close all
clc

tic;



%% -------------------------------
%% 🔁 Sweep parameters
%% -------------------------------
lambda_values = linspace(1060e-9, 1060e-9, 1);
Lval_values   = linspace(75e-6, 75e-6, 1);
RoC_values    = linspace(27.5e-3, 27.5e-3, 1);
NA_values     = linspace(0.0779, 0.0779, 1); 
Lz_values     = linspace(1e-2, 1e-2, 1);  

%% -------------------------------
%% 📦 Preallocate results
%% -------------------------------
total_runs = length(lambda_values)*length(Lval_values)* ...
             length(RoC_values)*length(NA_values)*length(Lz_values);

results = zeros(total_runs,6);
% [lambda, Lval, RoC, NA, Lz, power]

run_id = 1;

%% -------------------------------
%% 🔁 FULL PARAMETER SWEEP
%% -------------------------------
for il = 1:length(lambda_values)
for iv = 1:length(Lval_values)
for ir = 1:length(RoC_values)
for ina = 1:length(NA_values)
for iz = 1:length(Lz_values)

    lambda = lambda_values(il);
    Lval   = Lval_values(iv);
    RoC    = RoC_values(ir);
    NA     = NA_values(ina);
    Lz     = Lz_values(iz);

    fprintf('Run %d/%d | λ=%.0f nm | L=%.1f um | RoC=%.1f mm | NA=%.3f | Lz=%.1f mm\n',...
        run_id, total_runs,...
        lambda*1e9, Lval*1e6, RoC*1e3, NA, Lz*1e3);

    %% ---------------------------
    %% Initialize model
    %% ---------------------------
    P = BPMmatlab.model;

    P.name = 'FD_BPM_sweep';
    P.useAllCPUs = true;   % ✅ you CAN use all CPUs now
    P.useGPU = false;
    P.calcModeOverlaps = false;

    updatestepsize = 0.1e-3;

    %% Assign parameters
    P.lambda = lambda;
    P.bendingRoC = RoC;
    P.Lz = Lz;

    %% Grid
    P.Lx_main = Lval/sqrt(2);
    P.Ly_main = Lval/sqrt(2);
    P.Nx_main = round(200/50e-6*Lval);
    P.Ny_main = round(200/50e-6*Lval);
    P.padfactor = 1.5;
    P.dz_target = 0.5e-6;
    P.alpha = 3e14;

    %% Fiber parameters
    a_core = 0.5*(12E-6);
    n_clad = 1.4515;
    a_clad = 0.5*(80E-6);
    n_fa = 1.37;

    % ✅ NA sweep
    n_core = sqrt(n_clad^2 + NA^2);

    P.n_0 = n_core;
    P.n_background = n_clad;

    P = initializeRIfromFunction(P,@calcRI,{n_core,a_core,n_clad,a_clad});

    %% Modes
    P.bendDirection = 0;
    P = findModes(P,8);

    %% Inject LP01
    P.E = P.modes(1);

    %% Air cladding
    P.n_background = n_fa;
    P = initializeRIfromFunction(P,@calcRI,{n_core,a_core,n_clad,a_clad});

    %% Propagation
    P.updates = ceil(P.Lz/updatestepsize);

    %% Run BPM
    P = FD_BPM(P);

    %% Store results
    results(run_id,:) = [
        lambda*1e9,...
        Lval*1e6,...
        RoC*1e3,...
        NA,...
        Lz*1e3,...
        P.powers(end)
    ];

    run_id = run_id + 1;

end
end
end
end
end

%% -------------------------------
%% 📋 Convert to table
%% -------------------------------
results_table = array2table(results,...
    'VariableNames', { ...
    'Wavelength_nm',...
    'Lval_um',...
    'RoC_mm',...
    'NA',...
    'Lz_mm',...
    'FinalPower'});

%% -------------------------------
%% 💾 Save file
%% -------------------------------
filename = 'sweep_full.csv';
counter = 1;

while isfile(filename)
    filename = sprintf('sweep_full_%d.csv', counter);
    counter = counter + 1;
end

writetable(results_table, filename);

fprintf('\n✅ Sweep complete → %s\n', filename);

elapsedTime = toc;
fprintf('⏱️ Total time: %.2f seconds\n', elapsedTime);


%% -------------------------------
%% RI FUNCTION
%% -------------------------------
function n = calcRI(X,Y,n_background,nParameters)
n = n_background*ones(size(X));
n(X.^2 + Y.^2 < nParameters{4}^2) = nParameters{3};
n(X.^2 + Y.^2 < nParameters{2}^2) = nParameters{1};
end