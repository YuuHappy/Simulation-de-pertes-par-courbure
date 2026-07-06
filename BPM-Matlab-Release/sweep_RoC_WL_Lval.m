clear all
close all
clc

tic;

%% --- Sweep parameters ---
lambda_values = linspace(1060e-9, 1060e-9, 1);
Lval_values   = linspace(55e-6, 55e-6, 1);    
RoC_values    = linspace(27.5e-3, 27.5e-3, 1);      % meters

%% --- Preallocate ---
total_runs = length(lambda_values)*length(Lval_values)*length(RoC_values);

results = zeros(total_runs,4); % [lambda, Lval, RoC, power]

run_id = 1;

for il = 1:length(lambda_values)
for iv = 1:length(Lval_values)
for ir = 1:length(RoC_values)

    fprintf('Run %d/%d | lambda=%.0f nm | L=%.1f um | RoC=%.1f mm\n',...
        run_id, total_runs,...
        lambda_values(il)*1e9,...
        Lval_values(iv)*1e6,...
        RoC_values(ir)*1e3);

    %% --- Initialize model ---
    P = BPMmatlab.model;

    P.name = 'FD_BPM_sweep';
    P.useAllCPUs = true;
    P.useGPU = false;
    P.calcModeOverlaps = false;

    updatestepsize = 0.1e-3;

    %% ✅ Assign swept parameters
    P.lambda = lambda_values(il);
    Lval     = Lval_values(iv);
    P.bendingRoC = RoC_values(ir);

    %% Grid
    P.Lx_main = Lval/sqrt(2);
    P.Ly_main = Lval/sqrt(2);
    P.Nx_main = round(200/50e-6*Lval);
    P.Ny_main = round(200/50e-6*Lval);
    P.padfactor = 1.5;
    P.dz_target = 0.5e-6;
    P.alpha = 3e14;

    %% Fiber parameters
    NA = 0.0779;
    a_core = 0.5*(12E-6);
    n_clad = 1.4515;
    %P.n_clad = n_clad;
    a_clad = 0.5*(80E-6);
    n_fa = 1.37;

    n_core = sqrt(n_clad^2 + NA^2);

    P = initializeRIfromFunction(P,@calcRI,{n_core,a_core,n_clad,a_clad});

    %% Bend + modes
    P.bendDirection = 0;

    P = findModes(P,8);

    %% Inject LP01
    P.E = P.modes(1);
    P.n_eff = P.modes(1).neff;

    %% Air background
    P.n_background = n_fa;
    P = initializeRIfromFunction(P,@calcRI,{n_core,a_core,n_clad,a_clad});

    %% Propagation
    P.Lz = 1e-2;
    P.updates = ceil(P.Lz/updatestepsize);

    %% Run
    P = FD_BPM(P);

    %% ✅ Store results
    results(run_id,:) = [
        P.lambda*1e9,...
        Lval*1e6,...
        P.bendingRoC*1e3,...
        P.powers(end)
    ];

    run_id = run_id + 1;

end
end
end

%% --- Convert to table ---
results_table = array2table(results,...
    'VariableNames', {'Wavelength_nm','Lval_um','RoC_mm','FinalPower'});

%% --- Save ---
filename = 'sweep_3D.csv';
counter = 1;

while isfile(filename)
    filename = sprintf('sweep_3D_%d.csv', counter);
    counter = counter + 1;
end

writetable(results_table, filename);

fprintf('\nSweep complete. Saved to %s\n', filename);

elapsedTime = toc;
fprintf('Total execution time: %.2f seconds\n', elapsedTime);


%% --- USER RI FUNCTION ---
function n = calcRI(X,Y,n_background,nParameters)
n = n_background*ones(size(X));
n(X.^2 + Y.^2 < nParameters{4}^2) = nParameters{3};
n(X.^2 + Y.^2 < nParameters{2}^2) = nParameters{1};
end