clear all
close all
clc

sweepTimer = tic;


%% Ask user for output folder name

folderName = input('Enter output folder name: ','s');

if isempty(folderName)
    folderName = datestr(now,'yyyy-mm-dd_HH-MM-SS');
end

% Prevent overwriting existing folder
baseFolder = folderName;
counter = 1;

while exist(folderName,'dir')
    folderName = sprintf('%s_%d',baseFolder,counter);
    counter = counter + 1;
end

mkdir(folderName);

fprintf('Results will be saved in:\n%s\n\n',folderName);


%% -------------------------------
%% 🔁 Sweep parameters
%% -------------------------------
lambda_values = linspace(850e-9, 1100e-9, 6);
Lval_values   = linspace(55e-6, 55e-6, 1);
RoC_values    = linspace(20.0e-3, 20.0e-3, 1);
NA_values     = linspace(0.0779, 0.0779, 1); 
Lz_values     = linspace(1e-2, 1e-1, 2);  

%% -------------------------------
%% 📦 Preallocate results
%% -------------------------------
total_runs = length(lambda_values)*length(Lval_values)* ...
             length(RoC_values)*length(NA_values)*length(Lz_values);

results = zeros(total_runs,9);
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
    
    if run_id>1
        while exist(folderName,'dir')
        folderName = sprintf('%s_%d',baseFolder,counter);
        counter = counter + 1;
        end
        mkdir(folderName);
        fprintf('Results will be saved in:\n%s\n\n',folderName);
    end
    fprintf('Run %d/%d | λ=%.0f nm | L=%.1f um | RoC=%.1f mm | NA=%.3f | Lz=%.1f mm\n',...
        run_id, total_runs,...
        lambda*1e9, Lval*1e6, RoC*1e3, NA, Lz*1e3);

    %% ---------------------------
    %% Initialize model
    %% ---------------------------
    P = BPMmatlab.model;
    P.saveVideo = true;

    P.name = fullfile(folderName, folderName);
    P.useAllCPUs = true;   % ✅ you CAN use all CPUs now
    P.useGPU = false;
    P.calcModeOverlaps = false;

    updatestepsize = 0.1e-3;

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
    P.n_clad = n_clad;

    % ✅ NA sweep
    n_core = sqrt(n_clad^2 + NA^2);
    P.n_0 = n_core;
    P.n_background = n_fa;

    P = initializeRIfromFunction(P,@calcRI,{n_core,a_core,n_clad,a_clad});

    %% Paramètres initiaux

    P.bendingRoC = Inf;
    P.lambda = lambda;
    P.Lz = Lz;

    %% Modes
    P.bendDirection = 0;
    P = findModes(P,8);
    P.modes(1)
    P.n_eff = P.modes(1).neff;
    fprintf('N_eff is %.5f\n', P.n_eff);

    %% Inject LP01
    P.E = P.modes(1);

    %% Assign parameters
    P.bendingRoC = RoC;

    %% Air cladding
    %P.n_background = n_fa;
    P = initializeRIfromFunction(P,@calcRI,{n_core,a_core,n_clad,a_clad});


    %% =====================================================
    %% Adiabatic transition from straight fiber to target RoC
    %% =====================================================
    
    RoC_final = P.bendingRoC;
    
    % Initial straight launch section
    P.bendingRoC = Inf;
    P.bendDirection = 0;
    
    Nseg = 20; %number of segments
    Lz_straight = 1e-3; %initial straight launch
    Lz_ramp_total = 4e-3;

    P.Lz = Lz_straight;
    P.updates = ceil(P.Lz/updatestepsize);
    
    P.figTitle = 'Straight';
    
    P = FD_BPM(P);
    
    % Radius ramp (large -> final)
    RoC_ramp = 1 ./ linspace(0,1/RoC_final,Nseg+1);
    RoC_ramp(1) = Inf;
    RoC_ramp = RoC_ramp(2:end);

    % Total adiabatic length
    Lz_seg = Lz_ramp_total/Nseg;
    
    for k = 1:Nseg
    
        P.bendingRoC = RoC_ramp(k);
        P.Lz = Lz_seg;
        P.updates = ceil(P.Lz/updatestepsize);
        
        P.figTitle = sprintf('Ramp %d/%d | RoC=%.2f mm', ...
            k,Nseg,P.bendingRoC*1e3);

        P = FD_BPM(P);
    
    end

    %store data of power without insertion loss
    Powerafterinjection = P.powers(end);
    
    %% Main bent section
    P.bendingRoC = RoC_final;
    P.Lz = Lz-Lz_straight-Lz_ramp_total;
    P.updates = ceil(P.Lz/updatestepsize);
    
    
    P.figTitle = sprintf('Main bend | RoC=%.2f mm', P.bendingRoC*1e3);

    P = FD_BPM(P);
    

    if P.saveVideo && ~isempty(P.videoHandle)
        close(P.videoHandle);
    end


    %% Store results
    %Modify data for storage
    ZoneCalc = Lval*1e6/sqrt(2);
    Loss_dB = -10*log10(P.powers(end)/Powerafterinjection);
    dbsurm = Loss_dB/(Lz-Lz_straight-Lz_ramp_total);
    
    %store
    results(run_id,:) = [
        NA,...
        lambda*1e9,...
        Lz,...
        Lval*1e6,...
        ZoneCalc,...
        RoC*1e3,...
        P.powers(end),...
        dbsurm,...
        Powerafterinjection
    ];

    run_id = run_id + 1;

    hPower = figure('Visible','off');
    
    plot(P.z,P.powers,'LineWidth',2);
    hold on;
    
    xline(Lz_straight,'--k','Straight End');
    xline(Lz_straight + Lz_ramp_total,'--r','Ramp End');
    
    grid on;
    
    xlabel('Propagation distance [m]');
    ylabel('Relative power remaining');
    title('Relative Power Evolution');
    
    exportgraphics( ...
        hPower,...
        fullfile(folderName,[folderName '_Power.png']),...
        'Resolution',300);
    
    close(hPower);

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
    'NA',...
    'Wavelength_nm',...
    'Distance (m)',...
    'Lval_um',...
    'Zone calcul', ...
    'RoC_mm',...
    'FinalPower',...
    'db/m',...
    'PowerAfterInsertion'});

%% -------------------------------
%% 💾 Save file
%% -------------------------------

csvfile = fullfile(folderName,[folderName '_data.csv']);
writetable(results_table,csvfile)


fprintf('\n✅ Sweep complete → %s\n', folderName);

elapsedTime = toc(sweepTimer);
fprintf('⏱️ Total time: %.2f seconds\n', elapsedTime);


%% -------------------------------
%% RI FUNCTION
%% -------------------------------
function n = calcRI(X,Y,n_background,nParameters)
n = n_background*ones(size(X));
n(X.^2 + Y.^2 < nParameters{4}^2) = nParameters{3};
n(X.^2 + Y.^2 < nParameters{2}^2) = nParameters{1};
end