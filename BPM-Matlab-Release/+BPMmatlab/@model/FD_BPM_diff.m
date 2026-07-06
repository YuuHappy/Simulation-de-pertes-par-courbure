function P = FD_BPM(P)

format long
format compact

% --- Basic checks ---
if isempty(P.n.n)
  error('Refractive index not initialized.');
end

if isempty(P.E.field)
  error('Electric field not initialized.');
end

% --- Constants ---
k_0 = 2*pi/P.lambda;
dx = P.dx; dy = P.dy;
Nx = P.Nx; Ny = P.Ny;

x = P.x;
y = P.y;

[X,Y] = ndgrid(x,y);

% --- Normalize input field ---
powerFraction = 1/(1 + (P.xSymmetry ~= 0))/(1 + (P.ySymmetry ~= 0));

if ~P.priorData
  P.E.field = P.E.field / sqrt(sum(abs(P.E.field(:)).^2)/powerFraction);
end

% --- Interpolate E field ---
[Nx_source,Ny_source] = size(P.E.field);
dx_source = P.E.Lx/Nx_source;
dy_source = P.E.Ly/Ny_source;

x_source = getGridArray(Nx_source,dx_source,P.E.ySymmetry);
y_source = getGridArray(Ny_source,dy_source,P.E.xSymmetry);

[x_source,y_source,E_source] = calcFullField(x_source,y_source,P.E.field);

E = interpn(x_source,y_source,E_source,x,y.','linear',0);
E = E * sqrt(sum(abs(E_source(:)).^2)/sum(abs(E(:)).^2/powerFraction));

E = complex(single(E));

% --- Interpolate refractive index ---
[Nx_source,Ny_source,Nz_source] = size(P.n.n);
dx_source = P.n.Lx/Nx_source;
dy_source = P.n.Ly/Ny_source;

x_source = getGridArray(Nx_source,dx_source,P.n.ySymmetry);
y_source = getGridArray(Ny_source,dy_source,P.n.xSymmetry);

[x_source,y_source,n_source] = calcFullRI(x_source,y_source,P.n.n);

n = interpn(x_source,y_source,n_source,x,y.','linear',P.n_background);

% --- Z stepping ---
Nz = P.Nz;
dz = P.dz;

zUpdateIdxs = round((1:P.updates)/P.updates * Nz);

% --- BPM constants ---
ax = dz/(4i*dx^2*k_0*P.n_0);
ay = dz/(4i*dy^2*k_0*P.n_0);
d  = -dz*k_0;

% --- Absorber ---
xEdge = P.Lx_main*(1 + (P.ySymmetry ~= 0))/2;
yEdge = P.Ly_main*(1 + (P.xSymmetry ~= 0))/2;

multiplier = single(exp(-dz * max(0,max(abs(Y)-yEdge,abs(X)-xEdge)).^2 * P.alpha));

% --- Initialize power storage ---
P.powers = zeros(1, length(zUpdateIdxs));

% --- MEX parameters ---
mexParameters = struct( ...
  'dx',single(dx),'dy',single(dy),'dz',single(dz), ...
  'taperPerStep',single((1-P.taperScaling)/Nz), ...
  'twistPerStep',single(P.twistRate*P.Lz/Nz), ...
  'multiplier',single(multiplier), ...
  'n_mat',complex(single(n)), ...
  'dz_n',single(0), ...
  'd',single(d), ...
  'n_0',single(P.n_0), ...
  'ax',single(ax),'ay',single(ay), ...
  'useAllCPUs',P.useAllCPUs, ...
  'RoC',single(P.bendingRoC), ...
  'rho_e',single(P.rho_e), ...
  'bendDirection',single(P.bendDirection), ...
  'inputPrecisePower',powerFraction, ...
  'xSymmetry',uint8(P.xSymmetry), ...
  'ySymmetry',uint8(P.ySymmetry));

mexParameters.iz_start = int32(0);
mexParameters.iz_end   = int32(zUpdateIdxs(1));

% --- Main propagation loop ---
for updidx = 1:length(zUpdateIdxs)

  if updidx > 1
    mexParameters.iz_start = int32(zUpdateIdxs(updidx-1));
    mexParameters.iz_end   = int32(zUpdateIdxs(updidx));
  end

  checkMexInputs(E,mexParameters);

  if P.useGPU
    [E,~,precisePower] = FDBPMpropagator_CUDA(E,mexParameters);
  else
    [E,~,precisePower] = FDBPMpropagator(E,mexParameters);
  end

  mexParameters.inputPrecisePower = precisePower;

  P.powers(updidx) = precisePower / powerFraction;
end

% --- Print final result ---
final_power = P.powers(end);
fprintf('Final relative power remaining: %.6f\n', final_power);

% --- Store back results ---
P.E.field = E;
P.n.n = n;
P.priorData = true;

end


% --- Input validation ---
function checkMexInputs(E,P)
assert(all(isfinite(E(:))));
assert(~isreal(E));
end
