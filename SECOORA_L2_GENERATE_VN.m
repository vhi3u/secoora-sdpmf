% L2 data processing for SECOORA gliders
% Victor Hieu Nguyen 
% UNC-Chapel Hill

clc
clear;

% load CTD data

load('salacia_2021_09_CTD_L1.mat')

outdir = sprintf('/Users/victornguyen/Documents/MATLAB/dataSECOORA/level2/%s/', config.glider_name);

fprintf('Starting L2 CTD processing... \nGlider name: %s \nGen: %s \nisPumped = %i\n', ...
    config.glider_name, config.gen, config.isPumped);

% load ECO data



%%
% create grid
dt = 1/24;
dn_vec = ceil(min(ptime_datenum)):dt:max(ptime_datenum);
dz = 1;
zmax = 200; % vertical grid size
z_vec = 0:dz:zmax; % z vector
z_vec(1) = -eps; % bin edge includes zero; use eps. 

% bin center
dn_mean = (dn_vec(2:end) + dn_vec(1:end-1)) / 2;
z_mean = (z_vec(2:end) + z_vec(1:end-1)) / 2; 

% nan CTD
notnans = find(~isnan(tempCorrected));

% bindata

% main 
[tempG, yb, ystd, nnall] = bindata2_old(tempCorrected(notnans), ...
    ptime_datenum(notnans), depth(notnans), dn_vec, z_vec);

[salinG, yb, ystd, nnall] = bindata2_old(salinCorrected(notnans), ...
    ptime_datenum(notnans), depth(notnans), dn_vec, z_vec);

[densG, yb, ystd, nnall] = bindata2_old(densCorrected(notnans), ...
    ptime_datenum(notnans), depth(notnans), dn_vec, z_vec);

% uncorrected
[tempU, yb, ystd, nnall] = bindata2_old(temp(notnans), ...
    ptime_datenum(notnans), depth(notnans), dn_vec, z_vec);

[salinU, yb, ystd, nnall] = bindata2_old(salin(notnans), ...
    ptime_datenum(notnans), depth(notnans), dn_vec, z_vec);

[densU, yb, ystd, nnall] = bindata2_old(dens(notnans), ...
    ptime_datenum(notnans), depth(notnans), dn_vec, z_vec);

% grid
[dn, z] = meshgrid(dn_mean, z_mean);


%%
% test figure ? 

cs = {tempG', salinG', densG'};  % Transpose and store in cell array
tt = {'Temperature', 'Salinity', 'Density'};  % Cell array of titles

for i = 1:length(cs)
    figure;
    pcolor(dn, z, cs{i});
    shading interp;
    colorbar;
    set(gca, 'ydir', 'reverse');
    datetick('x', 'keeplimits');
    title(sprintf('Salacia September 2021 %s', tt{i}));
    xlabel('Time');
    ylabel('Depth (m)');
end


%% save

ofn = fullfile(outdir, sprintf('%s_%s_L2.mat', config.glider_name, config.date_tag));
