%{
SECOORA Glider Processing ECO L1 - Mostly PEACH for now. 
Takes in L0 .ebdasc and .dbdasc files and converts into an L1 mat file. 
Required files:
ebds.mat and dbds.mat files for respective deployment

Minmax filter for chlorophyll fluorescence, and 7 point median filter for
CDOM. 

% Current gliders supported: 

g3: angus franklin 1091
g2: modena
g1: ramses, pelagia, salacia, bass, sam

pumped: angus, franklin, 1091, modena
nonpumpted: ramses, pelagia, salacia, bass, sam

%}

% Add glider start and end date. refer to SECOORA sheet (Jennifer Dorton)
strStartDate = '11-Nov-2021';
strEndDate = '01-Dec-2021';

gliderName = 'angus';
dateTag = '2021_11';
dTag = 'November_21';

if ismember(gliderName, {'ramses', 'pelagia', 'salacia', 'bass', 'sam'})
    strGen = 'G1';
elseif ismember(gliderName, {'modena'})
    strGen = 'G2';
elseif ismember(gliderName, {'angus', 'franklin', 'unit_1091'}) 
    strGen = 'G3';
else
    error('%s is not a valid glider name', gliderName);
end

if ismember(strGen, {'G1'})
    isPumped = false;
else
    isPumped = true;
end

fprintf('Starting L1 ECO processing... \nGlider name: %s \nGen: %s \nisPumped = %i\n', gliderName, strGen, isPumped);

% LOAD DATA

ebddir = sprintf('/Users/victornguyen/Documents/MATLAB/dataSECOORA/level0/%s/%s',gliderName, dateTag);
dbddir = sprintf('/Users/victornguyen/Documents/MATLAB/dataSECOORA/level0/%s/%s',gliderName, dateTag);

outdir = sprintf('/Users/victornguyen/Documents/MATLAB/dataSECOORA/level1/%s/', gliderName);

ebd = fullfile(ebddir, sprintf('%s_%s_allebds.mat', dTag, gliderName));
dbd = fullfile(dbddir, sprintf('%s_%s_alldbds.mat', dTag, gliderName));

% use ebd.mat and dbds.mat instead of all*bds when we add angle of attack
% for velocity...

load(ebd); % sstruct : science structure
load(dbd); % fstruct : flight structure 

% DBD
ptime_dbd = fstruct.m_present_time;
horizontalVelocity = fstruct.m_speed;
altitude = fstruct.m_altitude;
depth = fstruct.m_depth;
pitch = fstruct.m_pitch;
avgDepthRate = fstruct.m_avg_depth_rate;
angleOfAttack = fstruct.u_angle_of_attack;
gpsLat = fstruct.m_gps_lat;
gpsLon = fstruct.m_gps_lon;

% EBD

pres = sstruct.sci_water_pressure;
oxyw_oxygen = 

% L0 BASIC QUALITY CONTROL CHECK (run once due to sort)


% sort ECO
[Y,I] = sort(ptime_ebd);
ptime_ebd = Y;
pres = pres(I);
chlor=chlor(I);
cdom_raw=cdom_raw(I);
scatter=scatter(I);
chlor_sig=chlor_sig(I);
cdom_sig=cdom_sig(I);
scatter_sig=scatter_sig(I);
chlor_ref=chlor_ref(I);
cdom_ref=cdom_ref(I);
scatter_ref=scatter_ref(I);
eco_time=eco_time(I);

[ptime_ebd, IA] = unique(ptime_ebd, 'stable');
pres = pres(IA);
chlor = chlor(IA);
cdom_raw = cdom_raw(IA);
scatter = scatter(IA);
chlor_sig = chlor_sig(IA);
cdom_sig = cdom_sig(IA);
scatter_sig = scatter_sig(IA);
chlor_ref = chlor_ref(IA);
cdom_ref = cdom_ref(IA);
scatter_ref = scatter_ref(IA); 
eco_time = eco_time(IA);

% use ptime_ebd as timebase
ptime = ptime_ebd;

% remove triple zeros
tri0 = find(chlor > 0 & cdom_raw > 0 & scatter > 0);
ptime = ptime(tri0);
pres = pres(tri0);
chlor = chlor(tri0);
cdom_raw = cdom_raw(tri0);
scatter = scatter(tri0);
chlor_sig = chlor_sig(tri0);
cdom_sig = cdom_sig(tri0);
scatter_sig = scatter_sig(tri0);
chlor_ref = chlor_ref(tri0);
cdom_ref = cdom_ref(tri0);
scatter_ref = scatter_ref(tri0);

% ptime to datenum format
ptime_datenum = ptime / 3600 / 24 + datenum(1970, 1, 1, 0, 0, 0);

% sort dbd
[ts, isort] = sort(ptime_dbd);
ptime_dbd = ts;
depth = depth(isort);
gpsLat = gpsLat(isort);
gpsLon = gpsLon(isort);

i = find(abs(gpsLat) <= 90.0);
gpsLat = gpsLat(i);  gpsLon = gpsLon(i);  ptime_dbd = ptime_dbd(i);
i = find(abs(gpsLon) <= 180.0);
gpsLat = gpsLat(i);  gpsLon = gpsLon(i);  ptime_dbd = ptime_dbd(i);

[ptime_dbd, Iunique] = unique(ptime_dbd, 'stable');
gpsLat = gpsLat(Iunique);
gpsLon = gpsLon(Iunique);
depth = depth(Iunique);  % optional: if used later for something else

gpsLat = interp1(ptime_dbd, gpsLat, ptime);
gpsLon = interp1(ptime_dbd, gpsLon, ptime);

depth = sw_dpth(pres*10, gpsLat);

% Minmax filter for chlor
% Estimate lower envelope of bb using two-stage minimum filter
% to isolate transient spikes from background scattering

ichl = 11;
xchl = chlor;
chlor_npts = ichl;
chlor_nside = floor(chlor_npts / 2);
chlor_notnans = find(~isnan(xchl));
chlor_min = nan(size(xchl));    
chlor_maxmin = nan(size(xchl));   
chlor_min(1+chlor_nside:end-chlor_nside) = minmaxfilt1(xchl(chlor_notnans), chlor_npts, 'min', 'valid');
chlor_maxmin(1+chlor_nside:end-chlor_nside) = minmaxfilt1(chlor_min, chlor_npts, 'min', 'valid');

spikes_chlor = xchl - chlor_maxmin;
chlor_small = chlor_maxmin;
chlor_large = spikes_chlor;

% 7 pt median filter for cdom
cdom = medfilt1(cdom_raw, 7);

% minmax filter for backscatter
ibb = 11;
xbb = scatter;
scatter_npts = ibb;
scatter_nside = floor(scatter_npts / 2);
notnansbb = find(~isnan(xbb));
scatter_min = nan(size(xbb));
scatter_maxmin = nan(size(xbb));
scatter_min(1+scatter_nside:end-scatter_nside) = minmaxfilt1(xbb(notnansbb), scatter_npts, 'min', 'valid');
scatter_maxmin(1+scatter_nside:end-scatter_nside) = minmaxfilt1(scatter_min, scatter_npts, 'min', 'valid');

spikes_scatter = xbb - scatter_maxmin;
scatter_small = scatter_maxmin;
scatter_large = spikes_scatter;

% save variables

ofn = fullfile(outdir, sprintf('%s_%s_ECO_L1.mat', gliderName, dateTag));


units = struct('gpsLat', 'decimal degrees',...
               'gpsLon', 'decimal degrees',...
               'ptime_ebd', 'seconds since 0000-01-01T00:00',...
               'ptime_datenum', 'days since 1970-01-01T00:00',...
               'depth', 'm',...
               'chlor', '10e-6 g/l',...
               'chlor_small', '10e-6 g/l',...
               'chlor_large', '10e-6 g/l', ...
               'cdom', 'ppb',...
               'scatter', '1/mm',...
               'scatter_small', '1/mm',...
               'scatter_large', '1/mm',...
               'chlor_sig', 'count',...
               'cdom_sig', 'count',...
               'scatter_sig', 'count',...
               'chlor_ref', 'count',...
               'cdom_ref', 'count',...
               'scatter_ref', 'count');

variable_description = struct('depth', 'depth calculated as function of pressure and position latitude', ...
                              'gpsLat', 'position latitude measured by glider GPS',...
                              'gpsLon', 'position longitude measured by glider GPS',... 
                              'ptime', 'ptime',...
                              'ptime_datenum', 'serial date number string',...
                              'chlor', 'chlorophyll',...
                              'chlor_small','chlorophyll after QCed', ...
                              'chlor_large', 'chlorophyll after QCed',...
                              'cdom', 'chromophoric dissolved organic matter',...
                              'scatter', 'backscatter',...
                              'scatter_small', 'scatter after QCed',...
                              'scatter_large', 'scatter after QCed',...
                              'chlor_sig','chlorophyll raw count',...
                              'cdom_sig', 'chromophoric dissolved organic matter raw count',...
                              'scatter_sig', 'backscatter raw count',...
                              'chlor_ref','chlorophyll reference count',...
                              'cdom_ref', 'chromophoric dissolved organic matter reference count',...
                              'scatter_ref', 'backscatter reference count');
config = struct('glider_name', gliderName, ...
                'gen', strGen, ...
                'isPumped', isPumped, ...
                'date_tag', dateTag, ... 
                'start_date', strStartDate, ...
                'end_date', strEndDate, ...
                'var_description', variable_description, ...
                'var_units', units);

save(ofn,...
    'config',...
    'gpsLat',...
    'gpsLon',...
    'ptime_ebd',...
    'ptime_datenum',...
    'depth',...
    'chlor',...
    'chlor_small',...
    'chlor_large',...
    'cdom',...
    'scatter',...
    'scatter_small',...
    'scatter_large',...
    'chlor_sig',...
    'cdom_sig',...
    'scatter_sig',...
    'chlor_ref',...
    'cdom_ref',...
    'scatter_ref');

% diagnostic plots

figure
plot(ptime_datenum, cdom, 'DisplayName','cdom');
hold on
plot(ptime_datenum, cdom_raw, '.', 'DisplayName','raw cdom');
datetick;
xlabel('date')
ylabel('cdom')
title('cdom median filter check');

figure
plot(ptime_datenum, chlor_small, 'DisplayName','chlor (small)');
hold on
plot(ptime_datenum, chlor_large, 'DisplayName','chlor (large');
datetick;
xlabel('date');
ylabel('chlor')
title('chlor minmax filter check');

figure
plot(ptime_datenum, cdom_raw, '.', 'DisplayName', 'Raw CDOM', 'Color', [0.7 0.7 0.7]) % raw points
hold on
plot(ptime_datenum, cdom, 'r', 'LineWidth', 1.5, 'DisplayName', 'Median Filtered CDOM') % filtered line

datetick('x', 'dd-mmm-yyyy')
xlabel('Date')
ylabel('CDOM (ppb)')
title('Comparison of Raw vs Median Filtered CDOM')
legend('Location', 'best')
grid on