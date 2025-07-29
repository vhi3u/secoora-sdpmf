%{
SECOORA Glider Processing L1 (G1 GLIDERS) - Mostly PEACH for now. 
Takes in L0 .ebdasc and .dbdasc files and converts into an L1 mat file. 
Basic QC, Thermal Lag Correction. 
Required files:
ebds.mat and dbds.mat files for respective deployment

different versions being created later. 

% Thermal lag correction only applies to G1 gliders 

% Current gliders supported: 

g3: angus franklin 1091
g2: modena
g1: ramses, pelagia, salacia, bass, sam

pumped: angus, franklin, 1091, modena
nonpumpted: ramses, pelagia, salacia, bass, sam

%}
clc
clear;

% Add glider start and end date. refer to SECOORA sheet (Jennifer Dorton)
strStartDate = '25-Sep-2021';
strEndDate = '02-Nov-2021';

gliderName = 'salacia';
dateTag = '2021_09';
dTag = 'September_21';

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

fprintf('Starting L1 CTD processing... \nGlider name: %s \nGen: %s \nisPumped = %i\n', gliderName, strGen, isPumped);

%%

% load .ebdasc and .dbdasc files 

ebddir = sprintf('/Users/victornguyen/Documents/MATLAB/dataSECOORA/level0/%s/%s',gliderName, dateTag);
dbddir = sprintf('/Users/victornguyen/Documents/MATLAB/dataSECOORA/level0/%s/%s',gliderName, dateTag);

outdir = sprintf('/Users/victornguyen/Documents/MATLAB/dataSECOORA/level1/%s/', gliderName);

ebd = fullfile(ebddir, sprintf('%s_%s_allebds.mat', dTag, gliderName));
dbd = fullfile(dbddir, sprintf('%s_%s_alldbds.mat', dTag, gliderName));

% use ebd.mat and dbds.mat instead of all*bds when we add angle of attack
% for velocity...

load(ebd); % sstruct : science structure
load(dbd); % fstruct : flight structure

% prepare arrays for L0

% DBD
altitude = [];
ptime_dbd =[];
horizontalVelocity =[];
depth = [];
pitch=[];
avgDepthRate = [];
angleOfAttack = [];
gpsLat = [];
gpsLon = [];
num_inflections = [];

% EBD
temp = [];
cond =[];
pres =[];
ctd_time = [];
ptime_ebd =[];
% might want to add the ECO stuff later? 

% fill in DBD

ptime_dbd = fstruct.m_present_time;
horizontalVelocity = fstruct.m_speed;
altitude = fstruct.m_altitude;
depth = fstruct.m_depth;
pitch = fstruct.m_pitch;
avgDepthRate = fstruct.m_avg_depth_rate;
angleOfAttack = fstruct.u_angle_of_attack;
gpsLat = fstruct.m_gps_lat;
gpsLon = fstruct.m_gps_lon;

num_inflections = fstruct.m_tot_num_inflections;

% fill in EBD

temp = sstruct.sci_water_temp;
cond = sstruct.sci_water_cond;
pres = sstruct.sci_water_pressure;
ctd_time = sstruct.sci_ctd41cp_timestamp; % not sure if this is the right timestamp?
ptime_ebd = sstruct.sci_m_present_time;

%% L0 Quality Control

% monotonic time sort (flight files)

[Y,I] = sort(ptime_dbd);
ptime_dbd = Y;
horizontalVelocity = horizontalVelocity(I); length(horizontalVelocity)
depth = depth(I);
altitude = altitude(I);
pitch = pitch(I);
avgDepthRate = avgDepthRate(I);
angleOfAttack = angleOfAttack(I);

gpsLat = gpsLat(I);
gpsLon = gpsLon(I);

% monotonic time sort (science files)

[Y,I] = sort(ptime_ebd);
ptime_ebd = Y;
temp = temp(I);
cond = cond(I);
pres = pres(I);
ctd_time=ctd_time(I);

% remove disrepancies between ptime_ebd and ctd_time

iweird = find(ptime_ebd-ctd_time > 10);
ptime_ebd(iweird)=NaN; 
temp(iweird)=NaN; 
pres(iweird)=NaN; 
cond(iweird)=NaN;
ctd_time(iweird)=NaN;

% remove duplicate timestamps in flight files 
[unique_dbds, IA] = unique(ptime_dbd);
ptime_dbd = unique_dbds;
horizontalVelocity = horizontalVelocity(IA);
depth = depth(IA);
altitude = altitude(IA);
pitch = pitch(IA);
avgDepthRate = avgDepthRate(IA);
angleOfAttack = angleOfAttack(IA);
gpsLat = gpsLat(IA);
gpsLon = gpsLon(IA);
num_inflections = num_inflections(IA);

% remove nans from EBD data...
i = find(~isnan(temp) & ~isnan(pres) & ~isnan(cond));
ptime_ebd = ptime_ebd(i); temp = temp(i); cond = cond(i); pres = pres(i); 
ctd_time = ctd_time(i);

% remove conductivity values less than 1, must be at surface or bad
i = find(cond>=1);
ptime_ebd = ptime_ebd(i);  temp = temp(i);  pres = pres(i);  cond = cond(i); 
ctd_time = ctd_time(i);

% remove pressure values less than 0
i = find(pres>=0);
ptime_ebd = ptime_ebd(i);  temp = temp(i);  pres = pres(i);  cond = cond(i); 
ctd_time = ctd_time(i);

%% Thermal Lag Prep

% convert pitch and angle of attack into degrees
pitch = pitch*180/pi;
angleOfAttack = angleOfAttack*180/pi;

% compute glide angle
glideAngle = pitch + angleOfAttack;

% find glider profiles
[profile_index, profile_direction] = findProfiles(depth, 'stall', 1.5, ...
						  'inversion', 1.5, 'interrupt', inf);
num_profiles = floor(max(profile_index));
num_yos = floor(num_profiles);



% create arrays to store parameters (alpha_s, alpha_o, tau_s, tau_o)
params=[]; fvals=[]; residuals=[]; profiles=[]; pidx=[];


%% Science Corrections

% make copy of dbd time stamp vector for use in salinity/density correction...
ptime1_dbd = ptime_dbd;
% % test for duplicates
% [~, counts] = unique(ptime1_dbd(i));
% sum(counts)  % should equal length(ptime1_dbd(i)) if no duplicates
% length(ptime1_dbd(i)) - length(unique(ptime1_dbd(i)))
% 
% 
% duplicated_indices = find(histcounts(ptime1_dbd(i), unique(ptime1_dbd(i))) > 1);

% remove nans from dbd
i = find(~isnan(horizontalVelocity)); 
hv = interp1(ptime1_dbd(i), horizontalVelocity(i), ctd_time); 

i = find(~isnan(horizontalVelocity)&(horizontalVelocity>0.1 & horizontalVelocity<0.6));
horizontalVelocity = interp1(ptime1_dbd(i), horizontalVelocity(i), ctd_time);

i = find(~isnan(depth));
depth = interp1(ptime1_dbd(i), depth(i), ctd_time);

i = find(~isnan(altitude));
altitude = interp1(ptime1_dbd(i), altitude(i), ctd_time);

i = find(~isnan(pitch));
pitch = interp1(ptime1_dbd(i), pitch(i), ctd_time);

i = find(~isnan(angleOfAttack));
angleOfAttack = interp1(ptime1_dbd(i), angleOfAttack(i), ctd_time);

i = find(~isnan(avgDepthRate));
avgDepthRate = interp1(ptime1_dbd(i), avgDepthRate(i), ctd_time);

i = find(~isnan(glideAngle));
glideAngle = interp1(ptime1_dbd(i), glideAngle(i), ctd_time);

% last set of naning
i = find(~isnan(horizontalVelocity));
horizontalVelocity = horizontalVelocity(i); depth = depth(i); altitude = altitude(i); pitch = pitch(i);
avgDepthRate = avgDepthRate(i); glideAngle = glideAngle(i); 
ptime_ebd = ptime_ebd(i); temp = temp(i); cond = cond(i); pres = pres(i); 
ctd_time = ctd_time(i); hv = hv(i);

% scale up pressure
pres = pres*10;

% calculate salinity (without correction)...
salin = sw_salt(10*cond/sw_c3515, temp, pres);

% calculate density (without correction)...
dens = sw_pden(salin, temp, pres, 0);

% calculate glider velocity
gliderVelocity = sqrt(hv.^2 + avgDepthRate.^2);

% create datenum 
ptime_datenum = ctd_time/3600/24+datenum(1970, 1, 1, 0, 0, 0);

ptime_dbd_gps = ptime_dbd;

% convert lat and lon to digital degrees 
% UPDATE: Apparently this isn't needed... 
% gpsLat = ddmm2decdeg(gpsLat);
% gpsLon = ddmm2decdeg(gpsLon);

% eliminate outliers in gpsLat, gpsLon...
i = find(abs(gpsLat) <= 90.0);
gpsLat = gpsLat(i);  gpsLon = gpsLon(i);  ptime_dbd_gps = ptime_dbd_gps(i);
i = find(abs(gpsLon) <= 180.0);
gpsLat = gpsLat(i);  gpsLon = gpsLon(i);  ptime_dbd_gps = ptime_dbd_gps(i);

% remove NaNs in lat/lon
i = find(~isnan(gpsLat));
gpsLat = gpsLat(i);  gpsLon = gpsLon(i);  ptime_dbd_gps = ptime_dbd_gps(i);

% interpolate dbd data so it aligns with ebd lat/lon
gpsLat = interp1(ptime_dbd_gps, gpsLat, ctd_time);
gpsLon = interp1(ptime_dbd_gps, gpsLon, ctd_time);

depth_dbd = depth;

depth = sw_dpth(pres, gpsLat);

%% Thermal Lag Parameters (only for G1)

if ismember(strGen, 'G1')

    tic;
    np = 1; 
    
    % check to see if the first profile index is positive 1, if not, increment
    % to next one. 
    if ~(profile_index == np & profile_direction == 1)
        np = np + 1;
    end
    
    
    f = waitbar(0, 'Starting');
    
    % loop through and get thermal lag params for each profile (Sara Haines)
    while (np<num_profiles)
    
      waitbar(np/num_profiles, f, sprintf('Progress: %d %%', floor(np/num_profiles*100)));
    
      % ensure downcast followed by upcast
      subset1 = find(profile_index == np & profile_direction == 1);
      subset2 = find(profile_index == np+1 & profile_direction == -1);
      if (length(subset1)>=20 & length(subset2)>=20)
        time1 = ctd_time(subset1);
        time2 = ctd_time(subset2);
        t1 = ptime_datenum(subset1);
        t2 = ptime_datenum(subset2);
        cond1 = cond(subset1);
        cond2 = cond(subset2);
        temp1 = temp(subset1);
        temp2 = temp(subset2);
        pres1 = pres(subset1);
        pres2 = pres(subset2);
        flow1 = gliderVelocity(subset1);
        flow2 = gliderVelocity(subset2);
        
    
        % take out flow velocity values less than 1 cm/s to ensure stability of
        % findThermalLagParams (denominator close to zero)
        idx1 = find(flow1 < 0.01);
        flow1(idx1) = NaN;
        idx2 = find(flow2 < 0.01);
        flow2(idx2) = NaN;
        
        % knock down glider velocity to flow velocity using Morison relationship 
        % using polyfit that ensures cell flow goes to 0 when gv=0
        pval = [0.2906    0.1876   -0.0008];
        flow1 = polyval(pval, flow1);
        flow2 = polyval(pval, flow2);
        
        try
          warning off;
          % with graphics -- debug or analyze data issues
          % [param_set, exitflag, residual] = findThermalLagParams(time1,cond1,temp1,pres1,flow1,time2,cond2,temp2,pres2,flow2,'graphics','true','upper',[2 1 20 10]);
          [param_set, exitflag, residual] = findThermalLagParams(time1,cond1,temp1,pres1,flow1,time2,cond2,temp2,pres2,flow2,'graphics',false);
          % without graphs -- faster
          % [param_set, exitflag, residual] = findThermalLagParams(time1,cond1,temp1,pres1,flow1,time2,cond2,temp2,pres2,flow2);
          warning on;
          % pause(2)
          % close
        catch
          exitflag = -1;
        end
        if exitflag > 0
          params = [params; param_set];
          residuals = [residuals; residual];
          profiles = [profiles; np np+1];
          pidx = [pidx; subset2(1)]; % idx of start of upcast (middle of yo)
        end % if exitflag okay
      else
        fprintf('profile %d or %d less than 20 samples\n', np, np+1); 
      end % if subset lengths good
      np = np+2;
    end % while
    close(f)
    toc;
end

%% Store Parameters

% FIX INTO LOOP ABOVE 

% trim outliers
bounds = zeros(4,2);
for i = 1:4
    bounds(i,1) = prctile(params(:,i), 1);
    bounds(i,2) = prctile(params(:,i), 99);

end

valid_i = find(...
    params(:,1) >= bounds(1,1) & params(:,1) <= bounds(1,2) & ...
    params(:,2) >= bounds(2,1) & params(:,2) <= bounds(2,2) & ...
    params(:,3) >= bounds(3,1) & params(:,3) <= bounds(3,2) & ...
    params(:,4) >= bounds(4,1) & params(:,4) <= bounds(4,2));

params_trim = params(valid_i, :);

fprintf('Total profiles: %d\n', size(params,1));
fprintf('Profiles after trimming: %d\n', numel(valid_i));
fprintf('Profiles removed: %d\n', size(params,1) - numel(valid_i));

v_params = [];
v_params.alphao = params_trim(:, 1);
v_params.alphas = params_trim(:, 2);
v_params.tauo = params_trim(:, 3);
v_params.taus = params_trim(:, 4);

v1_params = [];
v1_params.alphao = params(:, 1);
v1_params.alphas = params(:, 2);
v1_params.tauo = params(:, 3);
v1_params.taus = params(:, 4);

%% calculate mean and median of thermal lag parameters

% ALSO FIX INTO LOOP ABOVE

thermal_lag_medians = [];
thermal_lag_medians(1, 1) = median(v1_params.alphao);
thermal_lag_medians(2, 1) = median(v1_params.alphas);
thermal_lag_medians(3, 1) = median(v1_params.tauo);
thermal_lag_medians(4, 1) = median(v1_params.taus);

thermal_lag_means = [];
thermal_lag_means(1, 1) = mean(v_params.alphao);
thermal_lag_means(2, 1) = mean(v_params.alphas);
thermal_lag_means(3, 1) = mean(v_params.tauo);
thermal_lag_means(4, 1) = mean(v_params.taus);

%% Apply thermal lag correction

% NaN out stalled profiles (where profile_index mod 2 != 0)

% FIX INTO LOOP ABOVE

valid = ~isnan(profile_index) & ~isnan(ptime1_dbd);
profile_index_int = interp1(ptime1_dbd(valid), profile_index(valid), ctd_time, 'nearest', NaN);

istall = abs(mod(profile_index_int, 1) - 0.5) <eps;

temp(istall) = NaN;

[temp_inside, cond_outside] = correctThermalLag(ctd_time,cond,temp,gliderVelocity,thermal_lag_medians);

tempCorrected = temp_inside;

salinCorrected = sw_salt(10*cond/sw_c3515, tempCorrected, pres);

% repopulate salinCorrected with raw salin values at ip index
salinCorrected(istall) = salin(istall); 

% add density
densCorrected = sw_pden(salinCorrected, tempCorrected, pres, 0);

% plot
figure;

ts1 = ccplot(salin, temp, ctd_time, [ctd_time(1) ctd_time(end)], '.',10);
colorbar;
clim([ctd_time(1) ctd_time(end)]);
xlabel('salinity (psu)')
ylabel('temp (c)');
title('original');
figure;
ts2 = ccplot(salinCorrected, tempCorrected, ctd_time, [ctd_time(1) ctd_time(end)], '.',10);
colorbar;
clim([ctd_time(1) ctd_time(end)]);
xlabel('salinity (psu)')
ylabel('temp (c)');
title('corrected');

figure;
tss1 = plot(salin, temp, '.', 'DisplayName','original');
hold on
tss2 = plot(salinCorrected, tempCorrected, '.', 'DisplayName','corrected');
hold off
xlabel('salinity (psu)')
ylabel('temp (c)');
title('original vs corrected (no time factor)');
legend;

%% Save
ofn = fullfile(outdir, sprintf('%s_%s_CTD_L1.mat', gliderName, dateTag));

 % 'ptime', 'seconds since 0000-01-01T00:00',...
%  'ptime', 'time vector reported by glider',...
%   'ptime',...

units = struct('altitude', 'm',...
               'angleOfAttack', 'decimal degrees',...
               'avgDepthRate', 'm/s',...
               'depth', 'm',...
               'horizontalVelocity', 'm/s',...
               'pitch', 'decimal degrees',...
               'ptime_datenum', 'days since 1970-01-01T00:00', ...
               'temp', 'deg C', ...
               'tempCorrected', 'deg C', ...
               'salin', 'psu', ...
               'salinCorrected', 'psu', ...
               'dens', 'kg m-3', ...
               'densCorrected', 'kg m-3');

variable_description = struct('altitude', 'altimeter measured distance from bottom',...
                              'angleOfAttack', 'difference between pitch and glider angle',...
                              'avgDepthRate', 'average rate of change of depth, >0 is down',...
                              'depth', 'depth calculated as function of pressure and position latitude',...
                              'horizontalVelocity', 'vehicle horizontal speed through water',...
                              'pitch', 'vehicle angle of inclination, >0 is nose up',...
                              'ptime_datenum', 'Serial Date Number string', ...
                              'temp', 'temperature measured from CTD', ...
                              'tempCorrected', 'temperature from thermal lag correction (Garau)', ...
                              'salin', 'salinity measured from calculation of temp and cond', ...
                              'salinCorrected', 'adjusted salinity with thermal lag correction', ...
                              'thermal_params', 'thermal lag coefficients based on Morison (1994)', ...
                              'dens', 'density measured', ...
                              'densCorrected', 'density corrected using salinCorrected');

config = struct('glider_name', gliderName, ...
                'gen', strGen, ...
                'isPumped', isPumped, ...
                'date_tag', dateTag, ... 
                'start_date', strStartDate, ...
                'end_date', strEndDate, ...
                'var_description', variable_description, ...
                'var_units', units, ...
                'thermal_params', thermal_lag_medians, ...
                'stall_index', istall);



save(ofn,...
     'config',...
     'altitude',...
     'angleOfAttack',...
     'avgDepthRate',...
     'depth',...
     'horizontalVelocity',...
     'pitch',...
     'ptime_datenum', ...
     'temp', ...
     'tempCorrected', ...
     'salin', ...
     'salinCorrected', ...
     'dens', ...
     'densCorrected');

