%{
SECOORA Glider Processing DO L1 - Mostly PEACH for now.  
Required files:
ebds.mat and dbds.mat files for respective deployment

Dissolved oxygen time lag for recalibration. 

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

fprintf('Starting L1 DO processing... \nGlider name: %s \nGen: %s \nisPumped = %i\n', gliderName, strGen, isPumped);

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
oxyw_oxygen = sstruct.sci_oxy4_oxygen;
oxyw_saturation = sstruct.sci_oxy4_saturation;
oxyw_temp = sstruct.sci_oxy4_temp;
oxyw_dphase = sstruct.sci_oxy4_calphase;