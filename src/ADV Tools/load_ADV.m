
%% Load ADV data and populate return structure
function [anc, config, data, time] = load_ADV(file_path, file_name, config)

temp = load([file_path file_name '.dat']);

%extract data from temporary structure

data.rvel_x = temp(:,3); %raw (non-despiked) velocity (m/s)
data.rvel_y = temp(:,4);
data.rvel_z = temp(:,5);

anc.amp1 = temp(:,6);   %echo amplitude (counts)
anc.amp2 = temp(:,7);
anc.amp3 = temp(:,8);

anc.snr1 = temp(:,9);   %signal to noise ratio (SNR)
anc.snr2 = temp(:,10);
anc.snr3 = temp(:,11);

anc.corr1 = temp(:,12); %correlation
anc.corr2 = temp(:,13);
anc.corr3 = temp(:,14);

data.p = temp(:,15);    %pressure

%load header data
temp = importdata([file_path file_name '.hdr'],'\t');

%extract sample rate
config.sample_rate = temp{12,1};
config.sample_rate = strrep(config.sample_rate,'Sampling rate                         ','');
config.sample_rate = strrep(config.sample_rate,' Hz','');
config.sample_rate = str2double(config.sample_rate);

%extract start time
config.t_start = temp{7,1};
config.t_start = strrep(config.t_start,'Time of first measurement             ','');
config.t_start = datenum(config.t_start);

%generate time stamps
time.t = [0:1/config.sample_rate:(length(data.p)-1)/config.sample_rate]';   %seconds
time.date = config.t_start/(24*3600) + time.t;                              %days

end
