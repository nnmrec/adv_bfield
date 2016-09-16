%Brian Polagye
%July 13, 2013

%Description: routine to process ADV data from flume tests

function batch_flume(parameters)

%%Unwrap parameters

file_name = parameters.file_name;
file_path = parameters.file_path;

config.x_pos = parameters.x_pos;
config.y_pos = parameters.y_pos;
config.z_pos = parameters.z_pos;

config.description = parameters.description;

%load ADV data and populate main structure
[anc, config, data, time] = load_ADV(file_path, file_name, config);

%rotate x,y,z data onto principal axis

%despike velocity data
[data.vel_x, data.vel_y, data.vel_z, ip] = ...
    func_despike_phasespace3d_3var(data.rvel_x, data.rvel_y, data.rvel_z, 2);
anc.num_spikes = length(ip)/length(data.p);
anc.spike_pts = ip;

%calculate velocity magnitude and perturbation
[data.s, data.smean, data.sprime] = calculate_turbulence(data.vel_x, data.vel_y, data.vel_z);

data.I = std(data.sprime)/data.smean;

%generate velocity spectra
[spectra.Pxx, spectra.Pyy, spectra.Pzz, spectra.f] = ... 
    calculate_spectra(data.vel_x, data.vel_y, data.vel_z, time.t, parameters.overlap);

%save data
save([file_path file_name '.mat'],'anc','config','data','time','spectra')

end

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

%% Calculate turbulence in flow
function [s, smean, sprime] = calculate_turbulence(vel_x, vel_y, vel_z)

%velocity magnitude
s = (vel_x.^2+vel_y.^2+vel_z.^2).^(0.5);

%mean magnitude
smean = mean(s);

%magnitude perturbation
sprime = s - smean;

end

%% Calculate velocity spectra
function [Pxx, Pyy, Pzz, f] = calculate_spectra(vel_x, vel_y, vel_z, t, overlap)

debug=1;

%calculate number of points per overlapping window
num_points = 2^(floor(log2(length(vel_x)))-3);    %include 6 DOF per sample

[Pxx, f] = speed_fft(vel_x,t,num_points,overlap);
[Pyy, ~] = speed_fft(vel_y,t,num_points,overlap);
[Pzz, ~] = speed_fft(vel_z,t,num_points,overlap);

if debug==1
   figure(1)
   clf
   
   loglog(f,Pxx,'-b','linewidth',2)
   hold on
   loglog(f,Pyy,'-r','linewidth',2)
   loglog(f,Pzz,'-g','linewidth',2)
   grid on
   xlabel('Frequency (Hz)','fontweight','b')
   ylabel('P (m^s/s^2/Hz)','fontweight','b')
   
   legend('U_x','U_y','U_z')
   
end

end