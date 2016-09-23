%% ADV tools by Ben Strom ... later modded by Danny Sale (8/23/2016)
% this code was written to post-process the Vectrino ADV on the gantry
% system, in Bamfield flume 2015
%
% changes: Danny adds the estimation of turbulent length scale,
%          formats some figures, adds a sample dataset and testcase,
%          import the corresponding LabView data and STC with the ADV
%
% to-do: add surface plots of the u, rms, TI, and L_turb fields
% to-do: setup autosave of figures to files
% to-do: plot the positions of the ADVs relative to flume walls
% to-do: return only part of the times series, choose longest time between "de-spikes"
% 
% bugs: there is a re-computation of the time variable , but this only version that works
%
%
%% STARTUP
clear all
clc
close all

% this is path to misc Matlab utilities
addpath( genpath([pwd filesep 'src']) )
% this is the path to the file 'calc_hwa_corrspec.m'
addpath( genpath('/mnt/data-RAID-1/danny/FrenchNavalAcademy/PIV_Processing/') )

% NOTE: to suppress all figure, start matlab like:
% matlab -nodesktop -nodisplay
% or try 
% set(0,'defaultFigureVisible','off')

%% USER INPUTS
dat_folder = [pwd filesep 'test_cases/3d_transects_no_turbine_fs1p2/x_0'];
timeID = {'19407848'};
% timeID = {'19407848', ...
%           '19407903', ...
%           '19407951', ...
%           '19407999', ...
%           '19408046', ...
%           '19408114', ...
%           '19408161', ...
%           '19408209', ...
%           '19408257', ...
%           '19408306', ...
%           '19408354', ...
%           '19408402', ...
%           '19408450', ...
%           '19408498', ...
%           '19408545', ...
%           '19408594', ...
%           '19408641', ...
%           '19408689', ...
%           '19408738', ...
%           '19408786', ...
%           '19408834', ...
%           '19408882', ...
%           '19408931', ...
%           '19408979', ...
%           '19409026', ...
%           '19409075', ...
%           '19409122', ...
%           '19409170', ...
%           '19409218', ...
%           '19409266', ...
%           '19409314', ...
%           '19409363', ...
%           '19409411', ...
%           '19409459', ...
%           '19409506'};
      
%% this reads the ADV hardware config and ADV measurements, then applies filtering and de-spiking, returning your "cleaned" ADV data    
for n = 1:numel(timeID)
    % read the ADV data file
    [v0(n), v1(n), sRate(n)] = vectrinoCalcPoint(dat_folder, timeID{n});
    
    
    
    % danny: attempt to compute turbulent length scales
    % consider this part of the time series as correlated, 
	% you may want to only use a segment of the signal that did not contain any "de-spikes"
%     Li     = 1:50:500;
    n0 = 400;       % start sample of 'good data'
    n1 = 700;       % end sample of 'good data'
    fi = 50;        % number of bins in integration
%     Li     = 400:50:700;  
%     Li = n0:fi:n1;
    
    
    % Li should start at 1

    % use only the 'good part' of the signal
    vel_x = v0.U(n0:n1,1);
%     vel_y = v0(n).U(n0:n1,2);
%     vel_z = v0(n).U(n0:n1,3);
    % reconstruct time
    dt = 1/sRate;
    t  = [0:dt:dt*(numel(vel_x)-1)]';
    
    Li = 1:fi:numel(vel_x);
    
%     L      = calc_hwa_corrspec(v0.U(:,1), t, Li, true);
%     L      = calc_hwa_corrspec(v0.U(n0:n1,1), t(n0:n1), Li, true);
    Lt      = calc_hwa_corrspec(vel_x, t, Li, true);
    L_turb(n) = Lt(end);

end

%% plot the positions of the ADVs
figure
hold on
for n = 1:numel(timeID)
    xi = v0(n).pos(1);
    yi = v0(n).pos(2);
    zi = v0(n).pos(3);
    plot3(xi, yi, zi, 'ok')
    
    xo = v1(n).pos(1);
    yo = v1(n).pos(2);
    zo = v1(n).pos(3);
    plot3(xo, yo, zo, 'or')
    
end

title('ADV positions')
xlabel('streamwise, x [m]')
ylabel('cross flow, y [m]')
zlabel('depth, z [m]')
box on; grid on;
axis vis3d


%%
n_samples   = numel(v0.U(:,1));
t_max       = n_samples / sRate;
t           = linspace(0, t_max, n_samples);

% speed time series
v0.spd = sqrt(v0.U(:,1).^2 + ...
              v0.U(:,2).^2 + ...
              v0.U(:,3).^2);
v1.spd = sqrt(v1.U(:,1).^2 + ...
              v1.U(:,2).^2 + ...
              v1.U(:,3).^2);
           
% v0.ux_mean = nanmean(v0.U(:,1));
% v0.uy_mean = nanmean(v0.U(:,2));
% v0.uz_mean = nanmean(v0.U(:,3));
% v0.u_mean  = mean([v0.ux_mean, v0.uy_mean, v0.uz_mean]);
v0.mean  = nanmean(v0.spd);
v0.std   = nanstd(v0.spd);
% v0.TI      = v0.u_std / v0.u_mean;
v0.TI    = nanstd(v0.spd - v0.mean) / v0.mean;
v0.skew  = skewness(v0.spd,1);

v0.u_avg   = nanmean(v0.U(:,1));
v0.v_avg   = nanmean(v0.U(:,2));
v0.w_avg   = nanmean(v0.U(:,3));
v0.u_prime = v0.U(:,1) - v0.u_avg;
v0.v_prime = v0.U(:,2) - v0.v_avg;
v0.w_prime = v0.U(:,3) - v0.w_avg;
v0.u_std   = nanstd(v0.u_prime);
v0.v_std   = nanstd(v0.v_prime);
v0.w_std   = nanstd(v0.w_prime);

% Nick wrote it this way
v0.rms     = sqrt((1/3)*(sum(v0.u_prime.^2)/length(v0.u_prime) + ...
                         sum(v0.v_prime.^2)/length(v0.v_prime) + ... 
                         sum(v0.w_prime.^2)/length(v0.w_prime)));
v0.ti      = v0.rms / sqrt(v0.u_avg^2+v0.v_avg^2+v0.w_avg^2);

% Nick wrote the RMS this way
% u_avg(i)=mean(u);
% u_prime=u-mean(u);
% v_avg(i)=mean(v);
% v_prime=v-mean(v);
% w_avg(i)=mean(w);
% w_prime=w-mean(w);
% 
% rms(i)=sqrt((1/3)*(sum(u_prime.^2)/length(u_prime)...
%     +sum(v_prime.^2)/length(v_prime)+sum(w_prime.^2)/length(w_prime)));
% TI(i)=rms(i)/(sqrt(u_avg(i)^2+v_avg(i)^2+w_avg(i)^2));



%% basic plot of the velocity components
figure()



subplot(4,1,1)
hold on
plot(t, v0.U(:,1), '-k');
plot(t, v1.U(:,1), '-r');
title(['ADV signal from gantry, timeID' timeID{1}] )
ylabel('u_x [m/s]')
box on; grid on;

set(findall(gca,'-property','FontSize'),'FontSize',14)
set(findall(gca,'-property','LineWidth'),'LineWidth',1)

subplot(4,1,2)
hold on
plot(t, v0.U(:,2), '-k');
plot(t, v1.U(:,2), '-r');
ylabel('u_y [m/s]')
box on; grid on;

set(findall(gca,'-property','FontSize'),'FontSize',14)
set(findall(gca,'-property','LineWidth'),'LineWidth',1)

subplot(4,1,3)
hold on
plot(t, v0.U(:,3), '-k');
plot(t, v1.U(:,3), '-r');
ylabel('u_z [m/s]')
box on; grid on;

set(findall(gca,'-property','FontSize'),'FontSize',14)
set(findall(gca,'-property','LineWidth'),'LineWidth',1)

subplot(4,1,4)
hold on
plot(t, v0.spd, '-k');
plot(t, v1.spd, '-r');
ylabel('speed [m/s]')
xlabel('time, t [s]')
legend('inboard', 'outboard')
box on; grid on;

set(findall(gca,'-property','FontSize'),'FontSize',14)
set(findall(gca,'-property','LineWidth'),'LineWidth',1)


figure()
hist(v0.spd - v0.mean, 500)
title('histogram of fluctuating velocity')
ylabel('count')
xlabel('fluctuating velocity, u'' [m/s]')
box on; grid on;

set(findall(gca,'-property','FontSize'),'FontSize',14)
set(findall(gca,'-property','LineWidth'),'LineWidth',1)


%% compute spectra
win_pts = 128;
overlap = 0.5;

% note: must deal with the NaNs before running this
% [Pxx, f, checksum] = speed_fft(v0.U(:,1), t, win_pts, overlap);
% [Pyy, f, checksum] = speed_fft(v0.U(:,2), t, win_pts, overlap);
% [Pzz, f, checksum] = speed_fft(v0.U(:,3), t, win_pts, overlap);
[Pxx, f, checksum] = speed_fft(v0.U(:,1), t, win_pts, overlap);
[Pyy, f, checksum] = speed_fft(v0.U(:,2), t, win_pts, overlap);
[Pzz, f, checksum] = speed_fft(v0.U(:,3), t, win_pts, overlap);

figure()
hold on
plot(f,Pxx,'-b','linewidth',1)
plot(f,Pyy,'-r','linewidth',1)
plot(f,Pzz,'-g','linewidth',1)
set(gca,'xscale','log');
set(gca,'yscale','log');
box on; grid on;
% legend('u_x','u_y','u_z')
legend('u''_x','u''_y','u''_z')
% xlabel('')
% ylabel('')

title('ADV speed, power spectra')
xlabel('freq. [Hz]')
ylabel('power spectral density [units of variance]')

set(findall(gca,'-property','FontSize'),'FontSize',14)
set(findall(gca,'-property','LineWidth'),'LineWidth',1)

x       = [0.5*max(f), 0.9*max(f)];
ynoll   = 0.2*max(Pxx);
y(1)    = ynoll;
y(2)    = y(1)*(x(2)/x(1))^(-5/3);
plot(x,y,'k-','linew',4)










%% compute turbulent length scales
% % 
% % this works, yes
% Li = 1:50:500;
% % Li = 400:50:700;
% L = calc_hwa_corrspec(v0.U(:,1), t, Li, true);
% 
% % this does not work??
% Li2 = 1:fi:numel(vel_x);
% L = calc_hwa_corrspec(v0.U(:,1), t, Li, true);
% 
% turb_length_scale = L(end);
% 













%% compute the average power available using definition in Neary et al.
r       = 0.225;        % rotor radius
u_ref   = 0.9;          % inflow speed
rho     = 1000;         % water density


