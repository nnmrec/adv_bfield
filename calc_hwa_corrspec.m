function L_turb = calc_hwa_corrspec(vx,t, Li, Nt, flag_plot)
%% calculate correlation / spectrum of HWA data
% the phuck method, modified by dsale
%
% questions for peter: in plotting, why [1:numel(cvxt)].*1e-4.*vxm
%                      instead of 1:numel(cvxt) ?
%                      why multiply by some scaling of the mean 1e-4.*vxm ?

%% converting data in python
%
% make sure to open python in the present working directory directory folder when saving
%
% python
% import numpy as np
% data = np.load('/Users/peterhuck/Desktop/HWA_Seattle/40M_center-Velocity.npy')
% np.savetxt("Numpytotxt.txt",data)
%
% if reading directly from a CSV file
% [labels,vx,y] = readColData('20M_vert-Velocity.txt',1);
%
%
%% filtering of the PIV signal
% need to remove any NANS
vx = naninterp(vx);
% maybe filter the RMS signal, remove outliers by +/- 3*STD

%% danny: look at some spectra
% returns the power spectrum E(W), where E is the
%   energy density and W the pulsation 'omega'.  W is *NOT* the frequency:
%   the frequency is W/(2*pi).
% If T is considered as a space coordinate,
%   W is a wave number (usually noted K = 2*PI/LAMBDA, where LAMBDA is a
%   wavelength).
%
%
% compare different ways to get the spectra
% 
% benefit of getting spectra from the spatials correlations is that 
% gives spectra in wave numbers, gets you E(k) in Pope as opposed to E(f)
% in freq space (what you get if you do fft of time series to get spectra)
% 
% to get E(k) take fft of the correction variable 'c' with space coordinate dr
% for E(k) start reading at Pope eqn 6.204
% compute the epsilon 6.191

if flag_plot
    
    % plot the time series
    figure
    hold on
    plot(1:numel(vx), vx, '-r');
    title('time series from PIV bin')
    xlabel('time, t [s]')
    ylabel('u_x [m/s]')
    box on; grid on;

    set(findall(gca,'-property','FontSize'),'FontSize',14)
    set(findall(gca,'-property','LineWidth'),'LineWidth',1)


    
    figure
    % subplot(1,4,1)
    subplot(1,3,1)
        hold on
        [Wf,Ef,hf] = ezfft(t,vx,'freq','hann','disp');
        set(gca, 'XScale', 'log', 'YScale', 'log')
        set(hf, 'LineStyle', '-', 'Color', 'r');
        title('"ezfft" method')
        box on; grid on;

        % plot with -5/3 line
        x       = [0.1*max(Wf), 0.9*max(Wf)];
        ynoll   = 0.001*max(Ef);
        y(1)    = ynoll;
        y(2)    = y(1)*(x(2)/x(1))^(-5/3);
        plot(x,y,'k-','linew',4)

    % subplot(1,4,2)
    subplot(1,3,2)
        hold on
        [Ws,Es,hs] = ezfft(t,vx,'space','hann','disp');
        set(gca, 'XScale', 'log', 'YScale', 'log')
        set(hs, 'LineStyle', '-', 'Color', 'k');
        title('"ezfft" method')
        box on; grid on;

        % plot with -5/3 line
        x       = [0.1*max(Ws), 0.9*max(Ws)];
        ynoll   = 0.001*max(Es);
        y(1)    = ynoll;
        y(2)    = y(1)*(x(2)/x(1))^(-5/3);
        plot(x,y,'k-','linew',4)

    % subplot(1,4,3)
    subplot(1,3,3)
        hold on
        win_pts = 96;
        overlap = 0.5;
        [Pxx, f, checksum] = speed_fft(vx, t, win_pts, overlap);
        plot(f, Pxx, 'r-o', 'linew', 2)
        set(gca, 'XScale', 'log', 'YScale', 'log')
        xlabel('frequency bands');
        ylabel('power spectral density (units of variance)')
        title('"speed-fft" method')
        box on; grid on;

        % plot with -5/3 line
        x       = [0.5*max(f), 0.9*max(f)];
        ynoll   = 0.9*max(Pxx);
        y(1)    = ynoll;
        y(2)    = y(1)*(x(2)/x(1))^(-5/3);
        plot(x,y,'k-','linew',4)

    % subplot(1,4,4)
    %     % If F is an array of fields, the average spectrum is returned
    % %     F1=loadVec(eval(PIV_locations{i}));    % assuming v in m/s and r in mm
    %     F1 = loadVec('/mnt/data-RAID-1/danny/backup_Rivendell/from_HAL3/Data/STC/Clipped_7D_DOWN_TSR_6.mat');
    %     for j=1:numel(F1)
    %       if isempty(F1(j).vx)==0;
    %           F(counter)=F1(j);
    %           counter=counter+1;
    %       else
    %           continue
    %       end
    %     end
    %     


%     %     sp = specf(truncf(F,2));
%     %     sp = specf(vx);
% %         sp = specf(truncf(vx));
% %         sp = specf(vx,'w');
%         ff = loadvec(
%         F.w = vx;
% %         sp = specf(f.vx,'w');
%         sp = specf(F.w,'w');
%         loglog(sp.k, sp.el,'LineWidth',2); 
%         set(gca,'FontSize',14,'LineWidth',0.75)
%         xlabel('k (m^{-1})'); 
%         ylabel('E(k)  (m^3 s^{-2})');
    
  
end




%% danny: version of the 'phuck' method
% use the entire time series and estimate where the signal becomes uncorrelated
% % L       = numel(vx);
% % Lt      = floor(numel(vx)./L);
% % Lmax    = L;
% % cvx     = zeros(1, 2*L-1);
% % % calc stats
% % vxm     = mean(vx);
% % % vx2     = mean(vx.^2);
% % vxstd   = std(vx);
% % 
% % % st = 1:Lt:numel(vx);
% % % define part of signal to use
% % vxt = vx(1:end-1);
% % % reduced sigmal, remove mean, unit variance
% % vxtr = (vxt-vxm)./vxstd;
% % c    = xcorr(vxtr,'none');
% % 
% % figure
% % plot([1:numel(cvx)].*1e-4.*vxm,c./c(1),'-','color',col(i,:))
% % plot([1:numel(cvx)] .* vxm * 1e-4, c./c(1),'r-')
% % 
% % Nvxtr = numel(vxtr);
% % jj=numel(vxt);
% % 
% % c = xcorr(vxtr,'none');
% % % cvx(Lmax-jj+1:Lmax+jj-1)=cvx(Lmax-jj+1:Lmax+jj-1)+c';
% % cvx(Lmax-Nvxtr+1:Lmax+Nvxtr-1) = cvx(Lmax-Nvxtr+1:Lmax+Nvxtr-1)+c';
% % 
% % 
% % % 
% % % % define part of signal to use
% % %         vxt = vx(me:me+L-1);
% % %         % reduced sigmal, remove mean, unit variance
% % %         vxtr=(vxt-vxm)./vxstd;
% % %         jj=numel(vxt);
% % %         
% % %         c=xcorr(vxtr,'none');
% % %         cvx(Lmax-jj+1:Lmax+jj-1)=cvx(Lmax-jj+1:Lmax+jj-1)+c';
% % %     
% % % % now step through the signal in increments, to get multiple estimates
% % % 
% % 
% % % use the entire time series and estimate where the signal becomes uncorrelated
% % L       = numel(vx);
% % Lt      = floor(numel(vx)./L);
% % Lmax    = L;
% % cvx     = zeros(1, 2*L-1);
% % % calc stats
% % vxm     = mean(vx);
% % % vx2     = mean(vx.^2);
% % vxstd   = std(vx);
% % 
% % % colormaps
% % % figure
% % % col = jet(numel(vx));
% % 
% % % me=0;
% % n  = 0;
% % st = 1:Lt:numel(vx);
% % for i =1:5:numel(st)-1;
% % 
% %     n=n+1;
% %     me=st(i);
% %     if me+L>numel(vx)
% %         continue
% %     end
% %     % define part of signal to use
% %     vxt=vx(me:me+L-1);
% %     % reduced sigmal, remove mean, unit variance
% %     vxtr=(vxt-vxm)./vxstd;
% %     jj=numel(vxt);
% % 
% %     c=xcorr(vxtr,'none');
% %     cvx(Lmax-jj+1:Lmax+jj-1)=cvx(Lmax-jj+1:Lmax+jj-1)+c';
% %         
% %     hold on
% %     plot([1:numel(cvx)].*1e-4.*vxm,c./c(1),'-','color',col(i,:))
% % end


%% note: by Taylor frozen turbulence, can convert time to length
% L=[800:200:8000];       % the length of a chunk of data where data is expected to be correlated
% L = 50:50:numel(vx);          % consider the entire data set correlated
% Ni = 10;
% Nt = 10;                         % the integration interval, time lag
% Li = 1:Ni:numel(vx);          % consider the entire data set correlated
% Li = 1:Ni:200;          % consider only first segment of the dataset correlated
% Li = 1:5:50;          % consider a specific segment of the dataset correlated
% Li = 50:50:numel(vx);          % consider the entire data set correlated

% col=jet(numel(Li));
% colc=0;

% L is the max length of the trajectory pieces we are using
if flag_plot
    hfig1 = figure;
    hfig2 = figure;
end

col1 = lines(numel(Li));    % pick colors for plots
% for L=[7200];
% for L=50;
% for L = L  
m = 1;
for L = Li 
%     L
%     Li

    Lt   = floor(numel(vx)./L);
    Lmax = L;
    cvx  = zeros(1,2*L-1);
    % calc stats
    vxm   = mean(vx);
%     vx2   = mean(vx.^2);
    vxstd = std(vx);
    
    
%     me=0;
    n    = 0;
%     Nstart = 50
%     Nend   = 150
%     st = Nstart:Lt:Nend;
    st   = 1:Lt:numel(vx);      % windows in the time signal
    col2 = lines(numel(st)-1);  % choose colors
%     Nt  = 5;
    for i =1:Nt:numel(st)-1;
        
        n  = n+1;
        me = st(i);
        if me+L>numel(vx)
            continue
        end
        % define part of signal to use
        vxt = vx(me:me+L-1);
        
        % reduced sigmal, remove mean, unit variance
        vxtr = (vxt-vxm)./vxstd;
        n    = numel(vxt);
        
        % Cross-correlation function estimates, of the reduced signal
        c = xcorr(vxtr,'none');
        cvx(Lmax-n+1:Lmax+n-1) = ...
        cvx(Lmax-n+1:Lmax+n-1) + c';
        
        if flag_plot
            figure(hfig1);
            hold on
%             plot([1:numel(cvx)].*1e-4.*vxm,c./c(1),'-','color',col(i,:))
%             plot((1:numel(cvx)).*vxm, c./c(1),'-','color',col2(i,:))
            plot(1:numel(cvx), c./c(1),'-','color',col2(i,:))
            xlabel('bin')
            ylabel('')
            title('cross-correlation of reduced signal')
            box on; grid on;
        end
    end
    % take average correlation
    cvxt = cvx./n;
    % take positive time correlation
%     cvxt = cvxt((numel(cvx)-1)/2+1:end);
    cvxtp = cvxt((numel(cvx)-1)/2+1:end);
    
    
    if flag_plot
        figure(hfig2)
        hold on;
%         plot([1:numel(cvxt)].*1e-4.*vxm,cvxt./cvxt(1),'-','color',col)
%         plot([1:numel(cvxt)].*vxm, cvxt./cvxt(1),'-','color','k')
%         plot((1:numel(cvxt)).*vxm, cvxt./cvxt(1),'-','color',col1(m,:))
%         plot(1:numel(cvxt), cvxt./cvxt(1),'-','color',col1(m,:))
        plot(1:numel(cvxtp), cvxtp./cvxtp(1),'-','color',col1(m,:))
        xlabel('bin')
        ylabel('correlation correction variable, c')
        title('correlation of reduced signal')

        box on; grid on;
        
        set(findall(gca,'-property','FontSize'),'FontSize',14)
        set(findall(gca,'-property','LineWidth'),'LineWidth',1)
                
    end
%     
    m = m+1;
end

%% integrate to calculate L_int
% dt=1/10000;
% dt=1/100;
dt = t(2) - t(1);

dr=dt*mean(vx);

% Lt=cumsum(cvxt./cvxt(1)).*dr;
% L_turb=cumsum(cvxt./cvxt(1)).*dr;
L_turb=cumsum(cvxtp./cvxtp(1)).*dr;

% refer to as integral length scale

% get the Taylor microscale by fitting a parabola or like 1/ ()see p 198)
% L_11 compute the paraboloa thing p(r) like eqn 6.54
% microscale is inbetween L and the Kolmnogorov length scale
if flag_plot
    figure()
    hold on;
%     plot([1:numel(cvxt)].*1e-4.*vxm,L_turb,'b-')
    plot(1:numel(cvxt), L_turb,'b-')
    xlabel('bin')
    ylabel('integral length scale, L_{turb}')
    title('integration of correlation function')
    box on; grid on;
    
    set(findall(gca,'-property','FontSize'),'FontSize',14)
    set(findall(gca,'-property','LineWidth'),'LineWidth',1)
end




end