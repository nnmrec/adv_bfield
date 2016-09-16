%Brian Polagye
%July 13, 2013

%Description: compare turbulence spectra

clear

%% Setup

compare_cases = [16,17];    %data files to compare

batch_file = 'FlumeData.xlsx';

%% Load batch data

%load batch data
[~,~,batch_list] = xlsread(batch_file);

file_paths = batch_list(2:end,2);           %vector (converted) file paths
file_names = batch_list(2:end,3);           %vector file names

%% Overlay spectra


for i = compare_cases
    load([file_paths{i} file_names{i} '.mat'])
    
    if i==compare_cases(1)
        figure(1)
        clf
        
        loglog(spectra.f,spectra.Pxx,'-b','linewidth',1)
        hold on
        loglog(spectra.f,spectra.Pyy,'-r','linewidth',1)
        loglog(spectra.f,spectra.Pzz,'-g','linewidth',1)
    else
                
        loglog(spectra.f,spectra.Pxx,'--b','linewidth',1)
        hold on
        loglog(spectra.f,spectra.Pyy,'--r','linewidth',1)
        loglog(spectra.f,spectra.Pzz,'--g','linewidth',1)
        
    end
    
end