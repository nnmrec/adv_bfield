%Brian Polagye
%July 13, 2013

%Description: routine to post-process ADV data from flume characterization
%studies

%Input - ASCII format Nortek Vector files (using Nortek file conversion
%utility)

%Output - Processed velocity files
%   - despiked x, y, z velocity
%   - pressure
%   - time stamps

clear

%% Setup

data_proc = [1:43];    %data files to process

batch_file = 'FlumeData.xlsx';

addpath('despiking_toolbox')

parameters.overlap = 0.5;

%% Load batch data

%load batch data
[~,~,batch_list] = xlsread(batch_file);

file_paths = batch_list(2:end,2);           %vector (converted) file paths
file_names = batch_list(2:end,3);           %vector file names
x_poss = cell2mat(batch_list(2:end,4));     %x positions
y_poss = cell2mat(batch_list(2:end,5));     %y positions
z_poss = cell2mat(batch_list(2:end,6));     %z positions
descriptions = batch_list(2:end,7);         %test descriptions

%% Process ADV data

for i = data_proc
    
    %wrap up input parameters into structure
    parameters.x_pos = x_poss(i);
    parameters.y_pos = y_poss(i);
    parameters.z_pos = z_poss(i);
    parameters.description = descriptions{i};
    parameters.file_name = file_names{i};
    parameters.file_path = file_paths{i};
    
    batch_flume(parameters)
    
end