%Brian Polagye
%July 13, 2013

%Description: report flume test statistics

clear

%% Setup

data_proc = [38:43];    %data files to process

batch_file = 'FlumeData.xlsx';

%% Load batch data

%load batch data
[~,~,batch_list] = xlsread(batch_file);

file_paths = batch_list(2:end,2);           %vector (converted) file paths
file_names = batch_list(2:end,3);           %vector file names

%% Generate stats

for i = data_proc
    load([file_paths{i} file_names{i} '.mat'])
    
    out_table(i,1) = anc.num_spikes;
    out_table(i,2) = data.smean;
    out_table(i,3) = data.I;
    
end