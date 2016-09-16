function [v0,v1,sRate]=vectrinoCalcPoint(dat_folder, timeID)
% Extracts the velocities, performs data checks and despiking for vectrino
% data taken with the v3 setup.  Currently assumes 2 vectrinos with
% filenames vectrinoData1_#fileIDnum#.mat and vectrinoData2_#fileIDnum#.mat

% Returns the velocity components and the sampling rate




% load('HWconfig.mat')
load([dat_folder filesep 'HWconfig.mat'])
% HWconfig = load([dat_folder filesep 'HWconfig.mat'])


% Load in point file
% load(strcat('point_',num2str(timeID),'.mat'));
load([dat_folder filesep strcat('point_',num2str(timeID),'.mat')]);

% Calculate positions of Vectrinos
v1.pos = pointDat.pos + vectrinoOrigin;
v0.pos = v1.pos + outboardOffset;


% Filecheck
% if isempty(dir(strcat('vectrinoData1_',num2str(timeID),'.mat')))
if isempty(dir([dat_folder filesep strcat('vectrinoData1_',num2str(timeID),'.mat')]))    
    error('Vectrino file %s not found. Have you converted the .ntk file?',strcat('vectrinoData1_',num2str(timeID),'.mat'))
end

%% Load/process outboard vectrino
% load(strcat('vectrinoData1_',num2str(timeID),'.mat'));
load([dat_folder filesep strcat('vectrinoData1_',num2str(timeID),'.mat')]);
% Extract the sample rate
sRate = Config.sampleRate;

% [v0.U,v0.badInds] = vectrinoDataScrub(Data,Config,0);
[v0.U,v0.badInds] = vectrinoDataScrub(Data,Config,2);
v0.U(:,3) = mean(v0.U(:,3:4),2);


% v0.U(:,1:2)=mtimesx(v0.U(:,1:2),[vectrinoRot.rotM0]');  % I had error with auto-compile on *nix, so the following is equavalent Matlab command:
v0.U(:,1:2) = v0.U(:,1:2) * [vectrinoRot.rotM0]';
v0.U(:,4)   = [];

%% Load/process inboard vectrino
load([dat_folder filesep strcat('vectrinoData2_',num2str(timeID),'.mat')]);
% CHeck that sample rate is the same
if Config.sampleRate~=sRate
    error('The two vectrino sample rates are not the same!')
end

[v1.U,v1.badInds]=vectrinoDataScrub(Data,Config,0);
v1.U(:,3)=mean(v1.U(:,3:4),2);
% v1.U(:,1:2)=mtimesx(v1.U(:,1:2),[vectrinoRot.rotM1]');
v1.U(:,1:2) = v1.U(:,1:2) * [vectrinoRot.rotM1]';
v1.U(:,4)=[];
