%%  Process ADV Tow Test
%   R. Cavagnaro, 8/31/12
%   This script takes in a .hdr and .dat file and processes velocity and
%   time for an ADV output.

clear all
close all

operation = 1; % 0 is to export velocity data, 1 is to output sen data

outPath = 'C:\Users\Rob\Desktop\Micropower\';

if operation == 0
    outFile = 'FwdADV_10_08_12';
    fileName{1} = 'FwdTurbineMount_VEC_08Oct2012.dat';  % .dat output by ADV
    fileName{2} = 'FwdTurbineMount_VEC_08Oct2012.hdr';  % .hdr output by ADV
    rawData = importdata(fileName{1});
    rawHdr = importdata(fileName{2},'\t');
    fileStart = rawHdr{7,1};
    fileStart = strrep(fileStart,'Time of first measurement             ','');
    fileStart = datenum(fileStart);
    sampRate = rawHdr{12,1};
    sampRate = strrep(sampRate,'Sampling rate                         ','');
    sampRate = strrep(sampRate,' Hz','');
    sampRate = str2num(sampRate);
    data.burstCounter = rawData(:,1);
    data.ensembleCounter = rawData(:,2);
    data.velocityEast = rawData(:,3);
    data.velocityNorth = rawData(:,4);
    data.velocityUp = rawData(:,5);
    data.ampBeam1 = rawData(:,6);
    data.ampBeam2 = rawData(:,7);
    data.ampBeam3 = rawData(:,8);
    data.SNR1 = rawData(:,9);
    data.SNR2 = rawData(:,10);
    data.SNR3 = rawData(:,11);
    data.corr1 = rawData(:,12);
    data.corr2 = rawData(:,13);
    data.corr3 = rawData(:,14);
    data.pressure = rawData(:,15);
    data.checksum = rawData(:,18);
    data.magVel = (data.velocityEast.^2+data.velocityNorth.^2).^0.5;
    data.dirVel = 180/pi*atan2(data.velocityNorth, data.velocityEast);
    clear rawData

    t = 1/sampRate;
    nSamp = length(data.burstCounter);
    time = zeros(nSamp,1)+fileStart;
    dt = t:t:(nSamp*t);
    dt = dt(:)./86400;
    time = time + dt;

    save([outPath outFile], 'data', 'time');
else
    outFile = 'FwdADV_Sen_10_08_12';
    fileName{1} = 'FwdTurbineMount_VEC_08Oct2012.sen';  % .sen output by ADV
    rawSen = importdata(fileName{1});
    time = datenum(rawSen(:,3),rawSen(:,1),rawSen(:,2),rawSen(:,4),rawSen(:,5),rawSen(:,6));
    data.soundSpeed = rawSen(:,10);
    data.heading = rawSen(:,11);
    data.pitch = rawSen(:,12);
    data.roll = rawSen(:,13);
    data.temperature = rawSen(:,14);
    
    save([outPath outFile], 'data', 'time');
end
