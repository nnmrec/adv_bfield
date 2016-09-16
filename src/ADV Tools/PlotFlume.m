%Brian Polagye
%July 13, 2013

%Description: color plot of flume chracteristics

clear

plot_case = [1:43];    %data poiunts to include in plot

batch_file = 'FlumeData.xlsx';

plot_parameter = 2;
    % 1 = mean velocity
    % 2 = turbulence intensity
    % 3 = % spikes

%% Load batch data

%load batch data
[~,~,batch_list] = xlsread(batch_file);

file_paths = batch_list(2:end,2);           %vector (converted) file paths
file_names = batch_list(2:end,3);           %vector file names

%% Load and store data

cnt = 1;

for i = plot_case
   
    load([file_paths{i} file_names{i} '.mat'])
    
    switch plot_parameter
        case 1
            plot_val(cnt) = data.smean;
        case 2
            plot_val(cnt) = data.I;
        case 3
            plot_val(cnt) = anc.num_spikes;
    end
    
    x(cnt) = config.x_pos;
    y(cnt) = config.y_pos;
    z(cnt) = config.z_pos;
    
    cnt = cnt + 1;
    
end

%% Plot results

figure(1)
clf

set(gcf,'position',[84 202 1398 596])

colormap(lbmap(100,'Blue'))

scatter3(x,y,z,100,plot_val,'fill')

xlabel('Distance from inlet (m)')
ylabel('Laterial distance (m)')
zlabel('Distance from bottom (m)')

grid on
axis equal

ca = colorbar;
axes(ca)

switch plot_parameter
    case 1
        ylabel('Velocity (m/s)')
    case 2
        ylabel('Turbulence Intensity')
    case 3
        ylabel('Fraction spikes')
    
end