%Brian Polagye
%July 13, 2013

%Description: Plot profiles along flume

clear

clear

plot_case = [1:43];    %data points to include in plot

batch_file = 'FlumeData.xlsx';

x_pos = [2, 4, 6, 8, 10];
y_pos = [-0.25, 0, 0.25];

plot_parameter = 1;
    % 1 = mean velocity
    % 2 = turbulence intensity

%% Load batch data

%load batch data
[~,~,batch_list] = xlsread(batch_file);

file_paths = batch_list(2:end,2);           %vector (converted) file paths
file_names = batch_list(2:end,3);           %vector file names

%% Load and store data

cnt = 1;

plot_val = zeros(length(plot_case),1);
x = zeros(size(plot_val));
y = zeros(size(plot_val));
z = zeros(size(plot_val));

for i = plot_case
   
    load([file_paths{i} file_names{i} '.mat'])
    
    switch plot_parameter
        case 1
            plot_val(cnt) = data.smean;
        case 2
            plot_val(cnt) = data.I;

    end
    
    x(cnt) = config.x_pos;
    y(cnt) = config.y_pos;
    z(cnt) = config.z_pos;
    
    cnt = cnt + 1;
    
end

%% Plot Results

figure(1)
clf

for i = 1:length(y_pos)
    for j = 1:length(x_pos)
        
        subplot(length(y_pos),length(x_pos),j+(i-1)*length(x_pos))
        
        pts = find(x==x_pos(j) & y==y_pos(i));
        temp = [plot_val(pts), z(pts)];
        temp = sortrows(temp,2);
        
        plot(temp(:,1),temp(:,2),'o--b','markerface','b')
        
        switch plot_parameter
            case 1
                set(gca,'XLim',[0.9 1.1])
                text(0.92,0.74,['x = ' num2str(x_pos(j)) ' m, y = ' num2str(y_pos(i)) ' m'],'fontweight','b')
                xlabel('U (m/s)','fontweight','b')
            case 2
                set(gca,'XLim',[0 0.1])
                text(0.01,0.74,['x = ' num2str(x_pos(j)) ' m, y = ' num2str(y_pos(i)) ' m'],'fontweight','b')
                xlabel('I','fontweight','b')
        end
        set(gca,'YLim',[0 0.8])
        
        grid on
        
        ylabel('z (m)','fontweight','b')
        
    end
end
