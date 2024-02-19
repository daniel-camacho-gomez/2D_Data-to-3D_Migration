clearvars; close all hidden
%...Load the data 
[file_sim, file_dir] = uigetfile();
full_file_sim        = fullfile(file_dir, file_sim);
[fpath, fname, fext] = fileparts(full_file_sim);
addpath(file_dir);
load(full_file_sim);

%% 3D trajectories
figure
plot3(0,0,0,'xr')
hold on
plot3(x_sim,y_sim,z_sim,'color',[0 0.4470 0.7410],'LineWidth',1.5)
hold on 
plot3(x_sim(end,:),y_sim(end,:),z_sim(end,:),'linestyle','none','marker','o',...
    'markersize',15,'markerfacecolor',[0 0.4470 0.7410],'color','k',...
     'LineWidth',1.5)
grid on
xlim([-300 300])
ylim([-300 300])
zlim([-300 300])

xlabel('$X\;(\mu m)$','Interpreter','Latex','FontSize',15)
ylabel('$Y\;(\mu m)$','Interpreter','Latex','FontSize',15)
zlabel('$Z\;(\mu m)$','Interpreter','Latex','FontSize',15)

xticks(-200:100:200);
yticks(-200:100:200);

set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',15);