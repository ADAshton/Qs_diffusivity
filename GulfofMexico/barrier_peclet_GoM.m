% calculate Peclet

load('GoM_Barrier_Mapped_Data.mat')
BarrierIsland_name = strrep(BarrierIsland, '_', ' ');
barrier_num = [1:length(BarrierIsland);1:length(BarrierIsland)]';
MeanMinDepthm(13:14) = 1.4;

ind_West_East = [2 6 17 15 16 7 8 10 1 11 12 4 3 14 9 13 5];

load('diffusivity_gulfbarriers.mat')
ind_neg = (diffusivity<0);
unusuble_barriers = BarrierIsland_name(ind_neg);
usable_barriers_ind = ~ind_neg;

% calc overwash flux
u_ow_max = OverwashfluxMaxm3my./(MeanElevationm+MeanMinDepthm); % max range
u_ow_min = Overwashfluxminm3my./(MeanElevationm+MeanMinDepthm); % min range
u_ow_meas = OverwashFluxm3my./(MeanElevationm+MeanMinDepthm); % measured/calculated from literature

% Alongshore length scale
L_ast = AlongshoreLengthm;
L_ast = AlongshoreSeagullLengthm;
L_ow = MeanWidthm;

% calculate peclet
Pe_barrier_max = u_ow_max.*L_ast.^2./L_ow./(diffusivity);
Pe_barrier_min = u_ow_min.*L_ast.^2./L_ow./(diffusivity);
Pe_barrier_meas = u_ow_meas.*L_ast.^2./L_ow./(diffusivity);

%% plot
% peclet range
figure()
hold on
% plot range of peclet given range of ow flud
for i = 1:length(barrier_num)
plot(barrier_num(i,:),[Pe_barrier_min(ind_West_East(i)) Pe_barrier_max(ind_West_East(i))],'k','LineWidth',2)
end
% plot measured peclet given ow flux
scatter(ind_West_East,Pe_barrier_meas,50,'ro')
ylim([0 max([Pe_barrier_max;Pe_barrier_min;Pe_barrier_meas])])
set(gca,'yscale','log')
yline(1,'LineWidth',2);
set(gca,'xtick',[1:length(barrier_num)],'xticklabel',BarrierIsland_name(ind_West_East))
xtickangle(45)
ylabel('Barrier Peclet')
set(gca,'FontSize',16)
%% peclet vs
figure(); 
% sinuosity
subplot(3,2,1)
scatter(AlongshoreSinuosity(usable_barriers_ind),Pe_barrier_min(usable_barriers_ind))
set(gca,'yscale','log')
yline(1);
ylabel('Barrier Peclet')
xlabel('Sinuosity')
set(gca,'FontSize',16)
% elevation
subplot(3,2,2)
scatter(MeanElevationm(usable_barriers_ind),Pe_barrier_min(usable_barriers_ind))
set(gca,'yscale','log')
yline(1);
ylabel('Barrier Peclet')
xlabel('Elevation')
set(gca,'FontSize',16)
% width
subplot(3,2,3)
scatter(MeanWidthm(usable_barriers_ind),Pe_barrier_min(usable_barriers_ind))
set(gca,'yscale','log')
yline(1);
ylabel('Barrier Peclet')
xlabel('Width')
set(gca,'FontSize',16)
% depth
subplot(3,2,4)
scatter(MeanMinDepthm(usable_barriers_ind),Pe_barrier_min(usable_barriers_ind))
set(gca,'yscale','log')
yline(1);
ylabel('Barrier Peclet')
xlabel('Depth')
set(gca,'FontSize',16)
% alongshore length
subplot(3,2,5)
scatter(AlongshoreLengthm(usable_barriers_ind),Pe_barrier_min(usable_barriers_ind))
set(gca,'yscale','log')
yline(1);
ylabel('Barrier Peclet')
xlabel('Alongshore length')
set(gca,'FontSize',16)

%% aspect ratio (not useful because the relationship is from the equation)
% figure()
% subplot(2,2,1)
% scatter3(AlongshoreLengthm(usable_barriers_ind)./MeanWidthm(usable_barriers_ind),Pe_barrier_min(usable_barriers_ind),[1:length(find(usable_barriers_ind))],[],[1:length(find(usable_barriers_ind))])
% set(gca,'yscale','log')
% view(2)
% yline(1);
% ylabel('Barrier Peclet')
% xlabel('Alongshore length / Width (aspect ratio)')
% set(gca,'FontSize',16)
% 
% subplot(2,2,2)
% scatter3(AlongshoreLengthm(usable_barriers_ind)./MeanWidthm(usable_barriers_ind),AlongshoreSinuosity(usable_barriers_ind),[1:length(find(usable_barriers_ind))],[],[1:length(find(usable_barriers_ind))])
% set(gca,'yscale','log')
% view(2)
% yline(1);
% ylabel('Sinuosity')
% xlabel('Alongshore length / Width (aspect ratio)')
% set(gca,'FontSize',16)
% 
% subplot(2,2,3)
% scatter3(AlongshoreLengthm(usable_barriers_ind)./MeanWidthm(usable_barriers_ind),MeanMinDepthm(usable_barriers_ind),[1:length(find(usable_barriers_ind))],[],[1:length(find(usable_barriers_ind))])
% set(gca,'yscale','log')
% view(2)
% yline(1);
% ylabel('Depth')
% xlabel('Alongshore length / Width (aspect ratio)')
% set(gca,'FontSize',16)

%% depth vs. plots 
% elevation - depth
figure(); 
subplot(2,2,1)
scatter(MeanMinDepthm(usable_barriers_ind),MeanElevationm(usable_barriers_ind))
ylabel('Elevation')
xlabel('Depth')
set(gca,'FontSize',16)
% depth - width
% figure(); 
subplot(2,2,2)
scatter(MeanMinDepthm(usable_barriers_ind),MeanWidthm(usable_barriers_ind))
ylabel('Width')
xlabel('Depth')
set(gca,'FontSize',16)
% depth - alongshore dist
subplot(2,2,3)
scatter(MeanMinDepthm(usable_barriers_ind),AlongshoreLengthm(usable_barriers_ind))
ylabel('Alongshore length')
xlabel('Depth')
set(gca,'FontSize',16)

% depth - sinuosity
subplot(2,2,4)
scatter(MeanMinDepthm(usable_barriers_ind),AlongshoreSinuosity(usable_barriers_ind))
ylabel('Sinuosity')
xlabel('Depth')
set(gca,'FontSize',16)
%% width vs. plots 
% elevation
figure(); 
subplot(2,2,1)
scatter(MeanWidthm(usable_barriers_ind),MeanElevationm(usable_barriers_ind))
ylabel('Elevation')
xlabel('Width')
set(gca,'FontSize',16)
% depth
% figure(); 
subplot(2,2,2)
scatter(MeanWidthm(usable_barriers_ind),MeanMinDepthm(usable_barriers_ind))
ylabel('Depth')
xlabel('Width')
set(gca,'FontSize',16)
% alongshore length
subplot(2,2,3)
scatter(MeanWidthm(usable_barriers_ind),AlongshoreLengthm(usable_barriers_ind))
ylabel('Alongshore length')
xlabel('Width')
set(gca,'FontSize',16)

% sinuosity
subplot(2,2,4)
scatter(MeanWidthm(usable_barriers_ind),AlongshoreSinuosity(usable_barriers_ind))
ylabel('Sinuosity')
xlabel('Width')
set(gca,'FontSize',16)

%% alongshore length vs. plots 
% elevation
figure(); 
subplot(2,2,1)
scatter(AlongshoreLengthm(usable_barriers_ind),MeanElevationm(usable_barriers_ind))
ylabel('Elevation')
xlabel('Alongshore length')
set(gca,'FontSize',16)

% width
subplot(2,2,2)
scatter(AlongshoreLengthm(usable_barriers_ind),MeanWidthm(usable_barriers_ind))
ylabel('Width')
xlabel('Alongshore length')
set(gca,'FontSize',16)

% depth
subplot(2,2,3)
scatter(AlongshoreLengthm(usable_barriers_ind),MeanMinDepthm(usable_barriers_ind))
ylabel('Depth')
xlabel('Alongshore length')
set(gca,'FontSize',16)

% sinuosity
subplot(2,2,4)
scatter(AlongshoreLengthm(usable_barriers_ind),AlongshoreSinuosity(usable_barriers_ind))
ylabel('Sinuosity')
xlabel('Alongshore length')
set(gca,'FontSize',16)

%% alongshore
figure()
% length vs. sinuosity
subplot(2,2,1)
scatter(AlongshoreLengthm(usable_barriers_ind),AlongshoreSinuosity(usable_barriers_ind))
ylabel('Sinuosity')
xlabel('Alongshore length')
set(gca,'FontSize',16)
% length vs. diffusivity
subplot(2,2,2)
scatter(AlongshoreLengthm(usable_barriers_ind),diffusivity(usable_barriers_ind))
ylabel('Diffusivity')
xlabel('Alongshore length')
set(gca,'FontSize',16)
% sinuosity vs. diffusivity
subplot(2,2,3)
scatter(AlongshoreSinuosity(usable_barriers_ind),diffusivity(usable_barriers_ind))
xlabel('Sinuosity')
ylabel('Diffusivity')
set(gca,'FontSize',16)