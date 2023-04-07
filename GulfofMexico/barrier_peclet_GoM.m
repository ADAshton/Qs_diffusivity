% calculate Peclet

load('GoM_Barrier_Mapped_Data.mat')
load('sl_az_v2.mat')
BarrierIsland_name = strrep(BarrierIsland, '_', ' ');
barrier_num = [1:length(BarrierIsland);1:length(BarrierIsland)]';
MeanMinDepthm(11:12) = 1.4;
BarrierIsland_name{11} = 'Scofield';
BarrierIsland_name{12} = 'Pelican';
BarrierIsland_name{15} = 'Timbalier West';
BarrierIsland_name{16} = 'Timbalier East';
BarrierIsland_name{10} = 'Isle Grande Terre';
BarrierIsland_name{6} = 'Whiskey';
BarrierIsland_name{17} = 'Trinity';

ind_West_East = [2 6 17 15 16 7 8 10 1 11 12 4 3 14 9 13 5];

% load('diffusivity_gulfbarriers_v2.mat') % diffusivity is output as -HTscale * psi (m^2/s)
load('diffusivity_040723.mat')
ind_neg = logical(double((diffusivity_m2s<0))+double(isnan(diffusivity_m2s)));
diffusivity_m2y = diffusivity_m2s*365*24*60*60; % (dsf = 10; K2 = 0.17; already included in net diffusivity) m^2/s-->m^2/yr
diffusivity_kmyr = diffusivity_m2y/1e6;
diffusivity = diffusivity_m2y;
unusuble_barriers = BarrierIsland_name(ind_neg);
usable_barriers_ind = ~ind_neg;

% calc overwash flux
u_ow_max = OverwashfluxMaxm3my./(MeanElevationm+MeanMinDepthm); % max range
u_ow_min = Overwashfluxminm3my./(MeanElevationm+MeanMinDepthm); % min range
u_ow_meas = OverwashFluxm3my./(MeanElevationm+MeanMinDepthm); % measured/calculated from literature
u_ow_mean = mean([OverwashfluxMaxm3my Overwashfluxminm3my],2)./(MeanElevationm+MeanMinDepthm); % min range
% Alongshore length scale
% L_ast = AlongshoreLengthm;
L_ast = AlongshoreSeagullLengthm;
L_ow = MeanWidthm;
L_ast_c = sqrt((abs(diffusivity).*L_ow)./(u_ow_mean));
% calculate peclet
Pe_barrier_max = (u_ow_max.*(L_ast.^2))./(L_ow.*(diffusivity));
Pe_barrier_min = (u_ow_min.*(L_ast.^2))./(L_ow.*(diffusivity));
Pe_barrier_mean = (u_ow_mean.*(L_ast.^2))./(L_ow.*(diffusivity));
Pe_barrier_meas = (u_ow_meas.*(L_ast.^2))./(L_ow.*(diffusivity));

%% plot
% peclet range
figure()
hold on
% plot range of peclet given range of ow flud
for i = 1:length(barrier_num)
    plot(barrier_num(i,:),[Pe_barrier_min(ind_West_East(i)) Pe_barrier_max(ind_West_East(i))],'k','LineWidth',2)
end
% plot measured peclet given ow flux
scatter(barrier_num(:,1),Pe_barrier_meas(ind_West_East),50,'ro')
ylim([0 max([Pe_barrier_max;Pe_barrier_min;Pe_barrier_meas])])
set(gca,'yscale','log')
yline(1,'LineWidth',2);
set(gca,'xtick',[1:length(barrier_num)],'xticklabel',BarrierIsland_name(ind_West_East))
xtickangle(45)
ylabel('Barrier Peclet')
set(gca,'FontSize',16)
%% peclet order least to most stable
ind_all_notneg = find(~ind_neg);
Pe_barrier_max_notneg = Pe_barrier_max(ind_all_notneg);
Pe_barrier_min_notneg = Pe_barrier_min(ind_all_notneg);
Pe_barrier_mean_notneg = Pe_barrier_mean(ind_all_notneg);
Pe_barrier_meas_notneg = Pe_barrier_meas(ind_all_notneg);
barrier_num_notneg = barrier_num(ind_all_notneg,:);
barrier_name_notneg = BarrierIsland_name(ind_all_notneg,:);

[~,ind_sort] = sort(Pe_barrier_max_notneg,'descend');

figure()
hold on
% plot range of peclet given range of ow flud
for i = 1:length(ind_sort)
    plot([i i],[Pe_barrier_min_notneg(ind_sort(i)) Pe_barrier_max_notneg(ind_sort(i))],'k','LineWidth',2)
    scatter(i,Pe_barrier_mean_notneg(ind_sort(i)),1000,'.','k')
    scatter(i,Pe_barrier_meas_notneg(ind_sort(i)),50,'ro')
end
% plot measured peclet given ow flux
ylim([0 max([Pe_barrier_max;Pe_barrier_min;Pe_barrier_meas_notneg])])
set(gca,'yscale','log')
yline(1,'LineWidth',2);
set(gca,'xtick',[1:length(barrier_num_notneg)],'xticklabel',barrier_name_notneg(ind_sort))
xtickangle(45)
ylabel('Barrier Peclet')
set(gca,'FontSize',16)
%% stable lengths
figure()
hold on
for i = 1:length(ind_sort)
    scatter(i,L_ast_c(ind_sort(i))/1000,1000,'.','r')
    scatter(i,L_ast(ind_sort(i))/1000,1000,'.','k')
end
ylabel('Alongshore length (km)')
set(gca,'xtick',[1:length(barrier_num_notneg)],'xticklabel',barrier_name_notneg(ind_sort))
xtickangle(45)
set(gca,'FontSize',16)
legend('L_a_s_t_,_c','L_a_s_t')

figure()
hold on
for i = 1:length(ind_sort)
    scatter(i,(L_ast(ind_sort(i))/1000)./(L_ast_c(ind_sort(i))/1000),1000,'.','k')
%     scatter(i,L_ast(ind_sort(i))/1000,1000,'.','k')
end
yline(1)
ylabel('L_a_s_t / L_a_s_t_,_c (km)')
set(gca,'xtick',[1:length(barrier_num_notneg)],'xticklabel',barrier_name_notneg(ind_sort))
xtickangle(45)
set(gca,'FontSize',16)
% legend('L_a_s_t_,_c','L_a_s_t')
%% peclet vs
figure();
% % sinuosity
% subplot(3,2,1)
% scatter(AlongshoreSinuosity(usable_barriers_ind),Pe_barrier_min(usable_barriers_ind))
% set(gca,'yscale','log')
% yline(1);
% ylabel('Barrier Peclet')
% xlabel('Sinuosity')
% set(gca,'FontSize',16)
% elevation
subplot(3,2,3)
elev = MeanElevationm(usable_barriers_ind);
[elev_sort,ind_elev]= sort(elev);
Pe_barrier_mean_scatter = Pe_barrier_mean(usable_barriers_ind);
%next: make a matrix with a column of ones and a column of the x values
X=[ones(length(elev_sort),1) elev_sort];    %X=(column of ones,column of x values)
Y=log10(Pe_barrier_mean_scatter(ind_elev));        %log(height), for the regression
%next line: the standard matrix equation for a linear regression
b=inv(X'*X)*X'*Y;       %b(1)=intercept, b(2)=slope
predPe=10.^(X*b);   %the heights predicted by the linear regresison
[~,~,~,~,stats] = regress(Y,X);
r_sq = stats(1);
semilogy(elev_sort,predPe,'k-');
hold on
scatter(elev_sort,Pe_barrier_mean_scatter(ind_elev),100,'.','k')
set(gca,'yscale','log')
yline(1);
ylabel('Barrier Peclet')
xlabel('Elevation (m)')
set(gca,'FontSize',16)
title(r_sq)
% width
subplot(3,2,2)
width = MeanWidthm(usable_barriers_ind);
[width_sort,ind_width]= sort(width);
%next: make a matrix with a column of ones and a column of the x values
X=[ones(length(width_sort),1) width_sort];    %X=(column of ones,column of x values)
Y=log10(Pe_barrier_mean_scatter(ind_width));        %log(height), for the regression
%next line: the standard matrix equation for a linear regression
b=inv(X'*X)*X'*Y;       %b(1)=intercept, b(2)=slope
predPe=10.^(X*b);   %the heights predicted by the linear regresison
[~,~,~,~,stats] = regress(Y,X);
r_sq = stats(1);
semilogy(width_sort,predPe,'k-');
hold on
scatter(MeanWidthm(usable_barriers_ind),Pe_barrier_mean(usable_barriers_ind),100,'.','k')
set(gca,'yscale','log')
yline(1);
ylabel('Barrier Peclet')
xlabel('Width (m)')
set(gca,'FontSize',16)
title(r_sq)
% depth
subplot(3,2,4)
depth = MeanMinDepthm(usable_barriers_ind);
[depth_sort,ind_depth]= sort(depth);
%next: make a matrix with a column of ones and a column of the x values
X=[ones(length(depth_sort),1) depth_sort];    %X=(column of ones,column of x values)
Y=log10(Pe_barrier_mean_scatter(ind_depth));        %log(height), for the regression
%next line: the standard matrix equation for a linear regression
b=inv(X'*X)*X'*Y;       %b(1)=intercept, b(2)=slope
predPe=10.^(X*b);   %the heights predicted by the linear regresison
[~,~,~,~,stats] = regress(Y,X);
r_sq = stats(1);
semilogy(depth_sort,predPe,'k-');
hold on
scatter(MeanMinDepthm(usable_barriers_ind),Pe_barrier_min(usable_barriers_ind),100,'.','k')
set(gca,'yscale','log')
yline(1);
ylabel('Barrier Peclet')
xlabel('Depth (m)')
set(gca,'FontSize',16)
title(r_sq)
% He
subplot(3,2,5)
he = MeanMinDepthm(usable_barriers_ind)+MeanElevationm(usable_barriers_ind);
[he_sort,ind_he]= sort(he);
%next: make a matrix with a column of ones and a column of the x values
X=[ones(length(he_sort),1) he_sort];    %X=(column of ones,column of x values)
Y=log10(Pe_barrier_mean_scatter(ind_he));        %log(height), for the regression
%next line: the standard matrix equation for a linear regression
b=inv(X'*X)*X'*Y;       %b(1)=intercept, b(2)=slope
predPe=10.^(X*b);   %the heights predicted by the linear regresison
[~,~,~,~,stats] = regress(Y,X);
r_sq = stats(1);
semilogy(he_sort,predPe,'k-');
hold on
scatter(MeanMinDepthm(usable_barriers_ind)+MeanElevationm(usable_barriers_ind),Pe_barrier_min(usable_barriers_ind),100,'.','k')
set(gca,'yscale','log')
yline(1);
ylabel('Barrier Peclet')
xlabel('Effective Height (m)')
set(gca,'FontSize',16)
title(r_sq)
% alongshore length
subplot(3,2,1)
Al = AlongshoreLengthm(usable_barriers_ind);
[Al_sort,ind_Al]= sort(Al);
%next: make a matrix with a column of ones and a column of the x values
X=[ones(length(Al_sort),1) Al_sort];    %X=(column of ones,column of x values)
Y=log10(Pe_barrier_mean_scatter(ind_Al));        %log(height), for the regression
%next line: the standard matrix equation for a linear regression
b=inv(X'*X)*X'*Y;       %b(1)=intercept, b(2)=slope
predPe=10.^(X*b);   %the heights predicted by the linear regresison
[~,~,~,~,stats] = regress(Y,X);
r_sq = stats(1);
semilogy(Al_sort/1000,predPe,'k-');
hold on
scatter(AlongshoreLengthm(usable_barriers_ind)/1000,Pe_barrier_min(usable_barriers_ind),100,'.','k')
set(gca,'yscale','log')
yline(1);
ylabel('Barrier Peclet')
xlabel('Alongshore Length (km)')
set(gca,'FontSize',16)
title(r_sq)
% OW flux
subplot(3,2,6)
Qow = OverwashFluxm3my(usable_barriers_ind);
[Qow_sort,ind_Qow]= sort(Qow);
ind_Qow = ind_Qow(1:3);
%next: make a matrix with a column of ones and a column of the x values
X=[ones(length(ind_Qow),1) Qow_sort(1:3)];    %X=(column of ones,column of x values)
Y=log10(Pe_barrier_mean_scatter(ind_Qow));        %log(height), for the regression
%next line: the standard matrix equation for a linear regression
b=inv(X'*X)*X'*Y;       %b(1)=intercept, b(2)=slope
predPe=10.^(X*b);   %the heights predicted by the linear regresison
[~,~,~,~,stats] = regress(Pe_barrier_mean_scatter(ind_Qow),X);
r_sq = stats(1);
pp = polyfit(Qow_sort(1:3),Pe_barrier_mean_scatter(ind_Qow),1);
plot(Qow_sort(1:3),polyval(pp,Qow_sort(1:3)),'k-');
hold on
scatter(Qow_sort(1:3),Pe_barrier_mean_scatter(ind_Qow),100,'.','k')
% set(gca,'yscale','log')
yline(1);
ylabel('Barrier Peclet')
xlabel('Overwash flux (m^3/m/yr)')
set(gca,'FontSize',16)
title(r_sq)
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

% %% depth vs. plots
% % elevation - depth
% figure();
% subplot(2,2,1)
% scatter(MeanMinDepthm(usable_barriers_ind),MeanElevationm(usable_barriers_ind))
% ylabel('Elevation')
% xlabel('Depth')
% set(gca,'FontSize',16)
% % depth - width
% % figure();
% subplot(2,2,2)
% scatter(MeanMinDepthm(usable_barriers_ind),MeanWidthm(usable_barriers_ind))
% ylabel('Width')
% xlabel('Depth')
% set(gca,'FontSize',16)
% % depth - alongshore dist
% subplot(2,2,3)
% scatter(MeanMinDepthm(usable_barriers_ind),AlongshoreLengthm(usable_barriers_ind))
% ylabel('Alongshore length')
% xlabel('Depth')
% set(gca,'FontSize',16)
% 
% % depth - sinuosity
% subplot(2,2,4)
% scatter(MeanMinDepthm(usable_barriers_ind),AlongshoreSinuosity(usable_barriers_ind))
% ylabel('Sinuosity')
% xlabel('Depth')
% set(gca,'FontSize',16)
% %% width vs. plots
% % elevation
% figure();
% subplot(2,2,1)
% scatter(MeanWidthm(usable_barriers_ind),MeanElevationm(usable_barriers_ind))
% ylabel('Elevation')
% xlabel('Width')
% set(gca,'FontSize',16)
% % depth
% % figure();
% subplot(2,2,2)
% scatter(MeanWidthm(usable_barriers_ind),MeanMinDepthm(usable_barriers_ind))
% ylabel('Depth')
% xlabel('Width')
% set(gca,'FontSize',16)
% % alongshore length
% subplot(2,2,3)
% scatter(MeanWidthm(usable_barriers_ind),AlongshoreLengthm(usable_barriers_ind))
% ylabel('Alongshore length')
% xlabel('Width')
% set(gca,'FontSize',16)
% 
% % sinuosity
% subplot(2,2,4)
% scatter(MeanWidthm(usable_barriers_ind),AlongshoreSinuosity(usable_barriers_ind))
% ylabel('Sinuosity')
% xlabel('Width')
% set(gca,'FontSize',16)
% 
% %% alongshore length vs. plots
% % elevation
% figure();
% subplot(2,2,1)
% scatter(AlongshoreLengthm(usable_barriers_ind),MeanElevationm(usable_barriers_ind))
% ylabel('Elevation')
% xlabel('Alongshore length')
% set(gca,'FontSize',16)
% 
% % width
% subplot(2,2,2)
% scatter(AlongshoreLengthm(usable_barriers_ind),MeanWidthm(usable_barriers_ind))
% ylabel('Width')
% xlabel('Alongshore length')
% set(gca,'FontSize',16)
% 
% % depth
% subplot(2,2,3)
% scatter(AlongshoreLengthm(usable_barriers_ind),MeanMinDepthm(usable_barriers_ind))
% ylabel('Depth')
% xlabel('Alongshore length')
% set(gca,'FontSize',16)
% 
% % sinuosity
% subplot(2,2,4)
% scatter(AlongshoreLengthm(usable_barriers_ind),AlongshoreSinuosity(usable_barriers_ind))
% ylabel('Sinuosity')
% xlabel('Alongshore length')
% set(gca,'FontSize',16)
% 
% %% alongshore
% figure()
% % length vs. sinuosity
% subplot(2,2,1)
% scatter(AlongshoreLengthm(usable_barriers_ind),AlongshoreSinuosity(usable_barriers_ind))
% ylabel('Sinuosity')
% xlabel('Alongshore length')
% set(gca,'FontSize',16)
% % length vs. diffusivity
% subplot(2,2,2)
% scatter(AlongshoreLengthm(usable_barriers_ind),diffusivity(usable_barriers_ind))
% ylabel('Diffusivity')
% xlabel('Alongshore length')
% set(gca,'FontSize',16)
% % sinuosity vs. diffusivity
% subplot(2,2,3)
% scatter(AlongshoreSinuosity(usable_barriers_ind),diffusivity(usable_barriers_ind))
% xlabel('Sinuosity')
% ylabel('Diffusivity')
% set(gca,'FontSize',16)