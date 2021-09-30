function calculateQs_Dif(inputFolder_QsDif,outputFolder_QsDif,SLdata_detrend,ind_open,RotationAngle,plot_on)

% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isdir(inputFolder_QsDif)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', inputFolder_QsDif);
  uiwait(warndlg(errorMessage));
  return;
end
if ~isdir(outputFolder_QsDif)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', outputFolder_QsDif);
  uiwait(warndlg(errorMessage));
  return;
end

% Get a list of all files in the folder with the desired file name pattern.
infilePattern = fullfile(inputFolder_QsDif, '*.mat'); % Change to whatever pattern you need.
input_files = dir(infilePattern);

Xrot = SLdata_detrend(:,1);
Yrot = SLdata_detrend(:,2);

for runs = 1 : length(input_files)
    tic;
    
    baseFileName = input_files(runs).name;
    fullFileName = fullfile(inputFolder_QsDif, baseFileName);
    %fprintf(1, 'Now reading %s\n', fullFileName);

    [pathstr,filename,ext] = fileparts(fullFileName);
        load(fullFileName);  
    savename = [outputFolder_QsDif filename '_qsdif.mat'];
    
    avewidth = 7; %cells per side - it "smooths" the result
    pi = 3.1415;
    degtorad = pi/180;
    
    SL_open = SLdata_detrend(ind_open,:);
    Xrot = SLdata_detrend(:,1);
    Yrot = SLdata_detrend(:,2);
    
    openl = length(SL_open);
    
    SLAng = zeros(openl,1);
    SLAve = zeros(openl,2);
    
    xliml = min(Xrot);
    xlimr = max(Xrot);
    
    ylimt = max(Yrot);
    ylimb = min(Yrot);
    
    crossmin = 0;
    crossmax = ylimt-ylimb;
    
    Xrotopen = Xrot(ind_open);
    Yrotopen = Yrot(ind_open);
    
    
    % wave data
    
    angles = ShAngs;
    Gamma = Gammas(:,1);
    Qs = Qss(:,1);
    Qsglobal = QsRaws(:,1);
    Diff = TrueDiffs(:,1);
    Diffglobal = TrueDiffs(:,1);
    
    %plotting
%     xlimr = Limxr;
%     xliml = Limxl;
    
    
    % smooth shoreline
    for i = 1+avewidth:openl-avewidth
        ii  = (i-avewidth): (i + avewidth);
        p = polyfit(Xrotopen(ii),Yrotopen(ii),1);
        xi = ii(floor(length(ii)/2));
        y = polyval(p,Xrotopen(xi));
        SLAve(i,:) = [Xrotopen(xi) y];
        SLAng(i) = p(1);
    end
    
    SLAng =  - atan(SLAng)/degtorad - RotationAngle;
    
    SLAng = mod(SLAng,360);
    
    spacing = 1;
    ii = 1:spacing:length(SLAve);
    SL = SLAve(ii,:);
    Ang = SLAng(ii,:);
    
    AngRot(1:length(Ang)) = 0;
    GammaSL = interp1(angles,Gamma,Ang);
    QsSL = interp1(angles,Qs,Ang);
    %QsSL = -QsSL/max(abs(QsSL));
    QsglobalSL = interp1(angles,-Qsglobal,Ang);
    DiffSL = interp1(angles,Diff,Ang);
    DiffglobalSL = interp1(angles,Diffglobal,Ang);
    
    
    % compute dQs/dx
    grad(1:length(QsSL)) = 0; %alongshore gradient
    
    for i=1:length(QsSL)-1
        grad(i)=(QsSL(i+1)-QsSL(i)/(Xrotopen(i+1)-Xrotopen(i)));
    end
    
    grady=zeros(1,length(grad))'; %for plotting
    
    slope=zeros(1,length(Xrotopen));
    for i = 1+avewidth:openl-avewidth
        ii  = (i-avewidth): (i + avewidth);
        p2 = polyfit(Xrotopen(ii),QsSL(ii),1);
        slope(i) = p2(1);
    end
    
    
    if plot_on
    
    % Plot of Diffs
    
    scr = get(0,'DefaultFigurePosition');
    scr = scr*2;
    figure ('Outerposition',scr)
    % plot of the shoreline
    subplot(4,1,1:2)
    plot(Xrot, Yrot,'b','Linewidth',0.25)
    hold on
    plot(Xrot(ind_open),Yrot(ind_open),'-r','Linewidth',2)
%     xlim([Limxl Limxr])
%     ylim([Limyd Limyu])
    %set(gca,'linewidth',2)
    set(gca,'fontweight','bold')
    set(gca,'fontsize',18)
    daspect([1 1 1])
    ylabel('Cross shore (km)')
    f = get(gca,'TickLength');
    set(gca,'TickLength',f*3);
    hold off
    %
    
    
    subplot(4,1,3)
    % plot ast
    plot(Xrotopen,QsSL,'g','linewidth',1.5)
    hold on
    %plot(Xint,QsLvR,'g','linewidth',1.5)
    line(Xrot,zeros(length(Xrot),1),'Color','k')
    set(gca,'fontweight','bold')
    set(gca,'fontsize',18)
    %legend('Net','LvR')
    ylim([-1 1])
%     xlim([xliml xlimr])
    ylabel('Qs')
    hold off
    
    subplot(4,1,4)
    % plot diffs
    plot(Xrotopen,GammaSL,'-m','linewidth',1.5)
    hold on
    %plot(Yp,DiffFake,'b','linewidth',1.5)
    line(Xrot,zeros(length(Xrot),1),'Color','k')
    %legend('Deep','Local')
%     xlim([xliml xlimr])
    ylim([-1 1])
    %ylim([-0.4 0.6])
    %set(gca,'linewidth',2)
    set(gca,'fontweight','bold')
    set(gca,'fontsize',14)
    ylabel('Inst')
    xlabel('Alongshore location (km)')
    f = get(gca,'TickLength');
    set(gca,'TickLength',f*3);
    
    %M(mi) = getframe(gcf);
    %mi = mi + 1;
    %
    % figure(99)
    % subplot(2,1,1,'align')
    % plot(Xrot, Yrot,'b','Linewidth',0.25)
    % hold on
    % plot(Xrot(good),Yrot(good),'-r','Linewidth',2)
    % daspect([1 1 1]);
    % title('LG Waves')
    % xlim([0 200])
    % %ylim([-3040 -3000])
    %
    % subplot(2,1,2,'align')
    % %plot(Xrot(good),SLAng)
    % %subplot(2,1,2)
    % plot(Xrot(good(ii)),GammaSL,'--b','Linewidth',1.5)
    % hold on
    % plot(Xrot(good(ii)),QsSL,'-g','Linewidth',1.5)
    % line(Xrot(good), zeros(length(good),1),'color','k','linestyle','-','linewidth',1)
    % xlim([0 350])
    % legend('Gamma','Qs')
    

    
%     figure()
%     scatter(Ang,GammaSL,'g','.')
%     %ylim([-1 1])
%     ylabel('Gamma')
%     xlabel('Ang')
%     hold off

    
    figure()
    scatter(Ang,QsSL,'g','.')
    %ylim([-1 1])
    ylabel('Qs normalized by total')
    xlabel('Ang')
    title('Qs normalized by total')

    
    figure()
    scatter(Ang,QsglobalSL,'r')
    ylabel('Qs Not Normalized')
    xlabel('Ang')
    title('Qs Not Normalized')
%     [sum(AST,2) ASTSum];
    end
    save(savename)
    
    time_run = zeros(length(input_files),1);
    time_run(runs) = toc./60 %time each run takes in minutes
    
end
