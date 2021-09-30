function computewaveclimate_CERC(inputFolder_rose,outputFolder_qsdif,plot_on)
% this version computes CERC formula of Qs


% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isdir(inputFolder_rose)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s', inputFolder_rose);
    uiwait(warndlg(errorMessage));
    return;
end
if ~isdir(outputFolder_qsdif)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s', outputFolder_qsdif);
    uiwait(warndlg(errorMessage));
    return;
end

% Get a list of all files in the folder with the desired file name pattern.
infilePattern_rose = fullfile(inputFolder_rose, '*.mat'); % Change to whatever pattern you need.
input_files_rose = dir(infilePattern_rose);

for runs = 1 : length(input_files_rose)
    runs
    tic;
    baseFileName_rose = input_files_rose(runs).name;
    fullFileName_rose = fullfile(inputFolder_rose, baseFileName_rose);
    %fprintf(1, 'Now reading %s\n', fullFileName);
    load(fullFileName_rose);
    [pathstr_rose,filename_rose,ext_rose] = fileparts(fullFileName_rose)
    savename = [outputFolder_qsdif filename_rose '_qsdif.mat'];
    
    
    % variables
        
    
    orients = 0;         % shoreline orientation of normal to data array
    orientdel = 5;
    orientend = 180;
    calcangles = (orientend - orients)/orientdel + 1;
    shornt = 0;           % degrees (clockwise = positive) difference between location and array
    g= 9.81;             % That's gravity, holmes
    pi = 3.1415;         % ref. Archimedes (-232)
    degtorad = pi/180;   % convert radians to degrees
    depth = 200;          % m - depth @ WIS station - not really used here
    binsize = 7.5;       % um, size of the bins
    binnumber = 360/binsize;
    binarray = -180+binsize/2:binsize:180-binsize/2;
    plotbinarray = -90:2*binsize:90;
    
    % load data
    
    ldata = length(HInV);
    k = 0;
    
    EnergyFluxarray = zeros(calcangles,360/binsize);
    ShAngs = zeros(calcangles, 1);
    QsRaws = zeros(calcangles,6);
    Qss = zeros(calcangles, 6);
    QssGross = zeros(calcangles, 6);
    ASTs = zeros(calcangles,6,binnumber);
    Diffs = zeros(calcangles, 6, binnumber);
    Gammas = zeros(calcangles, 6);
    TrueDiffs = zeros(calcangles, 6);
    
    % where loop will start
    
    for ShAng = orients:orientdel:orientend
        [runs ShAng]
        k = k+1;
        
        % zero some stuff
        ShAng;
        
        EnergyFlux = zeros(1,360/binsize);
        
        AST = zeros(2,360/binsize);
        ASTSum = zeros(2,1);
        Diff = zeros(2,360/binsize);
        DiffSum = zeros(2,1);
        
        Bad = 0; % bad data points
        NotRef = 0;
        
        for j = 1:ldata
            
            H = HDeepV(j);
            Ang = (mod(AngleDeepV(j)-ShAng+180,360)-180);
            T = TIn(j);
            
            if ((abs(Ang)<90))
                for i = 1:binnumber
                    
                    if ((((i-1)*binsize)-180 < Ang) && (Ang < i * binsize-180))
                        
                        asign = sign(Ang);
                        Angle = abs(Ang);
                        
                        % add up contribution to energy flux
                        HTscale = (H^(12/5)) * (T^(1/5));
                        EnergyFlux(i) = EnergyFlux(i) + HTscale;
                        
                        % CERC formula
                        qs = -asign * HTscale * (cos(Angle*degtorad)^(1.2)) * (sin(Angle*degtorad));
                        AST(1,i) = AST(1,i) + qs;
                        ASTSum(1) = ASTSum(1) + abs(qs);
                        
                        rc =  -HTscale * (6/5 * ((sin(Angle*degtorad)^2) * (cos(Angle*degtorad)^(.2))) - cos(Angle*degtorad)^2.2 );
                        Diff(1,i) = Diff(1,i) + rc;
                        DiffSum(1) = DiffSum(1) + abs(rc);
                        
                    end
                end
            end
            H = HInV(j);
            Ang = (mod(AngleInV(j)-ShAng+180,360)-180);
            
            if ((abs(Ang)<90))
                
                for i = 1:binnumber
                    if ((((i-1)*binsize)-180 < Ang) && (Ang < i * binsize-180))
                        asign = sign(Ang);
                        Angle = abs(Ang);
                        
                        % add up contribution to energy flux
                        HTscale = (H^(12/5)) * (T^(1/5));
                        EnergyFlux(i) = EnergyFlux(i) + HTscale;
                        
                        % CERC formula
                        qs = -asign * HTscale * (cos(Angle*degtorad)^(1.2)) * (sin(Angle*degtorad));
                        AST(2,i) = AST(2,i) + qs;
                        ASTSum(2) = ASTSum(2) + abs(qs);
                        
                        rc =  -HTscale * (6/5 * ((sin(Angle*degtorad)^2) * (cos(Angle*degtorad)^(.2))) - cos(Angle*degtorad)^2.2 );
                        Diff(2,i) = Diff(2,i) + rc;
                        DiffSum(2) = DiffSum(2) + abs(rc);
                    end
                    
                end
            end
            
            if mod(j,1000) == 0
                j;
            end
            
        end
        
        
        Qs = sum(AST,2) ./ ASTSum;
        QsRaw = sum(AST,2);
        Gamma = sum(Diff,2) ./ DiffSum;
        TrueDiff = sum(Diff,2);
        
        
        QsRaws(k,:) = QsRaw;
        ShAngs(k) = ShAng;
        Qss(k,:) = Qs;
        QssGross(k,:) = ASTSum;
        Gammas(k,:) = Gamma;
        TrueDiffs(k,:) = TrueDiff;
        
        ASTs(k,:,:) = AST(:,:);
        Diffs(k,:,:) = Diff(:,:);
        EnergyFluxarray(k,:) = EnergyFlux(:)';
        
    end
    
    highangle = (abs(binarray) > 45);
    
    if plot_on
    figure()
    plot(ShAngs,Qss(:,1),'g','linewidth',1.5)
    hold on
    %XLim([-90 90])
    %set(gca,'xtick',[plotbinarray])
    line(ShAngs,zeros(calcangles,1),'color','k')
    title(['Qs ', filename]);
    legend('CERC')
    hold off
    
    figure()
    plot(ShAngs,Gammas(:,1),'g','linewidth',1.5)
    hold on
    %XLim([-90 90])
    %set(gca,'xtick',[plotbinarray])
    line(ShAngs,zeros(calcangles,1),'color','k')
    title(['Diffs ', filename]);
    legend('CERC')
    hold off
    end
    
    clear AngleDeepV HInV HDeepV AngleInV TIn
    save(savename)
    
    
    time_run = zeros(length(input_files_rose),1);
    time_run(runs) = toc./60 %time each run takes in minutes
end
end