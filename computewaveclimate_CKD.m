function computewaveclimate(inputFolder_rose,outputFolder_qsdif,plot_on)
% this version computes CERC, Kamphius, and Deigaard formulas of Qs


% this is the "deigaard look up table"
load DeigLook


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
        
    
    %ShoreOrient = 72;
    %ShelfOrient = 72;   % for backrefraction computations - should backrefract compard to shelf, not shore orientation
    % Then, ewith deep variables, bring the wves back in.
    
    %cutang = 72; %angle to cut data around - don't go past 180
    orients = 0;         % shoreline orientation of normal to data array
    orientdel = 5;
    orientend = 180;
    calcangles = (orientend - orients)/orientdel + 1
    shornt = 0           % degrees (clockwise = positive) difference between location and array
    g= 9.81;             % That's gravity, holmes
    pi = 3.1415;         % ref. Archimedes (-232)
    degtorad = pi/180;   % convert radians to degrees
    depth = 200;          % m - depth @ WIS station - not really used here
    binsize = 7.5;       % um, size of the bins
    binnumber = 360/binsize;
    binarray = -180+binsize/2:binsize:180-binsize/2;
    plotbinarray = -90:2*binsize:90;
    
    % load data
    % DO OUTFILE NAME, FOO
    
    ldata = length(HInV);
    
    %ff = zeros(ldata,5);
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
        
        AST = zeros(6,360/binsize);
        ASTSum = zeros(6,1);
        Diff = zeros(6,360/binsize);
        DiffSum = zeros(6,1);
        
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
                        
                        
                        % Kamphuis
                        HTscale = (H^1.84) * (T ^ 1.82);
                        qs = -asign * HTscale * (cos(Angle*degtorad)^0.92) * (sin(Angle*degtorad)^0.6);
                        
                        AST(2,i) = AST(2,i) + qs;
                        ASTSum(2) = ASTSum(2) + abs(qs);
                        
                        % reduce the importance of infinite answers
                        AngleKamp = Angle;
                        if (AngleKamp < 1)
                            AngleKamp = 1;
                        end
                        
                        t1 = 0.92 * (sin(AngleKamp*degtorad)^1.6) * (cos(AngleKamp*degtorad)^-0.08);
                        t2 = 0.6 * (sin(AngleKamp*degtorad)^-0.4) * (cos(AngleKamp*degtorad)^1.6);
                        rc = (T ^ 1.82) * (H ^ 1.84) * (t2 - t1);
                        Diff(2,i) = Diff(2,i) + rc;
                        DiffSum(2) = DiffSum(2) + abs(rc);
                        
                        % Deigaard
                        
                        HTscale = (H^3.8) * (T^-0.5);
                        qs = -asign * HTscale * interp1(DeigASTAngles,DeigAST,Angle);
                        AST(3,i) = AST(3,i) + qs;
                        ASTSum(3) = ASTSum(3) + abs(qs);
                        
                        rc = HTscale * interp1(DeigDiffAngles, DeigLookDiff, Angle);
                        Diff(3,i) = Diff(3,i) + rc;
                        DiffSum(3) = DiffSum(3) + abs(rc);
                        
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
                        AST(4,i) = AST(4,i) + qs;
                        ASTSum(4) = ASTSum(4) + abs(qs);
                        
                        rc =  -HTscale * (6/5 * ((sin(Angle*degtorad)^2) * (cos(Angle*degtorad)^(.2))) - cos(Angle*degtorad)^2.2 );
                        Diff(4,i) = Diff(4,i) + rc;
                        DiffSum(4) = DiffSum(4) + abs(rc);
                        
                        
                        % Kamphuis
                        HTscale = (H^1.84) * (T ^ 1.82);
                        qs = -asign * HTscale * (cos(Angle*degtorad)^0.92) * (sin(Angle*degtorad)^0.6);
                        
                        AST(5,i) = AST(5,i) + qs;
                        ASTSum(5) = ASTSum(5) + abs(qs);
                        
                        % reduce the importance of infinite answers
                        AngleKamp = Angle;
                        if (AngleKamp < 1)
                            AngleKamp = 1;
                        end
                        
                        t1 = 0.92 * (sin(AngleKamp*degtorad)^1.6) * (cos(AngleKamp*degtorad)^-0.08);
                        t2 = 0.6 * (sin(AngleKamp*degtorad)^-0.4) * (cos(AngleKamp*degtorad)^1.6);
                        rc = (T ^ 1.82) * (H ^ 1.84) * (t2 - t1);
                        Diff(5,i) = Diff(5,i) + rc;
                        DiffSum(5) = DiffSum(5) + abs(rc);
                        
                        % Deigaard
                        
                        HTscale = (H^3.8) * (T^-0.5);
                        qs = -asign * HTscale * interp1(DeigASTAngles,DeigAST,Angle);
                        AST(6,i) = AST(6,i) + qs;
                        ASTSum(6) = ASTSum(6) + abs(qs);
                        
                        rc = HTscale * interp1(DeigDiffAngles, DeigLookDiff, Angle);
                        Diff(6,i) = Diff(6,i) + rc;
                        DiffSum(6) = DiffSum(6) + abs(rc);
                        
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
        
        %        QsNB = sum(ASTNB)/ASTSumNB;
        %        QsGrossNB = ASTSum;
        %        DiffNB = sum(RateContNB)/RateSumNB;
        
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
    plot(ShAngs,Qss(:,2),'r','linewidth',1.5)
    plot(ShAngs,Qss(:,3),'b','linewidth',1.5)
    %XLim([-90 90])
    %set(gca,'xtick',[plotbinarray])
    line(ShAngs,zeros(calcangles,1),'color','k')
    title(['Qs ', filename]);
    legend('CERC','Kamp', 'Deig')
    hold off
    
    figure()
    plot(ShAngs,Gammas(:,1),'g','linewidth',1.5)
    hold on
    plot(ShAngs,Gammas(:,2),'r','linewidth',1.5)
    plot(ShAngs,Gammas(:,3),'b','linewidth',1.5)
    %XLim([-90 90])
    %set(gca,'xtick',[plotbinarray])
    line(ShAngs,zeros(calcangles,1),'color','k')
    title(['Diffs ', filename]);
    legend('CERC','Kamp', 'Deig')
    hold off
    end
    
    clear AngleDeepV HInV HDeepV AngleInV TIn
    save(savename)
    
    
    time_run = zeros(length(input_files_rose),1);
    time_run(runs) = toc./60 %time each run takes in minutes
end
end