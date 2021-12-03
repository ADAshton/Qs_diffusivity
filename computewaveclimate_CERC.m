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
infilePattern_rose = fullfile(inputFolder_rose, '*_rose.mat'); % Change to whatever pattern you need.
input_files_rose = dir(infilePattern_rose);

for runs = 1 : length(input_files_rose)
    runs
    tic;
    baseFileName_rose = input_files_rose(runs).name;
    fullFileName_rose = fullfile(inputFolder_rose, baseFileName_rose);
    %fprintf(1, 'Now reading %s\n', fullFileName);
    load(fullFileName_rose);
    [~,filename_rose,~] = fileparts(fullFileName_rose);
    savename = [outputFolder_qsdif '/' filename_rose '_qsdif.mat'];
    
    
    % variables
        
    
    orients = 0;         % shoreline orientation of normal to data array
    orientdel = 5;
    orientend = 180;
    calcangles = (orientend - orients)/orientdel + 1;
    shornt = 0;           % degrees (clockwise = positive) difference between location and array
    g= 9.81;             % That's gravity, holmes
    pi = 3.1415;         % ref. Archimedes (-232)
    degtorad = pi/180;   % convert radians to degrees
    binsize = 7.5;       % um, size of the bins
    binnumber = 360/binsize;
    binarray = -180+binsize/2:binsize:180-binsize/2;
    plotbinarray = -90:2*binsize:90;
    
    % load data
    
    ldata = length(HInV);
    k = 0;
    
    EnergyFluxarray = zeros(calcangles,360/binsize);
    ShAngs = zeros(calcangles, 1);
    QsRaws = zeros(calcangles,1);
    Qss = zeros(calcangles, 1);
    QssGross = zeros(calcangles, 1);
    ASTs = zeros(calcangles,1,binnumber);
    Diffs = zeros(calcangles, 1, binnumber);
    Gammas = zeros(calcangles, 1);
    TrueDiffs = zeros(calcangles, 1);
    
    % where loop will start
    
    for ShAng = orients:orientdel:orientend
        [runs ShAng]
        k = k+1;
        
        % zero some stuff
        ShAng;
        
        EnergyFlux = zeros(1,360/binsize);
        
        AST = zeros(1,360/binsize);
        ASTSum = zeros(1,1);
        Diff = zeros(1,360/binsize);
        DiffSum = zeros(1,1);
        
        Bad = 0; % bad data points
        NotRef = 0;
        
        for j = 1:ldata
            
            H = HInV(j);
%             Ang = (mod(AngleInV(j)-ShAng+180,360)-180);
            if AngleInV(j)>180 % to use the wave crest orientation -/+ 90 from the orientation reported by WIS...
                Ang = (mod(AngleInV(j)+90-ShAng+180,360)-180);
            else
                Ang = (mod(AngleInV(j)-90-ShAng+180,360)-180);
            end
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
                        AST(i) = AST(i) + qs;
                        ASTSum = ASTSum + abs(qs);
                        
                        rc =  -HTscale * (6/5 * ((sin(Angle*degtorad)^2) * (cos(Angle*degtorad)^(.2))) - cos(Angle*degtorad)^2.2 );
                        Diff(i) = Diff(i) + rc;
                        DiffSum = DiffSum + abs(rc);
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
    
    clear HInV AngleInV TIn
    save(join(savename,''))
    
    
end
end