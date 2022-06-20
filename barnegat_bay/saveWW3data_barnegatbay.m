function saveWW3data(inputFolder_WW3,outputFolder_WW3,ShelfOrient)

%this function takes wave watch data and outputs wave roses


% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isdir(inputFolder_WW3)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', inputFolder_WW3);
  uiwait(warndlg(errorMessage));
  return;
end
if ~isdir(outputFolder_WW3)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', outputFolder_WW3);
  uiwait(warndlg(errorMessage));
  return;
end

% Get a list of all files in the folder with the desired file name pattern.
infilePattern = fullfile(inputFolder_WW3, '*.mat'); % Change to whatever pattern you need.
input_files = dir(infilePattern);

for runs = 1 : length(input_files)
    baseFileName = input_files(runs).name;
    fullFileName = fullfile(inputFolder_WW3, baseFileName);
    %fprintf(1, 'Now reading %s\n', fullFileName);
    load(fullFileName);
    [pathstr,filename,ext] = fileparts(fullFileName)
    savename = [outputFolder_WW3 '/' filename '_rose.mat'];
    
    %loadname = ['C:/Andrew/WIS/WISDATA/WISn' num2str(stations(k)) '_80_99']
    %load(loadname);
    %load NDBC_44018_02_12(mod)

    %savename = ['climdata/rose' num2str(stations(k)) '_80_99_' 'th' num2str(Shores(k)) '_pk']
    ShoreOrient = ShelfOrient;

  
    
    %ShoreOrient = Shores(k);
    %ShelfOrient = Shores(k);   % for backrefraction computations - should backrefract compard to shelf, not shore orientation
                               % Then, ewith deep variables, bring the wves back in.

    shornt = 0;           % degrees (clockwise = positive) difference between location and array
    g= 9.81;             % That's gravity, holmes
    degtorad = pi/180;   % convert radians to degrees
    depth = abs(depth);
    binsize = 15;       % um, size of the bins
    binnumber = 360/binsize
    binarray = binsize/2:binsize:360-binsize/2;
    
    EnergyFluxDeep = zeros(1,360/binsize);
    EnergyFluxIn = zeros(1,360/binsize);

    ldata = length(data);

    AngleInV = zeros(ldata,1);
    AngleDeepV= zeros(ldata,1);
    HInV = zeros(ldata,1);
    HDeepV = zeros(ldata,1);
    TIn = zeros(ldata,1);  % ref coastal list

    % where loop will start

    ShAng = ShoreOrient;

    % zero some stuff
    ShAng
    ASTSum = 0;
    RateSum = 0;
    ASTSumNB = 0;
    RateSumNB = 0;
    ASTNB = zeros(1,360/binsize);
    RateContNB = zeros(1,360/binsize);
    EnergyFluxNB = zeros(1,360/binsize);
    AST = zeros(1,360/binsize);
    RateCont = zeros(1,360/binsize);
    EnergyFlux = zeros(1,360/binsize);
    HighAST = 0;
    HighASTSum = 0;
    HighRateAsym = 0;
    HighRateSum = 0;
    LowAST = 0;
    LowASTSum = 0;
    LowRateAsym = 0;
    LowRateSum = 0;

    Bad = []; % bad data points
    NotRef = [];

    for j = 1:ldata

        HIn = data(j,1);
        T = data(j,3);

        WaveAngleIn = data(j,2);
        AngleRot = (mod(data(j,2)-ShelfOrient+180,360)-180);

        % for all calculated, look only at what approaches a shore,
        % so +- 90 is all
        % add up contributions to
        if (abs(AngleRot)< 90)

            % backrefract
            Length =   g * T^2/2/pi*(tanh(((2* pi)^2 * depth / (T^2 * g))^(3/4)) ^ (2/3));
            % Ldeep = g * T^2/2/pi;
            C = g*T /2/pi * tanh(2*pi*depth/Length);
            Cdeep = g*T /2/pi;
            AngleDeepRot = sign(AngleRot)*asin(Cdeep/C*sin(abs(degtorad*AngleRot))) / degtorad;
            n = 0.5 * (1+ (4*pi*depth/Length)/sinh(4*pi*depth/Length));
            H = HIn * ((2 * n * C / Cdeep * cos(AngleRot*degtorad) / cos(AngleDeepRot*degtorad))^.5);

            if (~isreal(AngleDeepRot) | ~isreal(H))
                H = 0;
                AngleDeepRot = 995 * sign(AngleRot);
                Bad = [Bad j];
            end

            % However, we still need to count stuff that wasn't backrefracted.
            % NOTE THAT THIS COULD BE BAD IF NOT DEEP ALREADY
            % but not bad if shore-parallel but not deep, or whatevas, eh?
        else
            AngleDeepRot = AngleRot;
            H = HIn;
            NotRef = [NotRef j];
        end

        % Unrotate the angle so we can save it
        AngleDeep = mod(AngleDeepRot+ ShelfOrient,360);

        AngleDeepV(j) = AngleDeep;
        HDeepV(j) = H;
        % Save in data for comparison the backrefracted stuff - keep oriented
        % in general reference frame
        AngleInV(j) = WaveAngleIn;
        HInV(j) = HIn;
        TIn(j) = T;

        for i = 1:binnumber

            if ((((i-1)*binsize)< AngleDeep) & (AngleDeep <= i * binsize))
                % add up contribution to energy flux
                HTscale = (H^(12/5)) * (T^(1/5));
                EnergyFluxDeep(i) = EnergyFluxDeep(i) + HTscale;
            end
            if ((((i-1)*binsize) < WaveAngleIn) & (WaveAngleIn <= i * binsize))
                % add up contribution to energy flux
                HTscale = (HIn^(12/5)) * (T^(1/5));
                EnergyFluxIn(i) = EnergyFluxIn(i) + HTscale;
            end
        end

        if mod(j,1000) == 0
            j
        end
    end

    Bad;
    NotRef;
    BadData = length(Bad)
    NotRefData = length(NotRef)

    figure()
    plot(binarray,EnergyFluxDeep,'g','linewidth',1.5)
    hold on
    plot(binarray,EnergyFluxIn,'r','linewidth',1.5)
    legend('Averaged','Peak')
    title(filename)
    xlim([0 360])
    hold off

    figure()
    deepnormal = EnergyFluxDeep/sum(EnergyFluxDeep);
    polargeo([binarray binarray(1)]*degtorad, [deepnormal deepnormal(1)],'g')
    hold on
    innormal = EnergyFluxIn/sum(EnergyFluxIn);
    polargeo([binarray binarray(1)]*degtorad, [innormal innormal(1)],'b')
    
    title(filename)
    legend('Deep','In')
    h= findobj('Color','r');
    set(h,'MarkerFaceColor','r')
    h= findobj('Color','g');
    set(h,'LineWidth',2)
    hold off

    save(join(savename,''))
    %savename
end

end