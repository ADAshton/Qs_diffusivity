function saveWISdata(inputFolder_WIS,outputFolder_WIS)

% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isdir(inputFolder_WIS)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', inputFolder_WIS);
  uiwait(warndlg(errorMessage));
  return;
end
if ~isdir(outputFolder_WIS)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', outputFolder_WIS);
  uiwait(warndlg(errorMessage));
  return;
end

% Get a list of all files in the folder with the desired file name pattern.
infilePattern = fullfile(inputFolder_WIS, '*.mat'); % Change to whatever pattern you need.
input_files = dir(infilePattern);

for k = 1 : length(input_files)
    baseFileName = input_files(k).name;
    fullFileName = fullfile(inputFolder_WIS, baseFileName);
    %fprintf(1, 'Now reading %s\n', fullFileName);
    load(fullFileName);
    [pathstr,filename,ext] = fileparts(fullFileName)
    savename = [outputFolder_WIS filename '_rose.mat'];
    
    %loadname = ['C:/Andrew/WIS/WISDATA/WISn' num2str(stations(k)) '_80_99']
    %load(loadname);
    %load NDBC_44018_02_12(mod)

    %savename = ['climdata/rose' num2str(stations(k)) '_80_99_' 'th' num2str(Shores(k)) '_pk']

    shornt = 0;           % degrees (clockwise = positive) difference between location and array
    g= 9.81;             % That's gravity, holmes
    pi = 3.1415;         % ref. Archimedes (-232)
    degtorad = pi/180;   % convert radians to degrees
    binsize = 15;       % um, size of the bins
    binnumber = 360/binsize
    binarray = binsize/2:binsize:360-binsize/2;
    
    EnergyFluxIn = zeros(1,360/binsize);

    ldata = length(data);

    AngleInV = zeros(ldata,1);
    HInV = zeros(ldata,1);
    TIn = zeros(ldata,1);  % ref coastal list

    % where loop will start


    % zero some stuff
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

        % Save in data for comparison the backrefracted stuff - keep oriented
        % in general reference frame
        AngleInV(j) = WaveAngleIn;
        HInV(j) = HIn;
        TIn(j) = T;

        for i = 1:binnumber
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
    plot(binarray,EnergyFluxIn,'r','linewidth',1.5)
    legend('Averaged','Peak')
    title(filename)
    xlim([0 360])
    hold off

    figure()
    innormal = EnergyFluxIn/sum(EnergyFluxIn);
    polargeo([binarray binarray(1)]*degtorad, [innormal innormal(1)],'b')
    
    title(filename)
    h= findobj('Color','r');
    set(h,'MarkerFaceColor','r')
    h= findobj('Color','g');
    set(h,'LineWidth',2)
    hold off

    save(join(savename,''))
    %savename
end

end