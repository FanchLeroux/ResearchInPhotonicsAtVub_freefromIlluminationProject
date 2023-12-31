% fleroux - 10/25/2023

%% main programm

function [ r ] = finiteDistanceGaussianSource( args )

    clc; format compact; clear;

    if ~exist('args', 'var')
        args = [];
    end

    % Initialize the OpticStudio connection
    TheApplication = InitConnection();
    if isempty(TheApplication)
        % failed to initialize a connection
        r = [];
    else
        try
            r = BeginApplication(TheApplication, args);
            CleanupConnection(TheApplication);
        catch err
            CleanupConnection(TheApplication);
            rethrow(err);
        end   
    end
    
end

%% body

function [r] = BeginApplication(TheApplication, ~)

    % 8<----------------- Define directories and file names ----------------->8

    dirc = "D:\moi\vub\researchInPhotonics\zemax\zosApi\gaussianSource\circularTopHat\"; % directory where the ZEMAX .zos file will be generated
    resultDir = dirc + "results\";
    
    zemaxFileName = dirc + "inputFiniteDistanceGaussianSource_outputCircularUniformIrradiance.zos";    
    cfg1FileName = dirc + "config\geometricImageAnalysis1.cfg";
    cfg2FileName = dirc + "config\geometricImageAnalysis2.cfg";

    % 8<--------------------- Define system parameters ---------------------->8

    % wavelength
    lambda = 0.633; % wavelength [um]

    % input and output irradiance distributions
    w = 5.0; % gaussian input beam waist [mm]
    k = 25.0; % output beam radius [mm]

    % optimization
    nParMax = 6; % number of aspheric coefficients set as variables
    sample = 100; % pupil sampling for the ray-mapping function targets computations

    % analysis
    nRays = 5000000; % number of rays for geometrical image analysis (typical: 5000000)
    imageSize = 100; % image size for geometrical image analysis

    % system start design. note that all surfaces are flat before optimization
    distanceSourceLens = 50; % [mm] distance between Source and freeform lens entrance facet (works with 50)
    apertureAngle = 16.6; % [°] half aperture angle (works with 16.6)
    objectSpaceNA = sin(pi/180 * apertureAngle);
    apodizationFactor = 9; %1/(w/(entrancePupilDiameter/2))^2;
    backFocalLength = 70;

    % 8<-------------------------- Create Zemax file ------------------------>8

    import ZOSAPI.*;

    % creates a new API directory
    apiPath = System.String.Concat(TheApplication.SamplesDir, '\API\Matlab');
    if (exist(char(apiPath), 'dir') == 0) 
        mkdir(char(apiPath)); 
    end

    % Set up primary optical system
    TheSystem = TheApplication.PrimarySystem;

    % Make new file
    testFile = zemaxFileName;
    TheSystem.New(false);
    TheSystem.SaveAs(testFile);

    TheSystem.SystemData.MaterialCatalogs.AddCatalog('SCHOTT');

    % 8<---------------------- Build analysis tools ------------------------->8

    % analysis 1: spot diagram
    analysis1 = TheSystem.Analyses.New_Analysis(ZOSAPI.Analysis.AnalysisIDM.GeometricImageAnalysis);
    analysis1Settings = analysis1.GetSettings();
    analysis1Settings.ShowAs = ZOSAPI.Analysis.GiaShowAsTypes.FalseColor; % make sure to perform this step before saving the settings in a configuration file
    analysis1Settings.SaveTo(cfg1FileName);
    analysis1Settings.ModifySettings(cfg1FileName, 'IMA_KRAYS', string(nRays/1000));
    analysis1Settings.ModifySettings(cfg1FileName, 'IMA_IMAGESIZE', string(imageSize));
    analysis1Settings.LoadFrom(cfg1FileName);
    
    % analysis 2: crossX
    analysis2 = TheSystem.Analyses.New_Analysis(ZOSAPI.Analysis.AnalysisIDM.GeometricImageAnalysis);
    analysis2Settings = analysis2.GetSettings();
    analysis2Settings.ShowAs = ZOSAPI.Analysis.GiaShowAsTypes.CrossX;
    analysis2Settings.SaveTo(cfg2FileName);
    analysis2Settings.ModifySettings(cfg2FileName, 'IMA_KRAYS', string(nRays/1000));
    analysis2Settings.ModifySettings(cfg2FileName, 'IMA_IMAGESIZE', string(imageSize));
    analysis2Settings.LoadFrom(cfg2FileName);
    
    % analysis 3: 3D layout    
    analysis3 = TheSystem.Analyses.New_Analysis(ZOSAPI.Analysis.AnalysisIDM.Draw3D);

    % 8<--------------------- Build results containers ---------------------->8
    
    stdVect = zeros(nParMax,1);
    crossXvect = zeros(nParMax, imageSize);
    
    % 8<------------------ System explorer parameters ----------------------->8    

    % Aperture    
    TheSystemData = TheSystem.SystemData;
    TheSystemData.Aperture.ApertureType = ZOSAPI.SystemData.ZemaxApertureType.ObjectSpaceNA;
    TheSystemData.Aperture.ApertureValue = objectSpaceNA;
    TheSystemData.Aperture.ApodizationType = ZOSAPI.SystemData.ZemaxApodizationType.Gaussian;    
    TheSystemData.Aperture.ApodizationFactor = apodizationFactor;
    TheSystemData.Aperture.SetCurrentGCRSSurf(2); % Set Surface 2 as the Global Coordinate Reference Surface

    % Wavelength
    TheSystemData.Wavelengths.RemoveWavelength(1);
    TheSystemData.Wavelengths.AddWavelength(lambda, 1.0);

    % 8<------------------ Lens data editor parameters ---------------------->8   

    TheLDE = TheSystem.LDE;
    TheLDE.InsertNewSurfaceAt(2);
    TheLDE.InsertNewSurfaceAt(3);
    Surface_0 = TheLDE.GetSurfaceAt(0);
    Surface_1 = TheLDE.GetSurfaceAt(1);
    Surface_2 = TheLDE.GetSurfaceAt(2);
    Surface_3 = TheLDE.GetSurfaceAt(3);

    Surface_0.Thickness = distanceSourceLens-5; % -5 term because of the dummy surface between the lens and the object
    Surface_1.Thickness = 5.0;
    Surface_1.Comment = 'dummy';
    Surface_2.Thickness = 30;
    Surface_2.Comment = 'front of lens';
    Surface_2.Material = 'N-BK7';
    Surface_3.Thickness = backFocalLength;
    Surface_3.Comment = 'rear of lens';      

    % set surface 3 type as even aspheric
    SurfaceType_EvenAspheric = Surface_3.GetSurfaceTypeSettings(ZOSAPI.Editors.LDE.SurfaceType.EvenAspheric);
    Surface_3.ChangeType(SurfaceType_EvenAspheric);

    % set stop
    Surface_2.IsStop = true;

    % 8<----------- Build merit function using ray mapping function --------->8

    TheMFE = TheSystem.MFE;

    % get entrance pupil diameter
    entrancePupilDiameter = Surface_2.GetCellAt(6).Value;
    entrancePupilDiameter = 2*str2double(char(entrancePupilDiameter)); % fix data type issues and converting the radius in diameter    

    % ray-maping function
    for j = 1:sample

        Operand_j = TheMFE.InsertNewOperandAt(j);
        Operand_j.ChangeType(ZOSAPI.Editors.MFE.MeritOperandType.REAY);
        Operand_j.Weight = 1.0;

        normalizedPupilCoordinate = j/sample;
        pupilCoordinate = normalizedPupilCoordinate*entrancePupilDiameter/2; % points along the pupil radius
        target = -k*sqrt(1-exp(-2*pupilCoordinate^2/w^2)); % annalytical ray-mapping function for circular uniform illumination from gaussian input

        Operand_j.Target = target;

        Operand_1_SurfCell = Operand_j.GetCellAt(2);
        Operand_1_SurfCell.IntegerValue = 4;

        Operand_1_PyCell = Operand_j.GetCellAt(7);
        Operand_1_PyCell.Value = string(normalizedPupilCoordinate);
        
    end

    % constraints
    edgeConstraint = TheMFE.InsertNewOperandAt(sample+1);
    edgeConstraint.ChangeType(ZOSAPI.Editors.MFE.MeritOperandType.MNEG);

    edgeConstraint_Surf1cell = edgeConstraint.GetCellAt(2);
    edgeConstraint_Surf1cell.IntegerValue = 2;
    edgeConstraint_Surf2cell = edgeConstraint.GetCellAt(3);
    edgeConstraint_Surf2cell.IntegerValue = 3;

    edgeConstraint.Weight = sample;
    edgeConstraint.Target = 1;

    % 8<-------------------- Define variable parameters --------------------->8

    % set radius of surface 3 variable
    Surface_3.RadiusCell.MakeSolveVariable();

    % set conic constant of surface 3 variable
    Surface_3.ConicCell.MakeSolveVariable();

    % loop across nPar: optimization for different number of variables to define Aspheric surface
    for nPar = 1:nParMax
        
        % set the nPar first aspheric coefficients of surface 2 as variables
        for j = 2:nPar % Par1 is 2nd order aspheric coefficient, fixed to 0 because conic constant is variable
            Surface_3_2kThOrderTermCell = Surface_3.GetSurfaceCell(ZOSAPI.Editors.LDE.SurfaceColumn.("Par"+string(j)));
            Solver = Surface_3_2kThOrderTermCell.CreateSolveType(ZOSAPI.Editors.SolveType.Variable);
            Surface_3_2kThOrderTermCell.SetSolveData(Solver);
        end

        % 8<-------------------------- Optimization ------------------------->8
        
        now = tic();
        LocalOpt = TheSystem.Tools.OpenLocalOptimization();
        if ~isempty(LocalOpt)
            LocalOpt.Algorithm = ZOSAPI.Tools.Optimization.OptimizationAlgorithm.DampedLeastSquares;
            LocalOpt.Cycles = ZOSAPI.Tools.Optimization.OptimizationCycles.Automatic;
            LocalOpt.NumberOfCores = 8;
            fprintf('Local Optimization...\n');
            fprintf('Initial Merit Function %6.3f\n', LocalOpt.InitialMeritFunction);
            LocalOpt.RunAndWaitForCompletion();
            fprintf('Final Merit Function %6.3f\n', LocalOpt.CurrentMeritFunction);
            LocalOpt.Close();
        end
        optimizationTime = toc(now);
        fprintf('Optimization took %6.4f seconds\n', optimizationTime)
        
        % 8<---------------------------- Analysis --------------------------->8
        
        now = tic();
        analysis1.Terminate();
        analysis1.WaitForCompletion();
        analysis1time = toc(now);
        fprintf('Analysis 1 took %6.4f seconds\n', analysis1time)

        now = tic();
        analysis2.Terminate();
        analysis2.WaitForCompletion();
        analysis2time = toc(now);
        fprintf('Analysis 2 took %6.4f seconds\n', analysis2time)
        
        % 8<------------------------- Save results -------------------------->8
        
        % write into text file
        results1FileName = resultDir + ...
        "falseColorAnalisis_" + ...
        "nPar=" + string(nPar) + ...
        ".txt";
        results2FileName = resultDir + ...
        "crossX_" + ...
        "nPar=" + string(nPar) + ...
        ".txt";    
        results1 = analysis1.GetResults();
        results1.GetTextFile(results1FileName);
        results2 = analysis2.GetResults();
        results2.GetTextFile(results2FileName);
        
        % 8<------------------------- Read results -------------------------->8
        
        data2 = readmatrix(results2FileName);
        
        % 8<------------------------- Process results ----------------------->8
        
        dataCrossX = data2(:,2);
        crossXvect(nPar, :) = dataCrossX;
        standardDeviation = std(dataCrossX(26:75));
        stdVect(nPar) = standardDeviation;

    end
    
    % add 3D layout
    
    now = tic();
    analysis3.ApplyAndWaitForCompletion()
    analysis3time = toc(now);
    fprintf('Analysis 3 took %6.4f seconds\n', analysis3time)

    % read results

    data1 = readmatrix(results1FileName);
    data2 = readmatrix(results2FileName);

    % figures

    figure(1)
    imagesc(data1)
    axis equal
    title("Image plane irradiance map nPar="+"nParMax")

    % profile as a fuction of nPar
    figure(2)
    hold on
    legendCell = NaN(1,nParMax-1);
    for nPar = 2:nParMax
        plot(crossXvect(nPar,:))
        legendCell(nPar-1)=nPar;
    end
    legend(string(legendCell))
    title("CrossX irradiance profile as a fuction of nPar")
    
    % std as a fuction of nPar    
    figure(3)
    plot(stdVect(2:end), '+')
    title("std as a fuction of nPar")

    % 8<--------------------- Save and close Zemax file --------------------->8
    
    TheSystem.Save();
    r = [];
    
end

%% else

function app = InitConnection()

import System.Reflection.*;

% Find the installed version of OpticStudio.
zemaxData = winqueryreg('HKEY_CURRENT_USER', 'Software\Zemax', 'ZemaxRoot');
NetHelper = strcat(zemaxData, '\ZOS-API\Libraries\ZOSAPI_NetHelper.dll');
% Note -- uncomment the following line to use a custom NetHelper path
% NetHelper = 'C:\Users\Documents\Zemax\ZOS-API\Libraries\ZOSAPI_NetHelper.dll';
NET.addAssembly(NetHelper);

success = ZOSAPI_NetHelper.ZOSAPI_Initializer.Initialize();
% Note -- uncomment the following line to use a custom initialization path
% success = ZOSAPI_NetHelper.ZOSAPI_Initializer.Initialize('C:\Program Files\OpticStudio\');
if success == 1
    LogMessage(strcat('Found OpticStudio at: ', char(ZOSAPI_NetHelper.ZOSAPI_Initializer.GetZemaxDirectory())));
else
    app = [];
    return;
end

% Now load the ZOS-API assemblies
NET.addAssembly(AssemblyName('ZOSAPI_Interfaces'));
NET.addAssembly(AssemblyName('ZOSAPI'));

% Create the initial connection class
TheConnection = ZOSAPI.ZOSAPI_Connection();

% Attempt to create a Standalone connection

% NOTE - if this fails with a message like 'Unable to load one or more of
% the requested types', it is usually caused by try to connect to a 32-bit
% version of OpticStudio from a 64-bit version of MATLAB (or vice-versa).
% This is an issue with how MATLAB interfaces with .NET, and the only
% current workaround is to use 32- or 64-bit versions of both applications.

app = TheConnection.CreateNewApplication();
if isempty(app)
   HandleError('An unknown connection error occurred!');
end
if ~app.IsValidLicenseForAPI
    HandleError('License check failed!');
    app = [];
end

end

function LogMessage(msg)
disp(msg);
end

function HandleError(error)
ME = MXException(error);
throw(ME);
end

function  CleanupConnection(TheApplication)
% Note - this will close down the connection.

% If you want to keep the application open, you should skip this step
% and store the instance somewhere instead.
TheApplication.CloseApplication();
end