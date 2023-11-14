% fleroux - 10/25/2023

%% main programm
function [ r ] = finiteDistanceGaussianSource( args )

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

% 8<------------------------ Define directories and file names ------------------------------>8

dirc = 'D:\moi\vub\researchInPhotonics\zemax\zosApi\'; % directory where the ZEMAX .zos file will be generated
zemaxFileName = 'inputFiniteDistanceGaussianSource_outputCircularUniformIrradiance.zos';

cfg1FileName = "D:\moi\vub\researchInPhotonics\zemax\zosApi\config\geometricImageAnalysis1.cfg";
cfg2FileName = "D:\moi\vub\researchInPhotonics\zemax\zosApi\config\geometricImageAnalysis2.cfg";
results1FileName = "D:\moi\vub\researchInPhotonics\zemax\zosApi\results\results1.txt";
results2FileName = "D:\moi\vub\researchInPhotonics\zemax\zosApi\results\results2.txt";

% 8<---------------------------- Define system parameters ----------------------------------->8

    % system design
distanceSourceLens = 50; % [mm] distance between Source and freeform lens entrance facet 
apertureAngle = 16.6; % [°]
objectSpaceNA = sin(pi/180 * apertureAngle);
apodizationFactor = 9; %1/(w/(entrancePupilDiameter/2))^2;
backFocalLength = 70;

    % input and output irradiance distributions
w = 5.0; % gaussian input beam waist [mm]
k = 25.0; % output beam radius [mm]

    % wavelength
lambda = 0.633; % wavelength [um]

    % optimization
nPar = 8; % number of aspheric coefficients set as variables
sample = 80; % pupil sampling for the ray-mapping function targets computations

    % analysis
nRays = 5000; % number of rays for geometrical image analysis (typical: 5000000)
imageSize = 100; % image size for geometrical image analysis

% 8<----- ############################################################## ------->8

import ZOSAPI.*;
    
    % creates a new API directory
    apiPath = System.String.Concat(TheApplication.SamplesDir, '\API\Matlab');
    if (exist(char(apiPath), 'dir') == 0) 
        mkdir(char(apiPath)); 
    end
    
    % Set up primary optical system
    TheSystem = TheApplication.PrimarySystem;
        
    % Make new file
    testFile = System.String.Concat(dirc, zemaxFileName);
    TheSystem.New(false);
    TheSystem.SaveAs(testFile);
    
    TheSystem.SystemData.MaterialCatalogs.AddCatalog('SCHOTT');
    
    % Aperture
    TheSystemData = TheSystem.SystemData;
    TheSystemData.Aperture.ApertureType = ZOSAPI.SystemData.ZemaxApertureType.ObjectSpaceNA;
    TheSystemData.Aperture.ApertureValue = objectSpaceNA;
    TheSystemData.Aperture.ApodizationType = ZOSAPI.SystemData.ZemaxApodizationType.Gaussian;    
    TheSystemData.Aperture.ApodizationFactor = apodizationFactor;
    
    % Set Surface 2 as the Global Coordinate Reference Surface
    TheSystemData.Aperture.SetCurrentGCRSSurf(2);
    
    % Set Wavelength
    TheSystemData.Wavelengths.RemoveWavelength(1);
    TheSystemData.Wavelengths.AddWavelength(lambda, 1.0);
    
    % Lens data 
    TheLDE = TheSystem.LDE;
    TheLDE.InsertNewSurfaceAt(2);
    TheLDE.InsertNewSurfaceAt(3);
    Surface_0 = TheLDE.GetSurfaceAt(0);
    Surface_1 = TheLDE.GetSurfaceAt(1);
    Surface_2 = TheLDE.GetSurfaceAt(2);
    Surface_3 = TheLDE.GetSurfaceAt(3);
    
    % Changes surface cells in LDE
    
    Surface_0.Thickness = distanceSourceLens-5; % -5 term because of the dummy surface between the lens and the object
    Surface_1.Thickness = 5.0;
    Surface_1.Comment = 'dummy';
    Surface_2.Thickness = 30;
    Surface_2.Comment = 'front of lens';
    Surface_2.Material = 'N-BK7';
    Surface_3.Thickness = backFocalLength;
    Surface_3.Comment = 'rear of lens';      
    
    % set surface 2 type as even aspheric
    SurfaceType_EvenAspheric = Surface_2.GetSurfaceTypeSettings(ZOSAPI.Editors.LDE.SurfaceType.EvenAspheric);
    Surface_2.ChangeType(SurfaceType_EvenAspheric);
    
    % set stop
    Surface_2.IsStop = true;
    
    % get entrance pupil diameter
    entrancePupilDiameter = Surface_2.GetCellAt(6).Value;
    entrancePupilDiameter = 2*str2double(char(entrancePupilDiameter)); % fix data type issues and converting the radius in diameter
    
    % 8<--------------------------- Set variables ----------------------------------------------->8
    
    % set radius of surface 2 variable
    Surface_2.RadiusCell.MakeSolveVariable();
    
    % set conic constant of surface 2 variable
    Surface_2.ConicCell.MakeSolveVariable();
    
    % set the nPar first aspheric coefficients of surface 2 as variables
    for j = 2:nPar % Par1 is 2nd order aspheric coefficient, fixed to 0 because conic constant is variable
        Surface_2_2kThOrderTermCell = Surface_2.GetSurfaceCell(ZOSAPI.Editors.LDE.SurfaceColumn.("Par"+string(j)));
        Solver = Surface_2_2kThOrderTermCell.CreateSolveType(ZOSAPI.Editors.SolveType.Variable);
        Surface_2_2kThOrderTermCell.SetSolveData(Solver);
    end
    
    % merit function
    
    TheMFE = TheSystem.MFE;
    
        % ray-mapping
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
        
        edgeConstraint.Weight = 0;
        edgeConstraint.Target = 5;
        
    % optimize
    
    tic;
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
    toc;
    
    % analyse
    
        % analysis 1: spot diagram
    analysis1 = TheSystem.Analyses.New_Analysis(ZOSAPI.Analysis.AnalysisIDM.GeometricImageAnalysis);
        % change analysis settings
    analysis1Settings = analysis1.GetSettings();
    analysis1Settings.ShowAs = ZOSAPI.Analysis.GiaShowAsTypes.SpotDiagram; % make sure to perform this step before saving the settings in a configuration file
    analysis1Settings.SaveTo(cfg1FileName);
    analysis1Settings.ModifySettings(cfg1FileName, 'IMA_KRAYS', string(nRays/1000));
    analysis1Settings.ModifySettings(cfg1FileName, 'IMA_IMAGESIZE', string(imageSize));
    
    analysis1Settings.LoadFrom(cfg1FileName);
    
    tic;
    analysis1.Terminate();
    analysis1.WaitForCompletion();
    toc;
    
    % save results under text file
    results1 = analysis1.GetResults();
    results1.GetTextFile(results1FileName);
    
    % read results text file
    
    data1 = readmatrix(results1FileName);
    
        % analysis 2: CrossX
    analysis2 = TheSystem.Analyses.New_Analysis(ZOSAPI.Analysis.AnalysisIDM.GeometricImageAnalysis);
    % change analysis settings
    analysis2Settings = analysis2.GetSettings();
    analysis2Settings.ShowAs = ZOSAPI.Analysis.GiaShowAsTypes.CrossX;
    analysis2Settings.SaveTo(cfg2FileName);
    analysis2Settings.ModifySettings(cfg2FileName, 'IMA_KRAYS', string(nRays/2000));
    analysis2Settings.ModifySettings(cfg2FileName, 'IMA_IMAGESIZE', string(imageSize));
    analysis2Settings.LoadFrom(cfg2FileName);
    
    tic;
    analysis2.Terminate();
    analysis2.WaitForCompletion();
    toc;
    
    % save results under text file
    
    results2 = analysis2.GetResults();
    results2.GetTextFile(results2FileName);
    
    data2 = readmatrix(results2FileName);
    
    % figures
    
    figure(1)
    plot(data1(:,8), data1(:,9), '+')
    axis equal
    title("Image plane positional spot diagram")
    
    % Uniformity computation across a line (horizontal direction) - tbd
    dataCrossX = data2(:,2);
    figure(2)
    plot(dataCrossX, '+')
    
    % std computation
    
    standardDeviation = std(dataCrossX(25:76));
    
    % Save and close
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


