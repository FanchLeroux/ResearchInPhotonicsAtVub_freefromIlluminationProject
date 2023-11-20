% fleroux - 2023/11/20

% 8<----------------- Define directories and file names ----------------->8

dirc = "D:\moi\vub\researchInPhotonics\zemax\zosApi\lambertianPointSource\circularTopHat\";
resultDir = dirc + "results\";
rayMappingFunctionDir = dirc + "rayMappingFunctions\fromAlejandro20231120\";

name = "idkName";

% 8<------------------- Import ray mapping function --------------------->8

load(rayMappingFunctionDir+"M1.mat", "M1");
load(rayMappingFunctionDir+"ZM1.mat", "ZM1");
load(rayMappingFunctionDir+"M2.mat", "M2");
load(rayMappingFunctionDir+"ZM2.mat", "ZM2");
load(rayMappingFunctionDir+"Zx.mat", "Zx");
load(rayMappingFunctionDir+"Zy.mat", "Zy");

% 8<----------- Build merit function using ray mapping function --------->8

surf = 7;
oper = strings([2*size(ZM1,1),1]);
oper(1:size(ZM1,1),1) = 'REAY';
oper(size(ZM1,1)+1:end,1) = 'REAX';

vacio = zeros(2*size(ZM1,1),1);

% Real target size units in mm
Xhalf = 3*1000;
Yhalf = 3*1000;

% MTF: Oper, surf, wave,Hx,Hy,Px,Py,vacio,vacio,target,weight,vacio
T = table(oper,surf.*ones(2*size(ZM1,1),1),ones(2*size(ZM1,1),1),vacio ...
    , vacio, [Zx; Zx], [Zy; Zy], vacio, vacio, [Yhalf.*ZM2;Xhalf.*ZM1]...
    ,ones(2*size(ZM1,1),1), vacio);

writetable(T,strcat(resultDir,strcat('MTF_',name,'.dat')),'Delimiter','\t','WriteRowNames',false);
