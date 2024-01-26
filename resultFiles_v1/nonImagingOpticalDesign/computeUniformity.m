clear; format compact; clc;

% Optical Design non-Imaging - Results analysis

dirc = "D:\moi\vub\researchInPhotonics\zemax\zosApi\resultFiles_v1\nonImagingOpticalDesign\irradiancesDistributions\";
filenameDisk = "disk.txt";
filenameRectangle = "rectangle.txt";

irradianceDisk = readmatrix(dirc + filenameDisk);
irradianceDisk = irradianceDisk(2:end,2:end); % get rid of pixels indices on firt row / column

irradianceRectangle = readmatrix(dirc + filenameRectangle);
irradianceRectangle = irradianceRectangle(2:end,2:end); % get rid of pixels indices on firt row / column

diskUniformity = uniformity(irradianceDisk);
diskRectangle = uniformity(irradianceRectangle);

figure(1);
imagesc(irradianceDisk)
axis equal

figure(2);
imagesc(irradianceRectangle)
axis equal

function u = uniformity(irradianceDistribution)

    N = sum(sum(irradianceDistribution>0));
    meanIrradiance = sum(sum(irradianceDistribution))/N;
    
    u = 1 - sum(abs(irradianceDistribution(irradianceDistribution>0)-meanIrradiance))/(N*meanIrradiance);

end



