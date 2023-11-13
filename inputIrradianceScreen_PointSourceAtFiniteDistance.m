% fleroux - 2023/11/13

clear; clc; format compact

% This script aim at computing the irradiance pattern produced by a point source over a flat 
% surface (input plane) at finite distance.
% The obtained irradiance distribution can then be used to compute numerically a
% ray-mapping function for a given targeted irradiance distrubution in a given
% output plane, and perform a freeform lens optimization based on this
% ray-mapping function.

dirc = "D:\moi\vub\researchInPhotonics\zemax\zosApi\results\";

% Parameters

n=1; % index of refraction

intensity = 1; % [W/sr] intensity of the point source

theta=82; % ° half angle of the maximum cone of light emitted by the LED reaching the input plane 
distanceLedInputPlane = 6e-3; % [m]

inputPlaneSampling = 100; % nuber of pixels in one pupil diameter

% Consequences
inputPlaneDiameter = 2*tan(pi/180 * theta) * distanceLedInputPlane;
inputPlaneSamplingStep = inputPlaneDiameter/inputPlaneSampling; % physical length of one pixel in the input plane 
numericalAperture = n*sin(pi/180 * theta); % object-space numerical aperture

% Irradiance screen computation
r = mod(inputPlaneSampling, 2);
q = (inputPlaneSampling-r)/2;
[x,y] = meshgrid((-q:q-1+r),(q-1+r:-1:-q)); % if inputPlaneSmapling is even, the origin is the bottom right of the four central pixels

inputPlane = (x.^2+y.^2).^0.5; 
inputPlane = inputPlane * inputPlaneSamplingStep; % scaling to obtain a map of radial physical distances

inputPlane = intensity/distanceLedInputPlane^2 * cos(atan(inputPlane/distanceLedInputPlane)).^3; % [W.m-2] irradiance map
inputPlane(inputPlane>q) = 0; % circular pupil

% save

% as .mat file
filename = "irradianceMapDistance_" + string(distanceLedInputPlane) + "m_" + "Angle_" + string(theta) + "degree";
save(dirc+filename+".mat", "inputPlane")

% as fits file
fitswrite(inputPlane, dirc+filename+".fits")

% sanity check (does not work)

% si theta = 90°, la somme des irradiances devraient donner 2pi W pour une
% intensitée de 1

total = sum(sum(inputPlane))*inputPlaneSamplingStep^2 / (2*pi); 



