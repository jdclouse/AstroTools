%% HW 5 Master Script

%% Initialize
if ispc
    addpath('C:\Users\John\Documents\Astro\ASEN5070_SOD\tools')
end
clear all
clc

% Cell array to track what functions are used, so they can be published
% later
global function_list;
function_list = {};

% publishing options
pub_opt.format = 'pdf';
pub_opt.outputDir = '.\html';
pub_opt.imageFormat = 'jpeg';
pub_opt.figureSnapMethod = 'entireGUIWindow';
pub_opt.useNewFigure = true ;
pub_opt.maxHeight = Inf;
pub_opt.maxWidth = Inf;
pub_opt.showCode = true;
pub_opt.evalCode = true;
pub_opt.catchError = true;
pub_opt.createThumbnail = true;
pub_opt.maxOutputLines = Inf;

%% Run Problem scripts and publish them

% Problem 1
publish('HW5_P1', pub_opt);

% Problem 2
publish('HW5_P2', pub_opt);

% Problem 3
publish('HW5_P3', pub_opt);


%% Publishing tools and support code
pub_opt.outputDir = '.\tools';
pub_opt.evalCode = false;

%Publish all used functions
function_list = ...
    [function_list; 'C:\Users\John\Documents\Astro\ASEN5070_SOD\tools\fcnPrintQueue'];
for idx = 1:length(function_list)
    publish(function_list{idx}, pub_opt);
end