%% HW 3 Master Script

%% Initialize
if ispc
    addpath('C:\Users\John\Documents\ASEN5010\tools')
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
pub_opt.imageFormat = 'bmp';
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

% Problem 2: S&J, Problem 3.23
publish('HW3_P2', pub_opt);

% Problem 4: S&J, Problem 3.28
publish('HW3_P4', pub_opt);

% Problem 6: Find DCMs from sensor measurements using the triad method
publish('HW3_P6', pub_opt);

% Problem 7: Find DCMs from sensor measurements using the Davenport method
publish('HW3_P7', pub_opt);

% Problem 8: Find DCMs from sensor measurements using the QUEST method
publish('HW3_P8', pub_opt);

% Problem 9: S&J, Problem 2.12
publish('HW3_P9', pub_opt);

% Problem 12: S&J, Problem 4.3
publish('HW3_P12', pub_opt);

%% Publishing tools and support code
pub_opt.outputDir = '.\tools';
pub_opt.evalCode = false;

%Publish all used functions
function_list = ...
    [function_list; 'C:\Users\John\Documents\ASEN5010\tools\fcnPrintQueue'];
for idx = 1:length(function_list)
    publish(function_list{idx}, pub_opt);
end