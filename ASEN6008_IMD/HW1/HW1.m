%% HW 2 Master Script

%% Initialize
clear all
clc
hw_pub.toolsPath = @(x) ...
    strcat('C:\Users\John\Documents\Astro\ASEN5050\tools\',x);
if ispc
    addpath(hw_pub.toolsPath(''))
end

hw_pub.script_list = {};
hw_pub.script_list = [hw_pub.script_list 'HW1_P1'];
hw_pub.script_list = [hw_pub.script_list 'HW1_P2'];
% Cell array to track what functions are used, so they can be published
% later
global function_list;
function_list = {};

% publishing options
hw_pub.pub_opt.format = 'pdf';
hw_pub.pub_opt.outputDir = '.\pdf';
hw_pub.pub_opt.imageFormat = 'bmp';
hw_pub.pub_opt.figureSnapMethod = 'entireGUIWindow';
hw_pub.pub_opt.useNewFigure = true ;
hw_pub.pub_opt.maxHeight = Inf;
hw_pub.pub_opt.maxWidth = Inf;
hw_pub.pub_opt.showCode = true;
hw_pub.pub_opt.evalCode = true;
hw_pub.pub_opt.catchError = true;
hw_pub.pub_opt.createThumbnail = true;
hw_pub.pub_opt.maxOutputLines = Inf;

% helpful figure dimensions for better resolution in the pdf
hw_pub.figWidth = 1120; % pixels
hw_pub.figHeight = 840; % pixels
hw_pub.figPosn = [0, 0, hw_pub.figWidth, hw_pub.figHeight];
% Example: some_fig = figure('Position', hw_pub.figPosnfigure);
hw_pub.lineWidth = 2; % pixels
% Example: plot(x,y,'r--','LineWidth',hw_pub.lineWidth);
% Example: h1 = plot(x,y)
%          h2 = plot(x,z, 'r');
%          h3 = plot(x,a, 'k');
%          legend([h1 h3],{'hello','world'});

%% Run Problem scripts and publish them
for idx = 1:length(hw_pub.script_list)
    fprintf(['Publishing ' hw_pub.script_list{idx} '\n'])
    publish(hw_pub.script_list{idx}, hw_pub.pub_opt);
end

%% Publishing tools and support code
hw_pub.pub_opt.outputDir = '.\pdf\tools';
hw_pub.pub_opt.evalCode = false;

%Publish all used functions
if isempty(function_list)
    return % nothing to publish
end
function_list = ...
    [function_list hw_pub.toolsPath('fcnPrintQueue.m')];
for idx = 1:length(function_list)
    publish(function_list{idx}, hw_pub.pub_opt);
end

pdfCmdStr = 'pdftk ';
for idx = 1:length(hw_pub.script_list)
    pdfCmdStr = [pdfCmdStr 'pdf\' hw_pub.script_list{idx} '.pdf ']; %#ok<*AGROW>
end
for idx = 1:length(function_list)
    [fpath, filename, ext] = fileparts(function_list{idx});
%     tmp = strsplit(function_list{idx},'\');
%     pdfCmdStr = [pdfCmdStr 'pdf\tools\' tmp{end}(1:end-2) '.pdf '];
    pdfCmdStr = [pdfCmdStr 'pdf\tools\' filename '.pdf '];
end

[fpath, hw_num] = fileparts(pwd);
% tmp = strsplit(pwd,'\');
% hw_num = tmp{end};
pdfCmdStr = [pdfCmdStr 'output pdf\Clouse_ASEN6008_' hw_num '.pdf'];
if system(pdfCmdStr);
    fprintf('Failed to PDF!\n')
end