clear
close all
clc
addpath('./functions')
%% basic parameters
path2data='./dataGIT/';
frameRate=1;
%% data structure
field1='measurementID';field2='timeRange';field3='startValues';
%
values1={'mp1_600','mp1_1000','mp1_1500','mp1_1000_posB','mp1_ref', ...
         'mp2_600','mp2_1000','mp2_1500','mp2_1000_posB','mp2_ref', ...
         'mp3_600','mp3_1000','mp3_1500','mp3_1000_posB','mp3_ref','mp3_1000_noFilter', ...
         'mp4_600','mp4_1000','mp4_1500','mp4_1000_posB','mp4_ref', ...
         'mp5_600','mp5_1000','mp5_1500','mp5_1000_posB','mp5_ref', ...
         'mp6_600','mp6_1000','mp6_1500','mp6_1000_posB','mp6_ref', ...
         'mp3_1000_classroom_noPeople','mp3_1200_classroom_wPeople', ...
         'mp1_600_hallway','mp1_1000_hallway','mp1_1585_hallway','mp1_reference_hallway',...
         'mp2_600_hallway','mp2_1000_hallway','mp2_1585_hallway','mp2_reference_hallway'};
% time range for exponential fit
values2={[.4 .75],[0.2 0.8],[0.15 0.6],[0.6 1.0],[0.3 1.1], ...
         [.25 .8],[0.25 0.8],[0.2 0.6],[0.6 1.0],[0.3 1.1], ...
         [.25 .8],[0.2 0.7],[0.25 0.6],[0.6 1.0],[0.3 1.1],[0.15 0.5], ...
         [.3 .8],[0.15 0.7],[0.2 0.6],[0.6 1.0],[0.3 1.1], ...
         [.45 .85],[0.3 0.9],[0.1 0.6],[0.6 1.0],[0.3 1.1], ...
         [.25 .8],[0.3 0.9],[0.2 0.6],[0.6 1.0],[0 .5], ...
         [0.05 0.5],[0.05 0.5], ...
         [0.4 1.0],[0.25 0.9],[0.3 0.7],[0.21 0.5], ...
         [0.2 0.75],[0.4 0.9],[0.4 0.8],[0.2 0.5]};
%
values3={[10000 3],[10000 3],[10000 3],[10000 3],[10000 .2], ...
         [10000 3],[10000 3],[10000 3],[10000 3],[10000 .2], ...
         [10000 3],[10000 3],[10000 3],[10000 3],[10000 .2],[10000 .2], ...
         [10000 3],[10000 3],[10000 3],[10000 3],[10000 .2], ...
         [10000 3],[10000 3],[10000 3],[10000 3],[10000 .2], ...
         [10000 3],[10000 3],[10000 3],[10000 3],[10000 .2], ...
         [10000 3],[10000 3], ...
         [10000 3],[10000 3],[10000 3],[10000 3], ...
         [10000 3],[10000 3],[10000 3],[10000 3]};
%
paramStruct=struct(field1,values1,field2,values2,field3,values3);
%%
expCoefficients=cell(size(paramStruct,2),2);
for k1=1:size(paramStruct,2)
    
    %% particle image count, seeded air:
    path2info=[path2data paramStruct(k1).measurementID '.csv'];
    davisTable=readtable(path2info,'HeaderLines',1);
    noParticles=double(davisTable.Var2)';
    % noParticlesMean=movmean(noParticles,15);
    timeInfo=0:(1/(frameRate*3600)):size(noParticles,2)/(frameRate*3600)-(1/(frameRate*3600));
    %%
    ind=find(timeInfo>paramStruct(k1).timeRange(1) & timeInfo<paramStruct(k1).timeRange(2));
    [fitParam,~]=fitExpDecay(timeInfo(ind), noParticles(ind), paramStruct(k1).startValues, paramStruct(k1).measurementID);
    expCoefficients{k1,2}=[fitParam.a fitParam.b log(.5)/(-fitParam.b) 1/fitParam.b];
    expCoefficients{k1,1}=paramStruct(k1).measurementID;
       
end
%%
function [fitresult, gof] = fitExpDecay(timeInfo, noParticles, startValues, measurementID)

    [xData, yData] = prepareCurveData( timeInfo, noParticles );
    %%
    ft = fittype( 'a*exp(-b*x)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Algorithm = 'Levenberg-Marquardt';
    opts.Display = 'Off';
    opts.Robust = 'LAR';
    opts.StartPoint = [startValues(1) startValues(2)];

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    
    % Plot fit with data.
    plotOn=0;
    if plotOn==1
        figure( 'Name', measurementID );
        h = plot( fitresult, xData, yData);
        set(h,'LineWidth',3)
        set(gca,'yscale','log')
        legend( h, 'particle image count over time', 'fit of exp. decay', 'Location', 'NorthEast' );
        % Label axes
        xlabel timeInfo
        ylabel noParticles
        grid on
    end
    
end
