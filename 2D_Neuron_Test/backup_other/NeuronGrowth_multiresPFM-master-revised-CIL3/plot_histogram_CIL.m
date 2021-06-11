%% This code generates plot of histogram for the length and angle information
clc
clear all
close all

load 'lenm_CIL';
load 'segm_CIL';
load 'angm_CIL';

%% Plot total histogram
Y = lenm_CIL;
numOfBins = 10;
[histFreq, histXout] = hist(Y, numOfBins);

figure;
bar(histXout, histFreq/sum(histFreq)*100);
mn=mean(Y);
stdv=std(Y);
mnlabel=sprintf('Mean-- %0.2f', mn);
stdlabel=sprintf('S. D.-- %0.2f', stdv);
xlabel('Total Length');
ylabel('Frequency');
h=annotation('textbox',[0.58 0.75 0.1 0.1]);
set(h,'String',{mnlabel,stdlabel});

%% Plot angle histogram
Y = angm_CIL;
mn=mean(Y);
stdv=std(Y);
mnlabel=sprintf('Mean-- %0.2f', mn);
stdlabel=sprintf('S. D.-- %0.2f', stdv);
numOfBins = 10;
[histFreq, histXout] = hist(Y, numOfBins);

figure;
bar(histXout, histFreq/sum(histFreq)*100);
xlabel('Relative Angle');
ylabel('Frequency');
h=annotation('textbox',[0.58 0.75 0.1 0.1]);
set(h,'String',{mnlabel,stdlabel});


%% Plot segment length histogram
Y = segm_CIL;
mn=mean(Y);
stdv=std(Y);
mnlabel=sprintf('Mean-- %0.2f', mn);
stdlabel=sprintf('S. D.-- %0.2f', stdv);

numOfBins = 10;
[histFreq, histXout] = hist(Y, numOfBins);

figure;
bar(histXout, histFreq/sum(histFreq)*100);
xlabel('Segment Length');
ylabel('Frequency');
h=annotation('textbox',[0.58 0.75 0.1 0.1]);
set(h,'String',{mnlabel,stdlabel});