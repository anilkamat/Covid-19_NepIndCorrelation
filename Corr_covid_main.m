%% setups
clc; clear all; close all
addpath 'D:\Covid_19_Nepal_research\DataSet';
[Data_main, text, alldata] = xlsread('owid-covid-data.xlsx');
%% extract countries
i = 1; j =1;
len_text = size(text,1);
temp = text(:,3);
countries = [];
countries = temp(2);
while(i < len_text)
    if(~isequal(temp(i),temp(i+1)))
        countries{j+1} = temp(i+1);
        j = j+1;
    end
    i = i+1;
end
%%
row_nep = [];
row_ind = [];
i = 1; j = 1; k = 1;
Data_Nep = [];
Data_Ind = [];
while(i< len_text)
    i
    if(isequal(string(text(i,3)),'Nepal'))
        Data_Nep(j,:) = Data_main(i,:);
        j = j+1;
    elseif(isequal(string(text(i,3)),'India'))
        Data_Ind(k,:) = Data_main(i,:);
        k = k+1;
    end
    i = i+1;
end
Data_Nep = Data_Nep(1:size(Data_Nep,1)-1,:)  ; % -1 to remove the incomplete data
Data_Ind = Data_Ind(1:size(Data_Ind,1)-1,:)  ; % -1 to remove the incomplete data
%% Coorelation analysis
[h1,pvalue1] = adftest(Data_Nep(:,2))
[h2, pvalue2] = adftest(Data_Ind(:,2))
% scatter plot
len_min = min(size(Data_Nep,1), size(Data_Ind,1));
figure()
scatter(Data_Nep(1:len_min,2),Data_Ind(1:len_min,2),'+','*');

% Pearson correlation
[R,P,RL,RU] = corrcoef(Data_Nep(1:len_min,2),Data_Ind(1:len_min,2))
t = R(1,2)*sqrt((len_min-2)/(1-R(1,2)^2));
p1= 1-tcdf(t,(len_min-2));

% distance correlation --> non linear correlation
dcor = dist_correlation(Data_Nep(1:len_min,2),Data_Ind(1:len_min,2))
t = dcor*sqrt((len_min-2)/(1-dcor^2));
p1= 1-tcdf(t,(len_min-2));
%% finding lag 
% lag using distance corr --> no need 
% dcor = []; R = [];
% lags = 1:len_min-1;
% for i = 1:len_min-1
%      R = corrcoef(Data_Nep(lags(i):len_min,2),Data_Ind(1:len_min-lags(i)+1,2),'Rows','pairwise');
%      dcor(i) = R;
% end
% figure(12)
% stem(lags,dcor)

% lag using similarity measurement dot product
[c,lags] = xcorr(Data_Nep(1:len_min,2),Data_Ind(1:len_min,2),'normalized');
[~,col] = max(c);
col = lags(col);
M = c.^1;
% second maximum lag for second largest peak
% cc = c(1:400,:);
% [~,col3] = max(cc);
% col3 = lags(col3);
% Define the vertices: the points at (x, f(x)) and (x, 0)
N = length(lags);
verts = [lags(:), c(:); lags(:) zeros(N,1)];
% Define the faces to connect each adjacent f(x) and the corresponding points at y = 0.
q = (1:N-1)';
faces = [q, q+1, q+N+1, q+N];
patch('Faces', faces, 'Vertices', verts, 'FaceVertexCData', [M(:); M(:)], 'FaceColor', 'interp', 'EdgeColor', 'interp'); hold on;
plot(lags,c, 'linewidth',2, 'color',[0 0.5 0])
ylim([0 1])
xline(col,'-.',['Max Lag: ',num2str(col)],'color','black','linewidth',0.5);
xlabel('Lags in days')
ylabel('Correlation')
h = colorbar;
%set(h, 'ylim', [0 1])
title('Cross-Correlations')
hold off;

% find delay between two signal
t21 = finddelay(Data_Nep(:,2),Data_Ind(:,2));
%% plot of incidence rate
figure()
p = plot(Data_Nep(1:len_min,2),'linewidth',2,'Marker','o','Markersize',1); hold on;
plot(Data_Ind(1:len_min,2),'linewidth',2,'Marker','*','Markersize',1);
xlabel('Day')
ylabel('Recorded new cases')
legend('Nepal','India','Location','northwest')
set(gca, 'YScale', 'log')
% modified jet-colormap
cd = [uint8(jet(len_min)*255) uint8(ones(len_min,1))].';
drawnow
set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd)
%% Change point detection
set(groot,'defaultLineLineWidth',1.5)
cpt_Nep = findchangepts(Data_Nep(1:len_min,2),'Statistic','linear','MinThreshold',10000000);
figure()
findchangepts(Data_Nep(1:len_min,2),'Statistic','linear','MinThreshold',10000000);
ylabel('Incidence')
xlabel('Days')
title('Phase determination: Nepal')

cpt_Ind = findchangepts(Data_Ind(1:len_min,2),'Statistic','linear','MinThreshold',10000000000);
figure
findchangepts(Data_Ind(1:len_min,2),'Statistic','linear','MinThreshold',10000000000);
ylabel('Incidence')
xlabel('Days')
title('Phase determination:India')
%% Correlation in each wave
% first wave
X1_Nep =  Data_Nep(cpt_Nep(1):cpt_Nep(3),2);          % cases in first wave
X1_Ind = Data_Ind(cpt_Ind(1):cpt_Ind(3),2);
minLenX1 = min(size(X1_Nep,1), size(X1_Ind,1));
% distance correlation --> non linear correlation
dcor1 = dist_correlation(X1_Nep(1:minLenX1),X1_Ind(1:minLenX1));
t = dcor1*sqrt((minLenX1-2)/(1-dcor1^2));
p_wave1= 1-tcdf(t,(minLenX1-2));

% lag using similarity measurement dot product
[c1,lags1] = xcorr(X1_Nep(1:minLenX1),X1_Ind(1:minLenX1),'normalized');
[~,col1] = max(c1);
col1 = lags1(col1);
M1 = c1.^1;
% % second maximum lag for second largest peak
% cc = c1(1:400,:);
% [~,col3] = max(cc);
% col3 = lags(col3);
% Define the vertices: the points at (x, f(x)) and (x, 0)
N = length(lags1);
verts = [lags1(:), c1(:); lags1(:) zeros(N,1)];
% Define the faces to connect each adjacent f(x) and the corresponding points at y = 0.
q = (1:N-1)';
faces = [q, q+1, q+N+1, q+N];
patch('Faces', faces, 'Vertices', verts, 'FaceVertexCData', [M1(:); M1(:)], 'FaceColor', 'interp', 'EdgeColor', 'interp'); hold on;
plot(lags1,c1, 'linewidth',2, 'color',[0 0.5 0])
ylim([0 1])
xline(col1,'-.',['Max Lag: ',num2str(col1)],'color','black','linewidth',0.5);
xlabel('Lags in days')
ylabel('Correlation')
h = colorbar;
title('Cross-Correlations')
hold off;


% second wave
X2_Nep =  Data_Nep(cpt_Nep(4):end,2);          % cases in first wave
X2_Ind = Data_Ind(cpt_Ind(4):cpt_Ind(end),2);
minLenX2 = min(size(X2_Nep,1), size(X2_Ind,1));
% distance correlation --> non linear correlation
dcor2 = dist_correlation(X2_Nep(1:minLenX2),X2_Ind(1:minLenX2));
t = dcor2*sqrt((minLenX2-2)/(1-dcor2^2));
p_wave2= 1-tcdf(t,(minLenX2-2));

% lag using similarity measurement dot product
[c2,lags2] = xcorr(X2_Nep(1:minLenX2),X2_Ind(1:minLenX2),'normalized');
[~,col2] = max(c2);
col2 = lags2(col2);
M2 = c2.^1;
N = length(lags2);
verts = [lags2(:), c2(:); lags2(:) zeros(N,1)];
% Define the faces to connect each adjacent f(x) and the corresponding points at y = 0.
q = (1:N-1)';
faces = [q, q+1, q+N+1, q+N];
patch('Faces', faces, 'Vertices', verts, 'FaceVertexCData', [M2(:); M2(:)], 'FaceColor', 'interp', 'EdgeColor', 'interp'); hold on;
plot(lags2,c2, 'linewidth',2, 'color',[0 0.5 0])
ylim([0 1])
xline(col2,'-.',['Max Lag: ',num2str(col2)],'color','black','linewidth',0.5);
xlabel('Lags in days')
ylabel('Correlation')
h = colorbar;
title('Cross-Correlations')
hold off;

%% ######################### auto correlation 
[c,lags] = autocorr(Data_Nep(1:len_min,2),'NumLags',len_min-1);  % Nepal
figure()
stem(lags,c)
[c,lags] = autocorr(Data_Ind(1:len_min,2),'NumLags',len_min-1);  % India
figure()
stem(lags,c)

%% whitening the signal
auto_corr = autocorr(Data_Ind(1:len_min,2),'NumLags',200,'NumSTD',2);
figure()
stem(auto_corr)
title('autocorr-without prewhitening')

[mappedX, W, mu_X,X_preWhiten] = prewhiten(Data_Ind(1:len_min,2));

auto_corr = autocorr(X_preWhiten,'NumLags',200,'NumSTD',2);
figure()
stem(auto_corr)
title('autocorr-with prewhitening')

[XnewInd] = SpectralWhitening(Data_Ind(1:len_min,2), 1.1574e-05);  % sampling rate = 1/24/3600 = 1/86400
[XnewNep] = SpectralWhitening(Data_Nep(1:len_min,2), 1.1574e-05);  % sampling rate = 1/24/3600 = 1/86400

auto_corr = autocorr(X_preWhiten);
figure()
stem(auto_corr)

% time serie X
 F=abs(fft(Xnew)); 
 figure(10)
 plot(Xnew(1:end/2));    % shape?
% end of whitening the signal
%%

    