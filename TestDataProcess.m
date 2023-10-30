% This script is used to process test data and calculate estimation errors.
% Jingxuan Chen, 2023.10.30
clear
close all

%% Initialization
fileindex=1;
load('bandpassFIR.mat')
fs=2.5e9;
lambda=299792458/433e6;
%Frith transfer formula
Gt=0;%dBi
Gr=0;%dBi
Loss=0;%dB
Frith_L=@(d,lambda) 20*log10(4*pi*d/lambda)-Gr-Gt+Loss;

%% Load the data file and analyze parameters
loadfilename="C_2_5_5_25_5_-5_17_3.mat";
load("./data/FourArrayAndTwoSource/"+loadfilename)
temp=sscanf(loadfilename,"%c_%d_%d_%d_%d_%d_%d_%d_%d.mat");
Q=temp(2);
true_power=zeros(Q,1);
posx=zeros(Q,1);
posy=zeros(Q,1);
true_theta=zeros(Q,1);
source_distance=zeros(Q,1);

if temp(1)=='P'%polar coordinates
    for q=1:Q
        true_power(q)=temp(3+3*(q-1));
        true_theta(q)=temp(4+3*(q-1));
        source_distance(q)=temp(5+3*(q-1));
    end
elseif temp(1)=='C'%Cartesian coordinates
    for q=1:Q
        true_power(q)=temp(3+3*(q-1));
        posx(q)=temp(4+3*(q-1));
        posy(q)=temp(5+3*(q-1));
    end
    posx=posx*0.6;
    posy=posy*0.6;
    true_theta=atand(posy./posx);
    true_theta(true_theta<0)=true_theta(true_theta<0)+180;
    source_distance=sqrt(posx.^2+posy.^2);
end
[N,K,L]=size(Y_all);
est_DOA=zeros(L,Q);
est_RSS=zeros(L,Q);
est_source_power=zeros(L,Q);
est_Ant_Amp=zeros(N,L);
thetagrid=0:0.1:180;

[~,index]=sort(true_theta);
true_RSS=sqrt(2*50*10.^((((10.*log10(true_power)+300)-300)-Frith_L(source_distance(index),lambda))/10));
%% Analyze each snapshot one by one
for l=1:L
    Y=Y_all(:,:,l);
    %% Draw spectrum of each channel
        f = fs*(0:(K/2))/K;
        figure(1)
        subplot 411
        Y_fft = abs(fft(Y(1,:))/K);
        Y_fft = Y_fft(1:floor(K/2+1));
        Y_fft(2:end-1) = 2*Y_fft(2:end-1);
        Y_fft=Y_fft./2;% Due to the Hilbert transformation
        plot(f./1e6,abs(Y_fft))
        xlabel("Frequency(MHz)")
        ylabel("Amplitude(V)")
        grid on
        xlim([-fs/20/1e6,fs/2/1e6])
        title("Channel 1")

        subplot 412
        Y_fft = abs(fft(Y(2,:))/K);
        Y_fft = Y_fft(1:floor(K/2+1));
        Y_fft(2:end-1) = 2*Y_fft(2:end-1);
        Y_fft=Y_fft./2;% Due to the Hilbert transformation
        plot(f./1e6,abs(Y_fft))
        xlabel("Frequency(MHz)")
        ylabel("Amplitude(V)")
        grid on
        xlim([-fs/20/1e6,fs/2/1e6])
        title("Channel 2")

        subplot 413
        Y_fft = abs(fft(Y(3,:))/K);
        Y_fft = Y_fft(1:floor(K/2+1));
        Y_fft(2:end-1) = 2*Y_fft(2:end-1);
        Y_fft=Y_fft./2;% Due to the Hilbert transformation
        plot(f./1e6,abs(Y_fft))
        xlabel("Frequency(MHz)")
        ylabel("Amplitude(V)")
        grid on
        xlim([-fs/20/1e6,fs/2/1e6])
        title("Channel 3")

        subplot 414
        Y_fft = abs(fft(Y(4,:))/K);
        Y_fft = Y_fft(1:floor(K/2+1));
        Y_fft(2:end-1) = 2*Y_fft(2:end-1);
        Y_fft=Y_fft./2;% Due to the Hilbert transformation
        plot(f./1e6,abs(Y_fft))
        xlabel("Frequency(MHz)")
        ylabel("Amplitude(V)")
        grid on
        xlim([-fs/20/1e6,fs/2/1e6])
        title("Channel 4")

        %% Filtering and truncation
        filterK=round(K*0.8);
        filteredwave=filter(bandpassFIR,Y,2);
        cutofftime=K-filterK+1;
        filterY=filteredwave(:,cutofftime:end);
        Y=filterY;

        %% DOA/RSS joint estimation with MUSIC
        [detectedtheta,spectrum,amplitudeS]=MUSIC_Amp(Y,Q,thetagrid);
        est_DOA(l,:)=sort(detectedtheta,'ascend');
        est_RSS(l,:)=amplitudeS(find(sum(thetagrid==detectedtheta.',1)));
        disp("Snapshot index: "+string(l))
        disp("Estimate DOA: "+string(est_DOA(l,:).')+"    True DOA: "+string(sort(true_theta)))
        disp("Estimate RSS: "+string(est_RSS(l,:).'*1e3)+" mV    True RSS: "+string(true_RSS*1e3)+" mV")
        figure(2)
        h=plot(thetagrid,amplitudeS,'Linewidth',2);
        hold on
        for q=1:Q
            plot(true_theta(q).*ones(Q,2),[0,0.2],'--k')
        end
        hold off
        xlabel('DOA (degrees)');
        ylabel('RSS (V)');
        xlim([0,180])
        ylim([0,0.5])
        grid on;
        
        figure(3)
        plot(thetagrid,spectrum,'Linewidth',2);
        hold on
        for q=1:Q
            plot(true_theta(q).*ones(Q,2),[0,-20],'--k')
        end
        hold off
        xlabel('DOA (degrees)');
        ylabel('Normalized spatial spectrum (dB)');
        xlim([0,180])
        ylim([-20,0])
        grid on
        
        % pause
end

estDOA_mean=mean(est_DOA,1);
estDOA_median=median(est_DOA,1);
estDOA_std=std(est_DOA,1);
estRSS_mean=mean(est_RSS,1);
estRSS_median=median(est_RSS,1);
estRSS_std=std(est_RSS,1);
mu = 1.4796;%Loss coefficient for MUSIC
calibrated_est_source_RSS=mu.*est_RSS;

disp("DOA median: "+string(estDOA_median.')+"    std: "+string(estDOA_std.'))
disp("RSS median: "+string(estRSS_median.'*1e3)+" mV    std: "+string(estRSS_std.'*1e3)+" mV")

figure(4)
plot(sort(true_theta),'^','MarkerSize',10)
hold on
boxplot(est_DOA);
grid on
xlabel("Source Index")
ylabel("DOA estimated (degrees)")

figure(5)
plot(mag2db(true_RSS)+30,'^','MarkerSize',10)
hold on
boxplot(mag2db(calibrated_est_source_RSS)+30);
grid on
xlabel("Source Index")
ylabel("RSS estimated (dBmV)")

total_DOA_RMSE=(sum((sort(true_theta)-sort(estDOA_median.')).^2,1)./Q).^0.5;
mu = 1.4796;%Loss coefficient for MUSIC
calibrated_RSS=mu.*estRSS_median;
total_RSS_RMSE=(sum((sort(true_RSS)-sort(calibrated_RSS.')).^2,1)./Q).^0.5;

disp("RMSE of DOA: "+string(total_DOA_RMSE)+" degree");
disp("RMSE of RSS: "+string(total_RSS_RMSE*1e3)+" mV");