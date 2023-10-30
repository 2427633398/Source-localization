% This script is used to process test data to calibrate the array pattern.
% Jingxuan Chen, 2023.10.30
clear
close all

%% Initialization
load('bandpassFIR.mat')
fs=2.5e9;
Q=1;
lambda=299792458/433e6;
%Frith transfer formula
Gt=0;%dBi
Gr=0;%dBi
Loss=0;%dB
Frith_L=@(d,lambda) 20*log10(4*pi*d/lambda)-Gr-Gt+Loss;

filelist=dir("./data/PatternCalibration/");
filenum=length(filelist)-2;
est_source_power=zeros(100,filenum);
est_source_DOA=zeros(100,filenum);
est_source_RSS=zeros(100,filenum);
source_theta_vec=zeros(1,filenum);
true_RSS=zeros(1,filenum);
%% Load all data files in the folder
for fileindex=1:filenum
    loadfilename=filelist(fileindex+2).name;
    disp(loadfilename)
    load("./data/PatternCalibration/"+loadfilename)
    %Read data file
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
    est_Ant_Amp=zeros(N,L);
    true_RSS(fileindex)=sqrt(2*50*10.^((((10.*log10(true_power)+300)-300)-Frith_L(source_distance,lambda))/10));
    thetagrid=0:0.1:180;

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
        est_DOA(l,:)=detectedtheta;
        est_RSS(l,:)=max(amplitudeS);
        disp("Estimate DOA: "+string(detectedtheta)+"    True DOA: "+string(true_theta))
        disp("Estimate RSS: "+string(est_RSS(l,:)*1e3)+" mV    True RSS: "+string(true_RSS(fileindex)*1e3)+" mV")
        figure(2)
        h=plot(thetagrid,amplitudeS,'Linewidth',2);
        hold on
        plot(true_theta*ones(Q,2),[0,0.05])
        hold off
        xlabel('DOA (degrees)');
        ylabel('RSS (V)');
        xlim([0,180])
        ylim([0,0.5])
        grid on;
        
        figure(3)
        plot(thetagrid,spectrum,'Linewidth',2);
        hold on
        plot(true_theta*ones(Q,2),[0,-20])
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

    est_source_DOA(:,fileindex)=est_DOA;
    fprintf("DOA median:%.2f, std:%.2f\n",estDOA_median,estDOA_std)
    est_source_RSS(:,fileindex)=est_RSS;
    fprintf("RSS median:%.2f mV, std:%.2f mV\n",estRSS_median*1e3,estRSS_std*1e3)
    source_theta_vec(fileindex)=true_theta;
end
%% Sort results by distance
[true_source_theta,tempindex]=sort(source_theta_vec);
est_source_DOA=est_source_DOA(:,tempindex);
est_source_RSS=est_source_RSS(:,tempindex);
mu = 1.4796;%Loss coefficient for MUSIC
calibrated_est_source_RSS=mu.*est_source_RSS;
true_RSS=true_RSS(:,tempindex);
median(est_source_DOA)-true_source_theta.';
%% draw
figure(4)
boxplot(est_source_DOA);
xticklabels(true_source_theta)
xlabel("True DOA (degrees)")
ylabel("DOA estimated (degrees)")
hold on 
grid on
ylim([20,100])
yticks(20:10:100)
set(gca,'FontName','Times New Roman','FontSize',22,'LineWidth',1);

figure(5)
plot(mag2db(true_RSS)+30,'-k','LineWidth',1)
hold on
boxplot(mag2db(calibrated_est_source_RSS)+30);
xticklabels(true_source_theta)
xlabel("True DOA (degrees)")
ylabel("RSS estimated (dBmV)")
yticks(4:17)
grid on
set(gca,'FontName','Times New Roman','FontSize',22,'LineWidth',1);