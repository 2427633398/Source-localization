% Multiple signal classification (MUSIC) with RSS estimation.
function [detectedtheta,spectrum,amplitudeS]=MUSIC_Amp(Y,Q,thetagrid)
% Return detectedtheta indicates the detected angle
% Return spectrum to represent spatial spectrum
% Return amplitudeS for amplitude spectrum
% Parameter Y is the array receiving data matrix, N * K
% Parameter Q is the number of radiation sources
% Parameter H is the guiding matrix
% Jingxuan Chen, 2023.10.30

[N,K]=size(Y);
M=length(thetagrid);
R=Y*Y'/K;

dd = 0.5;
small_d=0:dd:(N-1)*dd;

[EV,D]=eig(R);                 
EVA=diag(D).';                   
[~,I]=sort(EVA);                
EV=fliplr(EV(:,I));               
En=EV(:,Q+1:N);                   

spectrum=zeros(length(thetagrid),1);
for i = 1:M
    a=exp(-1j*2*pi*small_d*cosd((thetagrid(i)))).';
    spectrum(i)=1/(a'*(En*En')*a);
end
spectrum=abs(spectrum)./max(abs(spectrum));
spectrum=10*log10(spectrum);

[pks,locs] = findpeaks(spectrum);
[~,I]=sort(pks,'descend');
try
    detectedtheta=thetagrid(locs(I(1:Q)));
catch
    detectedtheta=thetagrid(locs(I));
end

d=0:0.5:(N-1)*0.5;
A=exp(-1j*2*pi*d'.*cosd(detectedtheta));
amplitude=sqrt(real(diag((A'*A)\A'*R*pinv(A'))));
amplitudeS=zeros(M,1);
amplitudeS(locs(I(1:Q)))=amplitude;
end