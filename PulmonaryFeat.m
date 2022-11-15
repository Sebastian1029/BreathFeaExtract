% =========================================================================
% Features extraction
% Authors: Sebastian Escobar-Pajoy and J.P. Ugarte
% Version: 1.0
% 2022
% =========================================================================
% Description
% =========================================================================
% This script implemements the feature extraction scheme used in 
% "S. Escobar-Pajoy and J.P. Ugarte, Computerized analysis of pulmonary 
% sounds using uniform manifold projection (2022). The implementation are
% cross-referenced with the corresponding manuscript section number.
% =========================================================================
% Requirements
% =========================================================================
% The Matlab Audio Toolbox is required.
% =========================================================================
% =========================================================================

% Load the signal T in the workspace. The working sample rate is 4 KHz
%% section 3.1
% Fourier Transform.
x=(T-mean(T))./std(T); %Signal normalization
X=fft(x);
%% section 3.2 and 3.3
% 
% Three types of features are extracted from signal T: spectral features
% Mel frequency cepstral coefficients and chroma features.

hs=113;% Overlapped samples
L=125;% Window size
x=(T-mean(T))./std(T);
[s1, fs1]=spectrogram(x,hanning(L),hs,L,4000);

THETA=zeros(size(s1,2),36); %Features matrix
% =========================================================================
% Mel Frequency Cepstral Coefficients 
[Mfcc]=mfcc(s1,4000,'LogEnergy','ignore'); 
% =========================================================================
% Chroma vector
Chroma=ChromaVector(abs(s1).^2).'; 
% =========================================================================
THETA(:,1:13)=Mfcc;
THETA(:,14:25)=Chroma;
% =========================================================================
% The spectral features are described at 
% https://la.mathworks.com/help/audio/ug/spectral-descriptors.html#mw_rtc_SpectralDescriptorsExample_7E5C532B
% and in the book  "An Introduction to Audio Content Analysis 
% Applicationsin Signal Processing and Music Informatics" by Alexander Lerch.
% These features are extracted from the three following frequency bands: 
% low (60-300Hz), mid (300-800Hz) and high (800-1800Hz). 

% (In this example, the mid-band is characterized)
Fi=300;% Band lower frequency (60 low, 300 mid, 800 high)
Ff=800;% Band upper frequency (300 low, 800 mid, 1800 high)
Bi=round(Fi*L/4000);
Bf=round(Ff*L/4000);
s=s1(Bi:Bf,:);
fs=fs1(Bi:Bf);

%Spectral features
THETA(:,26)=spectralEntropy(abs(s).^2,fs);
[FSkewness,FSpread,FCentroid]=spectralSkewness(abs(s).^2,fs);
THETA(:,27)=FSkewness;
THETA(:,28)=FCentroid;
THETA(:,29)=FSpread;
THETA(:,30)=spectralSlope(abs(s).^2,fs);
THETA(:,31)=spectralCrest(abs(s).^2,fs);
THETA(:,32)=spectralFlux(abs(s).^2,fs);
THETA(:,33)=spectralKurtosis(abs(s).^2,fs);
THETA(:,34)=spectralFlatness(abs(s).^2,fs);
THETA(:,35)=spectralRolloffPoint(abs(s).^2,fs);
THETA(:,36)=spectralDecrease(abs(s).^2,fs);
%% =========================================================================
% Descriptive statistics, e.g., the section 3.3
Mean=mean(THETA);
StandardDeviation=std(THETA);
CoefVariance=StandardDeviation./Mean;
% =========================================================================
% Fractional state space reprensentation. The Grunwald–Letnikov definition 
% of fractional order derivative is implemented.
% (In this example we calculate the 0.7 derivative of one of the feauture)
% signals
h=12/4000; %Time diferential
alpha=0.7; % Fractional order, can take values between -1 and 2
n=100; %Size of the buffer for the derivative calculation
D=FracDer(THETA(:,19),h,alpha,n);
% =========================================================================

% =========================================================================
% Auxiliary functions
% =========================================================================
function [Chromagram] = ChromaVector(s)
    BW=0.5;
    FrecPrin=[61.41,69.3,73.42,77.78,82.41,87.31,92.50,98,103.83,110,116.54,123.47];
    LowEdge=2^(-BW/12).*FrecPrin;
    HighEdge=2^(BW/12).*FrecPrin;
    Num=size(s,2);
    Chromagram=zeros(12,Num);
    for j=1:Num   
        X=s(:,j); 
        for gg=1:length(FrecPrin)
            for oct=1:5
            BinLow=floor((length(X)./4000)*LowEdge(gg)*2^(oct-1));
            BinHigh=floor((length(X)./4000)*HighEdge(gg)*2^(oct-1));
            if BinHigh>=2000
                BinHigh=1999;
            end
            Chromagram(gg,j)=(mean(X(1+(BinLow:BinHigh))))+Chromagram(gg,j);
            end            
        end
        Chromagram(:,j)=(Chromagram(:,j))./sum(Chromagram(:,j));   
    end
end


function [D] = FracDer(x,h,alpha,n)
h_a = h^alpha;
N=length(x);
D=zeros(100,1);
for i = 1 : N
    sigma = 0;
    cm=1;
    for m = 0 : 1 : n
    if i - m < 1
        break;
    end 
        sigma = sigma + cm* x(i - m);
        cm=cm*(1-(1+alpha)/(m+1));
    end
    D(i) = sigma / h_a;
end

end