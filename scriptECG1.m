
%Analisis de señal ECG,visulaización, filtrado y segmentación.

load('ecg.mat')

Fs = 250; %Frecuencia de muestreo
G = 2000; %Numero de muestras



ecg = ecg/G;
ecg = (ecg - mean(ecg))/std(ecg);
t = (1:1:length(ecg))*(1/Fs);

%Graficamos el ECG original
figure;
plot(t,ecg)
xlim([0 4])
xlabel('Tiempo (s)')
ylabel('Amplitud (mV)')
title('ECG Dominio del Tiempo')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Filtros para eliminar artefactos de desvio de linea
%%EM%%%%%%%%%%%%%%%%%%%%



Fs = 250;  % Frecuencia muestreo

Fpass1 = 59;          % First Passband Frequency
Fstop1 = 59.5;        % First Stopband Frequency
Fstop2 = 60;          % Second Stopband Frequency
Fpass2 = 61;          % Second Passband Frequency
Apass1 = 0.5;         % First Passband Ripple (dB)
Astop  = 3;           % Stopband Attenuation (dB)
Apass2 = 0.5;         % Second Passband Ripple (dB)
match  = 'stopband';  % Band to match exactly


h  = fdesign.bandstop(Fpass1, Fstop1, Fstop2, Fpass2, Apass1, Astop, ...
                      Apass2, Fs);
Hd1 = design(h, 'butter', 'MatchExactly', match);
ecg_limpio1 =filter(Hd1,ecg);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Filtro notch%%%%%%%%%%%%%%



Fnotch = 60; % Frecuencia notch
BW = 120; % Ancho banda
Apass = 1; % Atenuacion de ancho banda
[b, a] = iirnotch (Fnotch/ (Fs/2), BW/ (Fs/2), Apass);
Hd2 = dfilt.df2 (b, a);
ecg_limpio2 =filter(Hd2,ecg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Filtro Chebyshev I%%%%%%%%%%%%%

Fs = 250;  % Sampling Frequency

Fpass1 = 59;          % First Passband Frequency
Fstop1 = 59.5;        % First Stopband Frequency
Fstop2 = 60;          % Second Stopband Frequency
Fpass2 = 61;          % Second Passband Frequency
Apass1 = 0.5;         % First Passband Ripple (dB)
Astop  = 3;           % Stopband Attenuation (dB)
Apass2 = 0.5;         % Second Passband Ripple (dB)
match  = 'passband';  % Band to match exactly


h  = fdesign.bandstop(Fpass1, Fstop1, Fstop2, Fpass2, Apass1, Astop, ...
                      Apass2, Fs);
Hd3 = design(h, 'cheby1', 'MatchExactly', match);
ecg_limpio3 =filter(Hd3,ecg);


%%%%%%%%%%%%%%%%%%%%%%%%%Filtro Chebyshev II%%%%%%%%%%%%%%%%%%%%%%%%

Fs = 250;  % Sampling Frequency

Fpass1 = 59;          % First Passband Frequency
Fstop1 = 59.5;        % First Stopband Frequency
Fstop2 = 60;          % Second Stopband Frequency
Fpass2 = 61;          % Second Passband Frequency
Apass1 = 0.5;         % First Passband Ripple (dB)
Astop  = 3;           % Stopband Attenuation (dB)
Apass2 = 0.5;         % Second Passband Ripple (dB)
match  = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its CHEBY2 method.
h  = fdesign.bandstop(Fpass1, Fstop1, Fstop2, Fpass2, Apass1, Astop, ...
                      Apass2, Fs);
Hd4 = design(h, 'cheby2', 'MatchExactly', match);
ecg_limpio4 =filter(Hd4,ecg);



%%%%%%%%%%%%%%%%%%%%%%%%Filtro Cauer%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Fs = 250;  % Sampling Frequency

Fpass1 = 59;      % First Passband Frequency
Fstop1 = 59.5;    % First Stopband Frequency
Fstop2 = 60;      % Second Stopband Frequency
Fpass2 = 61;      % Second Passband Frequency
Apass1 = 0.5;     % First Passband Ripple (dB)
Astop  = 3;       % Stopband Attenuation (dB)
Apass2 = 0.5;     % Second Passband Ripple (dB)
match  = 'both';  % Band to match exactly

h  = fdesign.bandstop(Fpass1, Fstop1, Fstop2, Fpass2, Apass1, Astop, ...
                      Apass2, Fs);
Hd5 = design(h, 'ellip', 'MatchExactly', match);

ecg_limpio5 =filter(Hd5,ecg);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graficar ECG filtrados
figure;
subplot(3,2,1)
plot(t,ecg_limpio1);
xlim([0 4])
title('ECG Filtrado con Butterworth n=4')
xlabel('Tiempo (s)')
ylabel('Amplitud (mV)')
subplot(3,2,2)
plot(t,ecg_limpio2);
xlim([0 4])
title('ECG Filtrado con Notch')
xlabel('Tiempo (s)')
ylabel('Amplitud (mV)')
subplot(3,2,3)
plot(t,ecg_limpio3);
xlim([0 4])
title('ECG Filtrado con Chebyshev tipo I')
xlabel('Tiempo (s)')
ylabel('Amplitud (mV)')
subplot(3,2,4)
plot(t,ecg_limpio4);
xlim([0 4])
title('ECG Filtrado con fitro Chebyshev tipo II')
xlabel('Tiempo (s)')
ylabel('Amplitud (mV)')
subplot(3,2,5)
plot(t,ecg_limpio5);
xlim([0 4])
title('ECG Filtrado con fitro Eliptico')
xlabel('Tiempo (s)')
ylabel('Amplitud (mV)')












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%Filtros para eliminar linea base%%%%%%%%%%%%%%

%%%%%%%%%%%Filtro notch de orden 3%%%%%%%%%%%%%%%

Fs = 250 % Sampling Frequency
Fnotch = 0.67; % Notch Frequency
BW = 5; % Bandwidth
Apass = 1; % Bandwidth Attenuation
[b, a] = iirnotch (Fnotch/ (Fs/2), BW/(Fs/2), Apass);
Hd1 = dfilt.df2 (b, a);
ecg_limpio1=filter(Hd1,ecg);



%%%%%%%%%%%%%%%%Butterworth pasa bajas orden 41%%%%%%


Fs = 250;  % Sampling Frequency

Fpass = 0.6;        % Passband Frequency
Fstop = 0.667;       % Stopband Frequency
Apass = .5;         % Passband Ripple (dB)
Astop = 20;           % Stopband Attenuation (dB)
match = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its CHEBY1 method.
h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, Fs);
Hd2 = design(h, 'butter', 'MatchExactly', match);




%%%%%%%%%%%%%%%%%%%%%Chebyshev I orden 8%%%%%%%%%%

Fs = 250;  % Sampling Frequency

Fpass = 0.65;        % Passband Frequency
Fstop = 0.667;       % Stopband Frequency
Apass = 0.5;         % Passband Ripple (dB)
Astop = 6;           % Stopband Attenuation (dB)
match = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its CHEBY1 method.
h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, Fs);
Hd3 = design(h, 'cheby1', 'MatchExactly', match);

ecg_limpio3 =filter(Hd3,ecg);



%%%%%%%%%%%%%%%%Filtro Chebyshev II orden 8%%%%%%%%%%%%

Fs = 250;  % Sampling Frequency

Fpass = 0.65;        % Passband Frequency
Fstop = 0.667;       % Stopband Frequency
Apass = 0.5;         % Passband Ripple (dB)
Astop = 6;           % Stopband Attenuation (dB)
match = 'stopband';  % Band to match exactly


h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, Fs);
Hd4 = design(h, 'cheby2', 'MatchExactly', match);

ecg_limpio4 =filter(Hd4,ecg);

%%%%%%%%%%%%%%%%Filtro Eliptico orden 3%%%%%%%%%%%%%



Fs = 250;  % Sampling Frequency

Fpass = 0.65;        % Passband Frequency
Fstop = 0.667;       % Stopband Frequency
Apass = 0.5;         % Passband Ripple (dB)
Astop = 20;           % Stopband Attenuation (dB)
match = 'stopband';  % Band to match exactly


h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, Fs);
Hd5 = design(h, 'ellip', 'MatchExactly', match);
ecg_limpio5 =filter(Hd5,ecg);


%ecg_limpio2=filter(Hd1,ecg);

ecg_limpio3=filter(Hd4,ecg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graficar ECG filtrados
figure;
subplot(3,2,1)
plot(t,ecg_limpio2);
xlim([0 4])
title('ECG Filtrado con Butterworth n=41')
xlabel('Tiempo (s)')
ylabel('Amplitud (mV)')
subplot(3,2,2)
plot(t,ecg_limpio1);
xlim([0 4])
title('ECG Filtrado con Notch')
xlabel('Tiempo (s)')
ylabel('Amplitud (mV)')
subplot(3,2,3)
plot(t,ecg_limpio3);
xlim([0 4])
title('ECG Filtrado con Chebyshev tipo I oden 8')
xlabel('Tiempo (s)')
ylabel('Amplitud (mV)')
subplot(3,2,4)
plot(t,ecg_limpio4);
xlim([0 4])
title('ECG Filtrado con fitro Chebyshev tipo II oden 8')
xlabel('Tiempo (s)')
ylabel('Amplitud (mV)')
subplot(3,2,5)
plot(t,ecg_limpio5);
xlim([0 4])
title('ECG Filtrado con fitro Eliptico orden 3')
xlabel('Tiempo (s)')
ylabel('Amplitud (mV)')



%%%%%%%%%%%%%FILTRADO FINAL TIPO A%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load('ecg.mat')

Fs = 250; %Frecuencia de muestreo
G = 2000; %Numero de muestras



ecg = ecg/G;
ecg = (ecg - mean(ecg))/std(ecg);
t = (1:1:length(ecg))*(1/Fs);

Fnotch = 0.67; % Notch Frequency
BW = 5; % Bandwidth
Apass = 1; % Bandwidth Attenuation
[b, a] = iirnotch (Fnotch/ (Fs/2), BW/(Fs/2), Apass);
Hd1 = dfilt.df2 (b, a);
ecg_limpio1=filter(Hd1,ecg);



Fpass1 = 59;          % First Passband Frequency
Fstop1 = 59.5;        % First Stopband Frequency
Fstop2 = 60;          % Second Stopband Frequency
Fpass2 = 61;          % Second Passband Frequency
Apass1 = 0.5;         % First Passband Ripple (dB)
Astop  = 3;           % Stopband Attenuation (dB)
Apass2 = 0.5;         % Second Passband Ripple (dB)
match  = 'stopband';  % Band to match exactly


h  = fdesign.bandstop(Fpass1, Fstop1, Fstop2, Fpass2, Apass1, Astop, ...
                      Apass2, Fs);
Hd2 = design(h, 'butter', 'MatchExactly', match);
ecg_limpio2 =filter(Hd2,ecg_limpio1);




%Filtro MA para eliminar componenente EM
windowSize = 5; 
b = (1/windowSize)*ones(1,windowSize);
a=1;
ecg_limpio3= filter(b, a, ecg_limpio2);
ECGA=ecg_limpio3;
subplot(2,1,1)
plot(t,ecg)
xlim([0 4])
title('A)')
xlabel('Time (s)')
ylabel('Amplitude (mV)')
subplot(2,1,2)
plot(t,ecg_limpio3)
xlim([0 4])
title('B)')
xlabel('Time (s)')
ylabel('Amplitude (mV)')


%Densidad espectral
periodogram(ECGA,rectwin(length(ECGA)),length(ECGA))
periodogram(ecg,rectwin(length(ecg)),length(ecg))





%%%%%%%%%%%%%FILTRADO FINAL TIPO B%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Fs = 250; %Frecuencia de muestreo
G = 2000; %Numero de muestras



ecg = ecg/G;
ecg = (ecg - mean(ecg))/std(ecg);
t = (1:1:length(ecg))*(1/Fs);

Fnotch = 0.67; % Notch Frequency
BW = 5; % Bandwidth
Apass = 1; % Bandwidth Attenuation
[b, a] = iirnotch (Fnotch/ (Fs/2), BW/(Fs/2), Apass);
Hd1 = dfilt.df2 (b, a);
ecg_limpio1=filter(Hd1,ecg);


Fnotch = 60; % Frecuencia notch
BW = 120; % Ancho banda
Apass = 1; % Atenuacion de ancho banda
[b, a] = iirnotch (Fnotch/ (Fs/2), BW/ (Fs/2), Apass);
Hd2 = dfilt.df2 (b, a);
ecg_limpio2 =filter(Hd2,ecg_limpio1);





windowSize = 5; 
b = (1/windowSize)*ones(1,windowSize);
a=1;
ecg_limpio3= filter(b, a, ecg_limpio2);
ECGB=ecg_limpio3;
subplot(2,1,1)
plot(t,ecg_limpio3)
xlim([0 4])
title('ECG Filtrado Final B')
xlabel('Tiempo (s)')
ylabel('Amplitud (mV)')
subplot(2,1,2)
plot(t,ecg)
xlim([0 4])
title('ECG Original')
xlabel('Tiempo (s)')
ylabel('Amplitud (mV)')









%%%%%%%%%%%%%FILTRADO FINAL TIPO C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Fs = 250; %Frecuencia de muestreo
G = 2000; %Numero de muestras



ecg = ecg/G;
ecg = (ecg - mean(ecg))/std(ecg);
t = (1:1:length(ecg))*(1/Fs);




Fpass = 0.65;        % Passband Frequency
Fstop = 0.667;       % Stopband Frequency
Apass = 0.5;         % Passband Ripple (dB)
Astop = 6;           % Stopband Attenuation (dB)
match = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its CHEBY1 method.
h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, Fs);
Hd3 = design(h, 'cheby1', 'MatchExactly', match);

ecg_limpio1 =filter(Hd3,ecg);



Fpass1 = 59;          % First Passband Frequency
Fstop1 = 59.5;        % First Stopband Frequency
Fstop2 = 60;          % Second Stopband Frequency
Fpass2 = 61;          % Second Passband Frequency
Apass1 = 0.5;         % First Passband Ripple (dB)
Astop  = 3;           % Stopband Attenuation (dB)
Apass2 = 0.5;         % Second Passband Ripple (dB)
match  = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its CHEBY2 method.
h  = fdesign.bandstop(Fpass1, Fstop1, Fstop2, Fpass2, Apass1, Astop, ...
                      Apass2, Fs);
Hd4 = design(h, 'cheby2', 'MatchExactly', match);
ecg_limpio2 =filter(Hd4,ecg_limpio1);





windowSize = 5; 
b = (1/windowSize)*ones(1,windowSize);
a=1;
ecg_limpio3= filter(b, a, ecg_limpio2);
ECGC=ecg_limpio3;
subplot(2,1,1)
plot(t,ecg_limpio3)
xlim([0 4])
title('ECG Filtrado Final C')
xlabel('Tiempo (s)')
ylabel('Amplitud (mV)')
subplot(2,1,2)
plot(t,ecg)
xlim([0 4])
title('ECG Original')
xlabel('Tiempo (s)')
ylabel('Amplitud (mV)')





%%%%%%%%%%%%%FILTRADO FINAL TIPO D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





Fs = 250; %Frecuencia de muestreo
G = 2000; %Numero de muestras



ecg = ecg/G;
ecg = (ecg - mean(ecg))/std(ecg);
t = (1:1:length(ecg))*(1/Fs);



Fpass = 0.65;        % Passband Frequency
Fstop = 0.667;       % Stopband Frequency
Apass = 0.5;         % Passband Ripple (dB)
Astop = 20;           % Stopband Attenuation (dB)
match = 'stopband';  % Band to match exactly


h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, Fs);
Hd5 = design(h, 'ellip', 'MatchExactly', match);
ecg_limpio1 =filter(Hd5,ecg);



Fpass1 = 59;      % First Passband Frequency
Fstop1 = 59.5;    % First Stopband Frequency
Fstop2 = 60;      % Second Stopband Frequency
Fpass2 = 61;      % Second Passband Frequency
Apass1 = 0.5;     % First Passband Ripple (dB)
Astop  = 3;       % Stopband Attenuation (dB)
Apass2 = 0.5;     % Second Passband Ripple (dB)
match  = 'both';  % Band to match exactly

h  = fdesign.bandstop(Fpass1, Fstop1, Fstop2, Fpass2, Apass1, Astop, ...
                      Apass2, Fs);
Hd5 = design(h, 'ellip', 'MatchExactly', match);

ecg_limpio2 =filter(Hd5,ecg_limpio1);




windowSize = 5; 
b = (1/windowSize)*ones(1,windowSize);
a=1;
ecg_limpio3= filter(b, a, ecg_limpio2);
ECGD=ecg_limpio3;
subplot(2,1,1)
plot(t,ecg_limpio3)
xlim([0 4])
title('ECG Filtrado Final D')
xlabel('Tiempo (s)')
ylabel('Amplitud (mV)')
subplot(2,1,2)
plot(t,ecg)
xlim([0 4])
title('ECG Original')
xlabel('Tiempo (s)')
ylabel('Amplitud (mV)')







%%%%%%%%%%%%%%%Espectrogramas de las senales%%%%%%%%%%%%%%%%%%%

Fs=250
Nw=256



 [ywinhat,fw2,t2,P2] = spectrogram(ecg,hann(Nw,'periodic'),Nw/2,Nw,Fs);

  S2 = P2*Fs/Nw;

  SdB2 = 10*log10(S2); % Converir el espectro a escal de decibeles
  subplot(2,1,1)
  image(t2,fw2,SdB2,'CDataMapping','scaled')  

  axis xy
  xlabel('tiempo [s]')
  ylabel('Frecuencia [Hz]')
  title('Espectro ECG sin filtrar')
  colorbar


  [ywinhat,fw2,t2,P2] = spectrogram(ecg_limpio2,hann(Nw,'periodic'),Nw/2,Nw,Fs);

  S2 = P2*Fs/Nw;

  SdB2 = 10*log10(S2); % Converir el espectro a escal de decibeles
  subplot(2,1,2)
  image(t2,fw2,SdB2,'CDataMapping','scaled')  

  axis xy
  xlabel('tiempo [s]')
  ylabel('Frecuencia [Hz]')
  title('Espectro ECG filtrado')
  colorbar

 

figure;
hold on;
plot(t,ECGA)
xlim([0 10])
xlabel('tiempo [s]')
ylabel('Frecuencia [Hz]')
title('ECG filtrada A')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%INTERVALOS RR%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Ubicar los intervalos RR 

% Identificar los RR del filtrado
%save('ecg.mat','ecg','pks','locs');

figure;
hold on;
plot(t,ECGA)
xlim([0 10])
xlabel('tiempo [s]')
ylabel('Frecuencia [Hz]')
title('Identificación de Picos R en señal filtrada A')

umbral_y=6*mean(abs(ECGA));
umbral_x=(60/200)*Fs;


[pksA,locsA]=findpeaks(ECGA,MinPeakHeight=umbral_y,MinPeakDistance=umbral_x);
scatter(t(locsA),pksA)

%Numero de R encontrados con filtado A
length(pksA)
%Calculamos la diferencia entre el numero de R encontrado con Kubios

length(pksA)-length(pks)
%Calculamos el porcentaje de R encontrados en el mismo tiempo


m=abs(locsA-locs);

nnz(~m)/length(pksA)


%Calculamos los RR
RR=diff(pks);
%Obtenemos el valor de SDNN
SDNN=std(RR)*100;




% Identificar los RR del filtrado B%%%%%%%%%%%%%%


figure;
hold on;
plot(t,ECGA)
xlim([0 10])
xlabel('tiempo [s]')
ylabel('Frecuencia [Hz]')
title('Identificación de Picos R en señal filtrada B')

umbral_y=6*mean(abs(ECGA));
umbral_x=(60/200)*Fs;


[pksA,locsA]=findpeaks(ECGA,MinPeakHeight=umbral_y,MinPeakDistance=umbral_x);
scatter(t(locsA),pksA)

%Numero de R encontrados con filtado A
length(pksA)
%Calculamos la diferencia entre el numero de R encontrado con Kubios

length(pksA)-length(pks)
%Calculamos el porcentaje de R encontrados en el mismo tiempo


m=abs(locsA-locs);

nnz(~m)/length(pksA)


%Calculamos los RR
RR=diff(pks);
%Obtenemos el valor de SDNN
SDNN=std(RR)*100;





% Identificar los RR del filtrado C



figure;
hold on;
plot(t,ECGB)
xlim([0 10])
xlabel('tiempo [s]')
ylabel('Frecuencia [Hz]')
title('Identificación de Picos R en señal filtrada C')

umbral_y=6*mean(abs(ECGB));
umbral_x=(60/200)*Fs;


[pksB,locsB]=findpeaks(ECGB,MinPeakHeight=umbral_y,MinPeakDistance=umbral_x);
scatter(t(locsB),pksB)

%Numero de R encontrados con filtado B
length(pksB)
%Calculamos la diferencia entre el numero de R encontrado con Kubios

abs(length(pksB)-length(pks))



%Calculamos los RR
RR=diff(pksB);
%Obtenemos el valor de SDNN
SDNN=std(RR)*100;







% Identificar los RR del filtrado D



figure;
hold on;
plot(t,ECGD)
xlim([0 10])
xlabel('tiempo [s]')
ylabel('Frecuencia [Hz]')
title('Identificación de Picos R en señal filtrada D')

umbral_y=6*mean(abs(ECGD));
umbral_x=(60/200)*Fs;


[pksD,locsD]=findpeaks(ECGD,MinPeakHeight=umbral_y,MinPeakDistance=umbral_x);
scatter(t(locsD),pksD)


%Calculamos la diferencia entre el numero de R encontrado con Kubios

abs(length(pksD)-length(pks))


%Calculamos los RR
RR=diff(pksD);
%Obtenemos el valor de SDNN
SDNN=std(RR)*100;






% Identificar los RR del filtrado D



figure;
hold on;
plot(t,ECGD)
xlim([0 10])
xlabel('tiempo [s]')
ylabel('Frecuencia [Hz]')
title('Identificación de Picos R en señal filtrada')

umbral_y=6*mean(abs(ECGD));
umbral_x=(60/200)*Fs;


[pksD,locsD]=findpeaks(ECGD,MinPeakHeight=umbral_y,MinPeakDistance=umbral_x);
scatter(t(locsD),pksD)

%Numero de R encontrados con filtado B
length(pksD)
%Calculamos la diferencia entre el numero de R encontrado con Kubios

length(pksD)-length(pks)


%Calculamos los RR
RR=diff(pksD);
%Obtenemos el valor de SDNN
SDNN=std(RR)*100;























%%%%%%%%%%%%Carga del banco de datos%%%%%%%%%%%%%%%%%%%






%save ecgdatabank.mat ecg_mv
%load ecgdatabank.mat

load ECGData.mat
x=ECGData.Data(1,:);



%%%%%%%%%%%%%%%%%Cambiamos la frecuencia de muesteo%%%%%%%%%%%

Fs =128 % Sampling Frequency
Fnotch = 0.67; % Notch Frequency
BW = 5; % Bandwidth
Apass = 1; % Bandwidth Attenuation
[b, a] = iirnotch (Fnotch/ (Fs/2), BW/(Fs/2), Apass);
Hd1 = dfilt.df2 (b, a);

Fpass1 = 59;          % First Passband Frequency
Fstop1 = 59.5;        % First Stopband Frequency
Fstop2 = 60;          % Second Stopband Frequency
Fpass2 = 61;          % Second Passband Frequency
Apass1 = 0.5;         % First Passband Ripple (dB)
Astop  = 3;           % Stopband Attenuation (dB)
Apass2 = 0.5;         % Second Passband Ripple (dB)
match  = 'stopband';  % Band to match exactly

h  = fdesign.bandstop(Fpass1, Fstop1, Fstop2, Fpass2, Apass1, Astop, ...
                      Apass2, Fs);
Hd2 = design(h, 'butter', 'MatchExactly', match);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

baja=0;
media=0;
alta=0;

%Creamos un arreglo donde guardar los valores SDNN
SDNN = 1:162;

%Filtramos nuestras ECG (modificando la FS a 500 Hz) del banco de datos con nuestros filtrso

for i=1:162
%ecg_limpio1=filter(Hd1,ECGData.Data(i,:));
ecg_limpio1=filter(Hd1,ECGData.Data(i,:));
ecg_limpio2 =filter(Hd2,ecg_limpio1);
windowSize = 5; 
b = (1/windowSize)*ones(1,windowSize);
a=1;
ecg_mv_filt{i}= filter(b, a, ecg_limpio2);
umbral_y=6*mean(abs(ecg_limpio3));
umbral_x=(60/200)*Fs;
[pks0,locs0]=findpeaks(ecg_mv_filt{i},MinPeakHeight=umbral_y,MinPeakDistance=umbral_x);
%Obtenemos el valor de SDNN
RR=diff(ecg_mv_filt{i})*1000;
SDNN(i)=std(RR);


if SDNN(i) <= 50
     baja=baja+1;
      hrvbaja{baja}=ecg_mv_filt{i};
        elseif SDNN(i) >= 100
             alta=alta+1;
      hrvalta{alta}=ecg_mv_filt{i};
        else
             media=media+1;
      hrvmedia{media}=ecg_mv_filt{i};
end




end



bajaind=find(SDNN<=50);
mediaind=find(SDNN>50 & SDNN<=100);
altaind=find(SDNN>100);



%Vemos cuantas ECG son de HRV baja
baja
%Vemos cuantas ECG son de HRV media
media
%Vemos cuantas ECG son de HRV alta
alta

mediaind
altaind
bajaind


%Generamos los 16 vectores caractersticos para clasificacion de los RR de
%las ECG con HRV baja

for i=1:baja
    n=length(hrvbaja{i})/16;
    m=floor(n);

    temp=hrvbaja{i}(1:j*m);
    RR=diff(temp)*1000;
    
    Vectbaja(i,1)=std(RR);

    for j=2:15
    temp=hrvbaja{i}((j-1):j*m);
    RR=diff(temp)*1000;
    
    Vectbaja(i,j)=std(RR);
        
    end


    temp=hrvbaja{i}(15*m:length(hrvbaja{i}));
    RR=diff(temp)*1000;
    Vectbaja(i,16)=std(RR);
    
 

end












%Generamos los 16 vectores caractersticos para clasificacion de los RR de
%las ECG con HRV media

for i=1:media
    n=length(hrvmedia{i})/16;
    m=floor(n);

    temp=hrvmedia{i}(1:j*m);
    RR=diff(temp)*1000;
    
    Vectmedia(i,1)=std(RR);

    for j=2:15
    temp=hrvmedia{i}((j-1):j*m);
    RR=diff(temp)*1000;
    
    Vectmedia(i,j)=std(RR);
        
    end


    temp=hrvmedia{i}(15*m:length(hrvmedia{i}));
    RR=diff(temp)*1000;
    Vectmedia(i,16)=std(RR);
    
 

end








%Generamos los 16 vectores caractersticos para clasificacion de los RR de
%las ECG con HRV alta

for i=1:alta
    n=length(hrvalta{i})/16;
    m=floor(n);

    temp=hrvalta{i}(1:j*m);
    RR=diff(temp)*1000;
    
    Vectalta(i,1)=std(RR);

    for j=2:15
    temp=hrvalta{i}((j-1):j*m);
    RR=diff(temp)*1000;
    
    Vectalta(i,j)=std(RR);
        
    end


    temp=hrvalta{i}(15*m:length(hrvalta{i}));
    RR=diff(temp)*1000;
    Vectalta(i,16)=std(RR);
    
 

end






for i=1:baja
    n=length(hrvbaja{i})/16;
    m=floor(n);

    temp=hrvbaja{i}(1:j*m);
    RR=diff(temp)*1000;
    
    Vectbaja(i,1)=std(RR);

    for j=2:15
    temp=hrvbaja{i}((j-1):j*m);
    RR=diff(temp)*1000;
    
    Vectbaja(i,j)=std(RR);
        
    end


    temp=hrvbaja{i}(15*m:length(hrvbaja{i}));
    RR=diff(temp)*1000;
    Vectbaja(i,16)=std(RR);
    
 

end




%%%%Coeficientes de predicion lineal%%%%%%%



for i=1:alta

    Vectaltalpc{i}=lpc(hrvalta{i},15);
        
end



for i=1:media
    
    Vectmedialpc{i}=lpc(hrvmedia{i},15);
       
end



for i=1:baja
 
    Vectbajalpc{i}=lpc(hrvbaja{i},15);
        
end



for i=1:60
    
            if i <21
            validar{i}=Vectalta(i,:);
        elseif i>40 
            validar{i}=Vectbaja(i,:);
        else
            validar{i}=Vectmedia(i,:);
            end


end



for i=1:60
    

            Validarmat(i,:)=validar{i};
     
end


%%%%%%%%%LPC envolvente de ECG de muestra%%%%%%%%%%%%



figure;
hold on;
plot(t,ECGA)
xlim([0 10])
xlabel('tiempo [s]')
ylabel('Frecuencia [Hz]')
title('ECG filtrada A')
hold off


x=floor(length(ECGA)/16);
segmentECGA=ECGA(1:x);


coefA=lpc(segmentECGA,16);



hold on
plot(t,ECGA)
xlim([0 10])
[up,lo] = envelope(ECGA);
coefA=up;
plot(t,up,t,lo,'linewidth',1.5)
xlim([0 10])
legend('q','up','lo')
xlabel('tiempo [s]')
ylabel('Frecuencia [Hz]')
title('Envolvente de ECG mediante coefiecientes LPC')
hold off
freqz(ECGA)
hold on
freqz(coefA)
legend('ECGA','Envolvente')
title('ECGA vs su envolvente')
lines = findall(gcf,'type','line');
lines(1).Color = 'red'
lines(2).Color = 'green'
hold off

%%%%%%%%Codebook y prediccion%%%%%%%%%%%%%%%%%%

 [codebook1,pesos1]=vqsplit(transpose(Vectalta),16);

  [codebook2,pesos2]=vqsplit(transpose(Vectmedia),16);


   [codebook3,pesos3]=vqsplit(transpose(Vectbaja),16);




%%%Para los vectores caracteristicos lpc

for i =1:27
    Vectaltalpcm(i,:)=Vectaltalpc{i};
end

[codebook1lpc,pesos1lpc]=vqsplit(transpose(Vectaltalpcm),16);


for i =1:66
    Vectmedialpcm(i,:)=Vectmedialpc{i};
end


[codebook2lpc,pesos2lpc]=vqsplit(transpose(Vectmedialpcm),16);



for i =1:69
    Vectbajalpcm(i,:)=Vectbajalpc{i};
end

[codebook3lpc,pesos3lpc]=vqsplit(transpose(Vectbajalpcm),16);


%Validamos nuestros codebooks (clusters) en la  clasificacion de las
%muestras de entrenamiento.


matconf=zeros(3);
for i=1:60

    distances1 = pdist2(validar{1}, codebook1);
    distances2 = pdist2(validar{i}, codebook2);
    distances3 = pdist2(validar{i}, codebook3);
    minDistance = [min(distances1),min(distances2),min(distances3)];
    [minim,ind]=min(minDistance);



                if i <21
                    clasificar{i}=[ind,1];
                    if ind==1
                        matconf(1,1)=matconf(1,1)+1;
                    elseif ind==2
                        matconf(1,2)=matconf(1,2)+1;
                    else
                        matconf(1,3)=matconf(1,3)+1;
                    end


        elseif i>40 
           clasificar{i}=[ind,3];
                               if ind==1
                        matconf(3,1)=matconf(3,1)+1;
                    elseif ind==2
                        matconf(3,2)=matconf(3,2)+1;
                    else
                        matconf(3,3)=matconf(3,3)+1;
                    end
        else
            clasificar{i}=[ind,2];
                                if ind==1
                        matconf(2,1)=matconf(2,1)+1;
                    elseif ind==2
                        matconf(2,2)=matconf(2,2)+1;
                    else
                        matconf(2,3)=matconf(2,3)+1;
                    end
            end

end




Y = {'Alta'; 'Media'; 'Baja'};


cm = confusionchart(matconf,Y,'Title','Clasificacion ECG a partir de HRV','RowSummary','row-normalized','ColumnSummary','column-normalized')


%%%Para los vectores caracteristicos de lpc%%%%%%%%%%%%%%%

matconf=zeros(3);
for i=1:60

    distances1 = pdist2(validar{1}, codebook1lpc);
    distances2 = pdist2(validar{i}, codebook2lpc);
    distances3 = pdist2(validar{i}, codebook3lpc);
    minDistance = [min(distances1),min(distances2),min(distances3)];
    [minim,ind]=min(minDistance);



                if i <21
                    clasificar{i}=[ind,1];
                    if ind==1
                        matconf(1,1)=matconf(1,1)+1;
                    elseif ind==2
                        matconf(1,2)=matconf(1,2)+1;
                    else
                        matconf(1,3)=matconf(1,3)+1;
                    end


        elseif i>40 
           clasificar{i}=[ind,3];
                               if ind==1
                        matconf(3,1)=matconf(3,1)+1;
                    elseif ind==2
                        matconf(3,2)=matconf(3,2)+1;
                    else
                        matconf(3,3)=matconf(3,3)+1;
                    end
        else
            clasificar{i}=[ind,2];
                                if ind==1
                        matconf(2,1)=matconf(2,1)+1;
                    elseif ind==2
                        matconf(2,2)=matconf(2,2)+1;
                    else
                        matconf(2,3)=matconf(2,3)+1;
                    end
            end

end




Y = {'Alta'; 'Media'; 'Baja'};


cm = confusionchart(matconf,Y,'Title','Clasificacion ECG a partir de LPC','RowSummary','row-normalized','ColumnSummary','column-normalized')





%%%%%%%%%%%%Random Forest%%%%%%%%%%%%%%%%5

%%%Preparamos el conjunto de datos y vector de respuesta

Randdatf=zeros(16);
Yres=(1:169);


for i =1:27
    Randdatf(i,:)=Vectalta(i,:);
    Yres(i)=1;
end


for i =1:66
    Randdatf(i+27,:)=Vectmedia(i,:);
    Yres(i+27)=2;
end


for i =1:69
    Randdatf(i+93,:)=Vectbaja(i,:);
    Yres(i+93)=3;
end



%load('modelo.mat')
 polyre = reshape(ecg_mv_filt,[162 1]);
ECG=cell2mat(polyre);

M=[]
