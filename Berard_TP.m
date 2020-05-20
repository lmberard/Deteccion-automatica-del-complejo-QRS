%TP SEÑALES Y SISTEMAS. 
%LUCIA BERARD 101213

%%
%ANALISIS DE LA SEÑAL

%Importo los datos de la señal original
datos = importdata('103m.mat');
marcas = importdata('marcas.txt');

datos = datos - mean(datos); %offset
datos = datos ./200; %ganancia

%Para graficarlo en funcion del tiempo
fs = 200;
cantidad = length(datos);
T_marcas = cantidad*100;
time_marcas = (0:T_marcas - 1)/fs;
vector_marcas = zeros (1,T_marcas);
time_marcas = time_marcas ./100;
marcas = marcas.*100;
marcas = round(marcas);

for i = 1:length(marcas)
    indice = marcas(i);
    vector_marcas(indice) = 0.1;
end

time = (0:cantidad-1)/fs;
figure %Grafico la señal de entrada entera en tiempo y en marcas
subplot(2,1,1); 
    plot (time,datos); 
    grid on;
    xlim ([0 1800]);
    ylabel("Amplitud [mV]");xlabel ("Tiempo [s]");%1800 segundos son 30 minutos
subplot(2,1,2); 
    plot(datos); 
    grid on;
    ylabel("Amplitud [mV]");xlabel ("Marcas [n]");

marcas = marcas./100;

%%Agarro solo una QRS
qrt = datos(48656:64656); %16.000 muestras 
%tengo que tomar mucho mas ciclos para poder apreciar mas adelante los armonicos cuando hago la transformada
%Cada ciclo tiene aprox 160 muestras con un ancho del qrs de aprox 20
%muestras => 16000 muestras son aprox 100 ciclos. Esto me asegura que al
%aplicarle la transformada de fourier se puedan observar mejor los
%armonicos y despues en el espectrograma voy a poder observar los ruidos
%estacionarios

%%TRANSFORMADA RAPIDA DE FOURIER
figure %2
    y = fft(qrt);
    %crear eje de tiempo
    freq = 200*linspace(-0.5, 0.5, length(y));
    %acomoda el eje y
    response = abs(fftshift(y./max(abs(y))));
    plot(freq,response);
    xlim([0, 100]) %para mostrar solo la mitad
    xlabel('Muestras[n]'); ylabel('ECG[mV]');
    grid on

%ESPECTROGRAMA------------------------------------
%frecuencia de muestreo (dato)
%hay que ver las lineas horizontales en 60hz y la de los armonicos!!! que
%tiene que estar relacionado con el fft que grafique antes

figure %3
    window = hann(1000);
    spectrogram(datos, window,1,[], fs, 'yaxis');
    %[s,f,t] = spectrogram(x,window,noverlap,f,fs)
    %title('Espectrograma');

%%el ruido que aparece en los 60hz es el de la electricidad = ruido
%%estacionario

%%el otro ruido que se da con los picos se debe a un error humano o algun
%%problema con el instrumento (no es periodico, es un ruido no
%%estacionario)

%%
%PREPROCESAMIENTO DE LA SEÑAL

%5 ciclos aprox de QRS
xmin = 649;
xmax = 653; 


%FILTRO PASABAJOS-----------------------------------------------
num_low = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
den_low = [1 -2 1];
N = 2048;

%Rta en frecuencia
figure %4
    freqz(num_low, den_low, N);
    
%Diagrama polos y ceros 
figure %5
    zplane(num_low, den_low);
    
%Rta al impulso
figure %6
    impz(num_low, den_low);
    grid on;
     
y_low = filter(num_low, den_low, datos);
pendiente_low = 5;   


figure %3
    window = hann(1000);
    spectrogram(y_low, window,1,[], fs, 'yaxis');

    
%FILTRO PASAALTOS-----------------------------------------------
num_high = [(-1/32) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 (1/32)];
den_high = [1 -1];

%Rta en frecuencia
figure %8
    freqz(num_high, den_high,N);

%Diagrama polos y ceros
figure %9
    zplane(num_high, den_high);

%Rta al impulso
figure %10
    impz(num_high, den_high);
    grid on;
    
%Efecto del filtro pasaaltos
y_high = filter(num_high, den_high, y_low); %Filtrado
pendiente_high = 16; %pendiente de la fase del diagrama de bode

figure %3
    window = hann(1000);
    spectrogram(y_low, window,1,[], fs, 'yaxis');

%DERIVADOR---------------------------
num_deriv = [1 2 0 -2 -1];
den_deriv = [8];

%Rta en frecuencia
figure %12
    freqz(num_deriv, den_deriv, N);

%Efecto del filtro derivador
y_derivador = filter(num_deriv, den_deriv, y_high);
pendiente_derivador = 2;   
%zero_deriv= zeros(1, pendiente_deriv);
%y_derivador = cat(2, y_derivador, zero_deriv);

    
%CUADRADO----------------------------------
y_cuad_t= y_derivador.*y_derivador;

figure 
    plot(time,y_cuad_t);
    xlim([xmin xmax]);
    grid on; 
    ylabel('ECG[mV]');xlabel('Time[s]');
    
figure %Grafico la respuesta en frecuencia de la señal al cuadrado 
    window = hann(1000); % ventana de Hanning idéntica a la anterior
    spectrogram(y_cuad_t, window, 1, [], fs, 'yaxis');


%INTEGRADOR---------------------------------
N_integrador = 30;
num_integrador = ones(1,N_integrador)/N_integrador;
den_integrador = [1];

%Rta en frecuencia
figure %15
    freqz(num_integrador, den_integrador, N);

pendiente_int = 17;   
y_integrador = filter(num_integrador, den_integrador, y_cuad_t);

pendiente_total = pendiente_int + pendiente_derivador + pendiente_low + pendiente_high;
%zeros_int= zeros(1, pendiente_int);
%y_integrador = cat(2, y_integrador, zeros_int); 

y_integrador = delayseq(y_integrador',-pendiente_total);

figure %3
    window = hann(1000);
    spectrogram(y_low, window,1,[], fs, 'yaxis');


figure
    hold on
    plot(time,y_integrador);
    plot(time,datos);
    xlim([xmin xmax]);
    grid on; 
    ylabel('Amplitud[mV]');xlabel('Tiempo[s]');
    legend('Señal preprocesada','Señal original');
    hold off

%EFECTO DE TODOS LOS FILTRADOS-------------------------------------------
figure %17
subplot(3,2,1); plot(time,datos); xlim([xmin xmax]);title ('Señal original'); grid on; ylabel('Amplitud[mV]');xlabel('Tiempo[s]');
subplot(3,2,2); plot(time,y_low); xlim([xmin xmax]);title ('Pasabajos'); grid on; ylabel('Amplitud[mV]');xlabel('Tiempo[s]');
subplot(3,2,3); plot(time,y_high); xlim([xmin xmax]);title ('Pasaaltos'); grid on; ylabel('Amplitud[mV]');xlabel('Tiempo[s]');
subplot(3,2,4); plot(time,y_derivador); xlim([xmin xmax]);title ('Derivador');grid on; ylabel('Amplitud[mV]');xlabel('Tiempo[s]');
subplot(3,2,5); plot(time,y_cuad_t); xlim([xmin xmax]);title ('Cuadrado');grid on; ylabel('Amplitud[mV]');xlabel('Tiempo[s]');
subplot(3,2,6); plot(time,y_integrador); xlim([xmin xmax]);title ('Integrador');grid on; ylabel('Amplitud[mV]');xlabel('Tiempo[s]');


%%
%DETECCION DE LOS QRS:

%FINDPEAKS-----------------------------------
[PKS,LOCS]= findpeaks(y_integrador);
%findpeaks me devuelve dos vectores, PKS donde estan los valores de los
%maximos, y LOCS la posicion en el tiempo de las posiciones

%creo un vector de ceros y en cada posicion de los maximos, agrego el valor
%correspondiente (recorro con un for)
y_final = zeros(1, length(y_integrador));

for i = 1:length(LOCS)
    posicion_vector = LOCS(i);
    y_final(posicion_vector) = PKS(i);
end


%GOLD STANDARD------------------------------
figure %17
    hold on
    %corrimiento_total = 42;   
    %zero_total= zeros(1, corrimiento_total);
    %datos = cat(2, datos, zero_total); 
    posiciones_qrs = qrs_detection(y_final);
    plot(time,datos);
    plot(time,y_integrador);
    plot(time,y_final);
    plot(time,posiciones_qrs);
    xlim([xmin xmax]);
    ylabel('Amplitud[mV]');xlabel('Tiempo[s]');
    legend('Señal original', 'Integrador','Findpeaks','Qrs Detection');
    grid on
    hold off

%VERIFICACION DETECTADOS--------------------
cant_marcas_reales = length(marcas)

[valor_detectados,ubic_detectados]= findpeaks(posiciones_qrs);
cant_marcas_detector = length(valor_detectados);

detectados = 0;
tolerancia = 18;

 for i=1:1:cant_marcas_reales
     for j=1:cant_marcas_detector
        if (abs(round(marcas(i)-ubic_detectados(j)))) <= tolerancia 
               detectados = detectados + 1;
        end
     end
 end
 
detectados
falsos_negativos = cant_marcas_reales - detectados
falsos_positivos = cant_marcas_detector - detectados 

