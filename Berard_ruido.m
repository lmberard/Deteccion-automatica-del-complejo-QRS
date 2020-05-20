%TP SEÑALES Y SISTEMAS. 
%LUCIA BERARD 101213

%%
%ANALISIS DE LA SEÑAL

%Importo los datos de la señal original
original = importdata('103m.mat');
marcas = importdata('marcas.txt');

original = original-mean(original); %offset
original = original ./200; %ganancia

SNR = 10;
datos = awgn(original, SNR);

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

%5 ciclos aprox de QRS
xmin = 649;
xmax = 653; 

figure 
    hold on
    plot(time,original); 
    xlim([xmin xmax]);
    plot(time,datos); 
    grid on;
    legend('Señal original', 'Señal con ruido');
    ylabel("Amplitud [mV]"); xlabel("Tiempo [s]");
    hold off
    
marcas = marcas./100;

%%
%PREPROCESAMIENTO DE LA SEÑAL

%FILTRO PASABAJOS-----------------------------------------------
num_low = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
den_low = [1 -2 1];
N = 2048;

y_low = filter(num_low, den_low, datos);
pendiente_low = 5;   
    
%FILTRO PASAALTOS-----------------------------------------------
num_high = [(-1/32) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 (1/32)];
den_high = [1 -1];

y_high = filter(num_high, den_high, y_low); %Filtrado
pendiente_high = 16;
  
%DERIVADOR---------------------------
num_deriv = [1 2 0 -2 -1];
den_deriv = [8];

y_derivador = filter(num_deriv, den_deriv, y_high);
pendiente_derivador = 2;   
   
%CUADRADO----------------------------------
y_cuad_t= y_derivador.*y_derivador;

%INTEGRADOR---------------------------------
N_integrador = 30;
num_integrador = ones(1,N_integrador)/N_integrador;
den_integrador = [1];

pendiente_int = 17;   
y_integrador = filter(num_integrador, den_integrador, y_cuad_t);

pendiente_total = pendiente_int + pendiente_derivador + pendiente_low + pendiente_high;
y_integrador = delayseq(y_integrador',-pendiente_total);

figure
    hold on
    plot(time,y_integrador);
    plot(time,datos);
    xlim([xmin xmax]);
    grid on; 
    ylabel('Amplitud[mV]');xlabel('Tiempo[s]');
    legend('Señal preprocesada','Señal original con ruido');
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

y_final = zeros(1, length(y_integrador));

for i = 1:length(LOCS)
    posicion_vector = LOCS(i);
    y_final(posicion_vector) = PKS(i);
end

%GOLD STANDARD------------------------------
figure 
    hold on 
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
