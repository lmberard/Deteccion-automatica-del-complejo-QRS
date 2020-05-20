%TP SEÑALES Y SISTEMAS. 
%LUCIA BERARD 101213

%%
%ANALISIS DE LA SEÑAL

%Importo los datos de la señal original
original = importdata('103m.mat');
marcas = importdata('marcas.txt');

fs = 200;
fs_nueva= 360; 

%Necesito la relacion para calcular la expansion y decimador asique la
%calculo asi que me separa numerador y denominador
[expansor,decimador] = rat(fs_nueva/fs); 
n = (-8000:8000);
h_interpolador = sinc(n/expansor);

%aumento la cantidad de muestras de la señal original
original = resample(original, expansor, decimador , h_interpolador); 
original = original-mean(original); %offset
original = original ./200; %ganancia

SNR = 30;
datos = awgn(original, SNR);

%porque la señal tiene 361112 muestras y datos tiene 650002 (1.8 veces mas)
xmin = 450 * 1.8;
xmax = 453 * 1.8;

datos_min = x_min*200;
datos_max = x_max*200;

T = length(datos);
time = (0:T-1)/new_freq;
T_marcas = T*100;
time_marcas =(0:T_marcas-1)/new_freq;
vector_marcas = zeros (1,T_marcas);
time_marcas = time_marcas ./100;
marcas = marcas.*100;
marcas = round(marcas);


for i = 1:length(marcas)
    indice = marcas(i);
    vector_marcas (indice) = 0.7;
end

figure
    hold on
    plot(time,datos);
    plot(time,original);
    %xlim([0 1800]);
    xlim([x_min x_max]);
    legend('Señal con ruido', 'Señal original');
    ylabel('ECG[mV]');xlabel('Tiempo[s]');
    grid on
    hold off
    
marcas = marcas./100;
marcas = marcas.*1.8;

%%
%PREPROCESAMIENTO DE LA SEÑAL

%FILTRO PASABAJOS-----------------------------------------------
num_low = [1 2 3 4 5 6 5 4 3 2 1]; 
den_low= [1]; 
N = 2048;
imp_low_viejo = impz(num_low,den_low);
imp_low_nuevo= resample(imp_low_viejo, expansor, decimador, h_interpolador);

y_low = filter(imp_low_nuevo,1,datos);

%FILTRO PASAALTOS-----------------------------------------------
num_high = [(-1/32) (-1/32) (-1/32) (-1/32) (-1/32) (-1/32) (-1/32) (-1/32) (-1/32) (-1/32) (-1/32) (-1/32) (-1/32) (-1/32) (-1/32) (-1/32) (31/32) (-1/32) (-1/32) (-1/32) (-1/32) (-1/32) (-1/32) (-1/32) (-1/32) (-1/32) (-1/32) (-1/32) (-1/32) (-1/32) (-1/32) (-1/32)]; 
den_high = [1]; 
imp_high_viejo = impz(num_high,den_high);
imp_high_nuevo = resample(imp_high_viejo, expansor, decimador, h_interpolador);


y_high = filter(num_high,den_high,y_low); 

%DERIVADOR---------------------------
num_deriv = [1/8 2/8 0 -2/8 -1/8];
den_deriv = [1];
y_derivador = filter(num_deriv, den_deriv, y_high);

%CUADRADO----------------------------------
y_cuad_t= y_derivador.*y_derivador;

%INTEGRADOR---------------------------------
N_integrador = 54;
num_integrador = ones(1,N_integrador);
num_integrador = num_integrador./N_integrador;
den_integrador = [1];
y_integrador = filter(num_integrador, den_integrador, y_cuad_t);
 


%%
%DETECCION DE LOS QRS:
pendiente_total = 60;
y_integrador = delayseq(y_integrador',-pendiente_total);

%FINDPEAKS-----------------------------------
[PKS,LOCS]= findpeaks(y_integrador);

y_final = zeros(1, T);

%Construyo un vector de ceros que tenga solo los picos
for i = 1:length(LOCS)
    posicion_vector = LOCS(i);
    y_final(posicion_vector) = PKS(i);
end

%GOLD STANDARD------------------------------
posiciones_qrs = qrs_detection(y_final);

figure 
    hold on
    plot(time,datos);
    plot(time,y_integrador);
    plot(time,posiciones_qrs);
    xlim([xmin xmax]);
    ylabel('Amplitud[mV]');xlabel('Tiempo[s]');
    legend('Señal con ruido','Señal preprocesada','Qrs Detection');
    grid on
    hold off

%VERIFICACION DETECTADOS--------------------
%cant_detectada = length(find(posiciones_qrs))
[valor,posicion] = findpeaks(posiciones_qrs);
cant_detectadas = length(posicion)
cant_real = length(marcas)

detectados = 0;
tolerancia = 18*360/200;

 for i=1:1:cant_real
     for j=1:cant_detectadas
        if (abs((marcas(i)-posicion(j)))) <= tolerancia 
               detectados = detectados + 1;
        end
     end
 end
 
detectados
falsos_negativos = cant_real - detectados
falsos_positivos = cant_detectadas - detectados 

%EFECTO DE TODOS LOS FILTRADOS-------------------------------------------
% figure
% subplot(3,2,1); plot(time,datos); xlim([xmin xmax]);title ('Señal original'); grid on; ylabel('Amplitud[mV]');xlabel('Tiempo[s]');
% subplot(3,2,2); plot(time,y_low); xlim([xmin xmax]);title ('Pasabajos'); grid on; ylabel('Amplitud[mV]');xlabel('Tiempo[s]');
% subplot(3,2,3); plot(time,y_high); xlim([xmin xmax]);title ('Pasaaltos'); grid on; ylabel('Amplitud[mV]');xlabel('Tiempo[s]');
% subplot(3,2,4); plot(time,y_derivador); xlim([xmin xmax]);title ('Derivador');grid on; ylabel('Amplitud[mV]');xlabel('Tiempo[s]');
% subplot(3,2,5); plot(time,y_cuad_t); xlim([xmin xmax]);title ('Cuadrado');grid on; ylabel('Amplitud[mV]');xlabel('Tiempo[s]');
% subplot(3,2,6); plot(time,y_integrador); xlim([xmin xmax]);title ('Integrador');grid on; ylabel('Amplitud[mV]');xlabel('Tiempo[s]');


