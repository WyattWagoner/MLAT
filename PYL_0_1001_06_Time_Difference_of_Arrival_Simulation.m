 %% Time Difference of Arrival Simulation Code
    % Creator name: Griffin Curtis
    % Date Created: Feb 2016
    % Last Date Modified: 2016-11-07
    %
    % Description:
    % Calculates the Time Difference of Arrival (TDOA) of a transmission
    % being recieved at 4 different receivers. The receivers consist of 3
    % ground receivers and 1 LEO satellite. The transmitter is Echostar 14.
    % TDOAs are then sent to a multilateration code to calculate the
    % position of the target satellite based off of time delays of between
    % receivers and receiver coordinates.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Close figures
    close all;  
    clear all;
%% Define global variables
    global est_tdoa_with_noise1
    global est_tdoa_with_noise2
    global est_tdoa_with_noise3
    global est_tdoa_with_noise4
    global est_tdoa_with_noise5
    global est_tdoa_with_noise6
    global R_a_lla
    global R_b_lla
    global R_c_lla
    global R_d_lla
    global R_a 
    global R_b 
    global R_c
    global R_d
    global GEO_coordinates
    global GEO_lla
    
%% Receiver and transmitter positions
%    R_a_lla = [58.37, -134.2, 450000];             %Over Alaska!
%    R_a_lla = [45.996, -116.917, 450000];          %South East WA
%    R_a_lla = [45.996, -116.917, 265];             %South East WA
%    R_a_lla = [46.284, -124.062, 450000];          %South West WA
%    R_a_lla = [48.9798, -123.224, 450000];         %North West WA
%    R_a_lla = [48.975, -117.0404, 450000];         %North East WA
%    R_a_lla = [47.62, -122.4, 450000];             %Over Seattle
%    R_a_lla = [38.612435, -98.575397, 450000];     %Over Kansas!
    R_a_lla = [21.13062, -156.85729, 450000];      %Hawaii Space
%    R_a_lla = [35.08533, -106.60555, 450000];      %Albuquerque Space
    R_b_lla = [25.769322, -80.213916, 1.8];         %Miami
%    R_c_lla = [32.720849, -117.162867, 283];       %San Diego
%    R_c_lla = [47.62, -122.4, 158];                %Seattle
    R_c_lla = [35.08533, -106.60555,1619];          %Albuquerque
    R_d_lla = [47.119482, -88.563732, 196];         %Houghton
 
    GEO_lla = [0, -119, 35786000];                                  %GEO position
    R_a = lla2ecef(R_a_lla, 'WGS84');                               %[1x3] 3D position of receiver "a" Auris
    R_b = lla2ecef(R_b_lla, 'WGS84');                               %[1x3] 3D position of receiver "b" Miami
    R_c = lla2ecef(R_c_lla, 'WGS84');                               %[1x3] 3D position of receiver "c" San Diego
    R_d = lla2ecef(R_d_lla, 'WGS84');                               %[1x3] 3D position of receiver "d" Houghton
    GEO_coordinates = lla2ecef(GEO_lla, 'WGS84');                   %[1x3] 3D position of Transmitter "GEO" 
  
%% LEO to GEO distances and times
    light= 299792458;                               %speed of light (m/s)
    D1=norm(GEO_coordinates - R_a);                 %m
    D2=norm(GEO_coordinates - R_b);                 %m
    D3=norm(GEO_coordinates - R_c);                 %m
    D4=norm(GEO_coordinates - R_d);                 %m

%% Data Rate Variables
    Mbps=50*3/2;
    bps=Mbps*1000000;
    Time_per_bit= 1/bps;

%% Calculates times of arrival, calculates tdoa, places into a time shift vector
    T1=D1/(light);                                                  %s
    T2=D2/(light);                                                  %s
    T3=D3/(light);                                                  %s
    T4=D4/(light);                                                  %s
    time_shift= [(T2-T1),(T3-T1),(T4-T1),(T3-T2),(T4-T2),(T4-T3)];  %s
    frac_err =.01;                                                  % fraction of bits to be randomly flipped in both signals    
    process_duration = .020;                                        % time, in seconds, to grab data
    timer_accuracy = 20e-9;                                         % time accuracy of system
%% Create Data vectors received at satellites
    sampleLength=ceil(process_duration/Time_per_bit+2*max(time_shift)/Time_per_bit);
    cen = ceil(sampleLength/2);
    Data = false(4,sampleLength);                          %4 x nbits matrix

    %set each element of row 1 to 0 or 1 each row is data from each reciever
    Data(1,:) = rand(sampleLength,1) >= 0.5; 
 
%% Shift Data 
    %Shift data according to time_shift vector. Round is used 
    %to create an integer number for data, and reflects our time accuracy
    %FIX THIS
    %timer_error = timer_accuracy * randn(1,6);
    iteration = 1;
    for i = 1:length(time_shift)
        shiftIndex =round((time_shift(i))/Time_per_bit);   
        if time_shift(i) >= 0
            Data(i+1,(shiftIndex + 1):sampleLength) = Data(1,1:(sampleLength-shiftIndex));
            Data(i+1,1:shiftIndex) = double(rand(shiftIndex,1) >= 0.5);
        else
            shiftIndex = abs(shiftIndex);
            Data(i+1,1:(sampleLength-shiftIndex)) = Data(1,(shiftIndex + 1):sampleLength);
            Data(i+1,(sampleLength-shiftIndex+1):sampleLength) = double(rand(shiftIndex,1) >= 0.5);
        end
    end

%% Add Noise, Take Sample, Correlate, Filter, MLAT
    %Noise is added by creating new vectors called NsyDat1->NsyDat4. This saves
    %the original data so that it can be used to produce different levels of 
    %Bit Error Rate. Vectors are used instead of arrays to allow the 
    %SampleSizeIndex loop to be put into a parfor loop for parrallel processing 
    NsyDat = Data;
    fliplist = rand(size(NsyDat, 1), size(NsyDat, 2)) <= frac_err;
    NsyDat1=double(xor(NsyDat(1,:), fliplist(1,:))); 
    NsyDat2=double(xor(NsyDat(2,:), fliplist(2,:))); 
    NsyDat3=double(xor(NsyDat(3,:), fliplist(3,:))); 
    NsyDat4=double(xor(NsyDat(4,:), fliplist(4,:))); 
    
    csvwrite('Receiver1.csv',NsyDat1);
    csvwrite('Receiver2.csv',NsyDat2);
    csvwrite('Receiver3.csv',NsyDat3);

   %Cross-correllation was replaced with finddelay. It does the cross correlation for us. Very useful. 
    est_tdoa_with_noise1 =finddelay(NsyDat1,NsyDat2)/bps;
    fprintf('input tdoa = %11.10f\n', time_shift(1))
    fprintf('estimate of tdoa with noise = %11.10f\n', est_tdoa_with_noise1)

    est_tdoa_with_noise2 = finddelay(NsyDat1,NsyDat3)/bps;
    fprintf('input tdoa = %11.10f\n', time_shift(2))
    fprintf('estimate of tdoa with noise = %11.10f\n', est_tdoa_with_noise2)

    est_tdoa_with_noise3 = finddelay(NsyDat1,NsyDat4)/bps;
    fprintf('input tdoa = %11.10f\n', time_shift(3))
    fprintf('estimate of tdoa with noise = %11.10f\n', est_tdoa_with_noise3)


    est_tdoa_with_noise4 = finddelay(NsyDat2,NsyDat3)/bps;
    fprintf('input tdoa = %11.10f\n', time_shift(4))
    fprintf('estimate of tdoa with noise = %11.10f\n', est_tdoa_with_noise4)

    est_tdoa_with_noise5 = finddelay(NsyDat2,NsyDat4)/bps;
    fprintf('input tdoa = %11.10f\n', time_shift(5))
    fprintf('estimate of tdoa with noise = %11.10f\n', est_tdoa_with_noise5)

    est_tdoa_with_noise6 = finddelay(NsyDat3,NsyDat4)/bps;
    fprintf('input tdoa = %11.10f\n', time_shift(6))
    fprintf('estimate of tdoa with noise = %11.10f\n', est_tdoa_with_noise6)

        
%% Run multilateration
    %All variables needing to be passed are globals right now. MLAT
    %uses nested functions, making it tricky to pass the proper
    %variables. 
    Transmitter=PYLMultilateration00();

