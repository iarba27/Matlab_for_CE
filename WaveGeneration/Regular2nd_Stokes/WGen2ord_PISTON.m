%This script compares the SPH measured velocities with theoretical ones for
%Regulas waves and Calculate the time series of a PISTON-type wavemaker for
%1st and 2nd order wave generation theory (Madsen approximation, 1971)

%%
%set the input parameters for SPH and CIEM

clc
clear all
close all

f_out = 40;        %sample SPH frequency (number of fps)
nADV = 1;                  %number of ADV in SPH
tmax = 60;                %duration SPH run
name='reg01';     
namefolder=strcat('CaseCIEM_Globos_',num2str(name));
% timef  to apply the ramp-up
ramp=2;
rf=round(ramp*f_out); %seconds of ramp-up
quadd=linspace(0,1,rf).^2; %quadratic ramp

% for elevate piston type, h0= lenght of the fixed part.
h0=0.8673; %CIEM pos enganche 2
% h0=0;

% H=0.02; %wave height
% T=1.5; %wave period
% d=1.51; %water depth

H=0.5; %wave height
T=4.041; %wave period
d=2.47; %water depth

L=L_lin(T,d);
kL=2*pi/L; %wave number
f=1/T; %frequency
x=0; %x coordinates ADV
z=0; %z coordinate ADV
g=9.81; %gravity acceleration
ph0=0; %initial phase

%time min and max for plots
tplotmin=0;
tplotmax=4*T;

%% for PISTON wave generation 

% m1=4*(sinh(kL*d)^2)/(sinh(2*kL*d)+2*kL*d); %1st order Biesel function S/H
m1=4*((sinh(kL*d)^2)-sinh(kL*h0)*sinh(kL*d))/(sinh(2*kL*d)+2*kL*d); %1st order Biesel function S/H
s0=H/m1; %piston max Stroke 
o2=(H^2)/(32*d)*(3*cosh(kL*d)/(sinh(kL*d)^3)-2/m1); %2nd order component

% check for Madsen (1971) theory Madcrit<Madlim
Madlim=8*pi*pi/3
Madcrit=H*L*L/d^3

%% set timeline for SPH
row=f_out*tmax;

for i=1:row;
    time(i)=i/f_out;
end

%% theoretical velcoities

%1st and 2nd order wave celerity
C=g*T/(2*pi)*tanh(kL*d); 

%3rd order wave celerity
% C=C*(1+(pi*H/L)^2*((5+2*cosh(2*kL*d)+2*cosh(2*kL*d)^2)/(8*sinh(kL*d)^4)));

for i=1:row;
    
        %% piston displacement 
        e1(i,1)=0.5*s0*sin(2*pi*f*time(i)+ph0); %1st term
        e2(i,1)=o2*sin(4*pi*f*time(i)+2*ph0); %2nd order term (Madsen, 1971)
        e(i,1)=e1(i,1)+e2(i,1); %total signal 
        
        %% old stroke and piston movement(linear theory)
%         s0_old=H/((2*sinh(kL*d)^2)/(sinh(kL*d)*cosh(kL*d)+kL*d));
%         e1_old(i,1)=0.5*s0_old*sin(2*pi*f*time(i)+ph0);
    
        %% incident
        % 2nd order Stokes incident
        eta_inc(i,1)=H/2*cos(2*pi*f*time(i)-kL*x+ph0)+(kL*H^2)/16*((3*((coth(kL*d))^3))-coth(kL*d))*cos(4*pi*f*time(i)-2*kL*x+2*ph0);   
        u_inc(i,1)=H/2*2*pi*f*cosh(kL*(d+z))*(cos(2*pi*f*time(i)-kL*x+ph0))/sinh(kL*d)+3/4*(pi*H/L)^2*C*cosh(2*kL*(d+z))*(cos(4*pi*f*time(i)-2*kL*x+2*ph0))/sinh(kL*d)^4;
        v_inc(i,1)=-H/2*2*pi*f*sinh(kL*(d+z))*(sin(2*pi*f*time(i)-kL*x+ph0))/sinh(kL*d)-3/4*(pi*H/L)^2*C*sinh(2*kL*(d+z))*(sin(4*pi*f*time(i)-2*kL*x+2*ph0))/sinh(kL*d)^4;
        
        % only second component Stokes incident for eta
        eta_inc2(i,1)=+(kL*H^2)/16*((3*((coth(kL*d))^3))-coth(kL*d))*cos(4*pi*f*time(i)-2*kL*x+2*ph0); 
        
        
        %2nd order Stokes total
        eta_2(i,1)=eta_inc(i,1);
        u_2(i,1)=u_inc(i,1);
        v_2(i,1)=v_inc(i,1);
        
        % 1st order (Airy)
        eta(i,1)=H/2*cos(2*pi*f*time(i)-kL*x+ph0);  
        u(i,1)=H/2*2*pi*f*cosh(kL*(d+z))*(cos(2*pi*f*time(i)-kL*x+ph0))/sinh(kL*d);
        v(i,1)=-H/2*2*pi*f*sinh(kL*(d+z))*(sin(2*pi*f*time(i)-kL*x+ph0))/sinh(kL*d); 
end


%% wave board force

A=2*2*pi*f*s0*sinh(kL*d)/(kL*(sinh(2*kL*d)+2*kL*d));

Fr_max=(1000*2*pi*f*A*sinh(kL*d)/kL);
Ft_max=2*Fr_max; %if water on both piston sides and neglecting evanescent modes


%% Figure
figure(1)
figsize1 = [100 100 1300 900];
figure(1); clf(figure(1)); set(gcf, 'color', 'white','Position', figsize1);


    subplot(2,2,1)
    plot(time,u(:,1),'b'); hold on 
    plot(time,u_2(:,1),'r')    
    xlabel('t [s]')
    ylabel('u [m/s]')
    xlim([tplotmin tplotmax])
    legend('Theoretical (Airy)','Theoretical (2nd order Stokes)' )
%     ylim([-0.5 0.5])
    grid on
    
    subplot(2,2,2)
    plot(time,v(:,1),'b');hold on 
    plot(time,v_2(:,1),'r')  
    xlabel('t [s]')
    ylabel('v [m/s]')
    xlim([tplotmin tplotmax])
%     ylim([-0.4 0.4])
    legend('Theoretical (Airy)','Theoretical (2nd order Stokes)' )
    grid on
    
    subplot(2,2,3)
    plot(time,eta(:,1),'b');hold on 
    plot(time,eta_2(:,1),'r')
%     plot(time,eta_incSPM(:,1),'--k')
%     plot(time,eta_inc2(:,1),'--m') 
    xlabel('t [s]')
    ylabel('\eta [m]')
    xlim([tplotmin tplotmax])
    ylim([-1*H 1*H])
    grid on
    legend('Theoretical (Airy)','Theoretical (2nd order Stokes)' )
    
    subplot(2,2,4)
    plot(time,e(:,1),'--k'); hold on
    plot(time,e1(:,1),'r'); 
    plot(time,e2(:,1),'b');
%     plot(time,e1_old(:,1),'--g');
    
    xlabel('t [s]')
    ylabel('e(t) [m]')
    xlim([tplotmin tplotmax])
    ylim([-0.7*s0 0.7*s0])
    grid on
    legend('Combined signal','1st-order term','2nd-order term')

%% write the output file with all the 1st and 2nd order time series
time1=time';
S={'Time','u','v','eta','u_2','v_2','eta_2','e1','e2','e'};
OUT=[time1 u v eta u_2 v_2 eta_2 e1 e2 e];
txt=sprintf('%s\t',S{:});
txt(end)='';
dlmwrite(strcat(name,'_theor.txt'),txt,'');
dlmwrite(strcat(name,'_theor.txt'),OUT,'-append','delimiter','\t','precision','%4.4f');

%% application of ramp to avoid discontinuities at the beginning and end of the timef series

e(1:rf,1)=e(1:rf,1).*quadd';
%     e(length(t_p)-rf+1:length(t_p),1)=e(length(t_p)-rf+1:length(t_p),1).*(1-quadd');

% generation of the timef series for CIEM wavemaker for NUMERICAL MODELLING
% (DualSPHysics)

% cd ..
cd CONTROL
mkdir(namefolder);
cd(namefolder)

MOVX=e;
MOVY=-MOVX*tan(pi/6);
MOV = [time1 MOVX MOVY];
firstline=[0,0,0];
MOV =[firstline; MOV];

filename_save = strcat('CIEM_SPH_',num2str(name),'.dat');
fid = fopen(filename_save, 'w');
fprintf(fid, '%7.4f\t %7.4f\t %7.4f\n', MOV');
fclose(fid);

% generation of the timef series for CIEM wavemaker for EXPERIMENTS (CIEM)

filename_save = strcat('CIEM_mov_PdP_',num2str(name),'.dat');
fid = fopen(filename_save, 'w');
fprintf(fid, '%7.4f\n', MOVX');
fclose(fid);

% cd ..
cd ..
