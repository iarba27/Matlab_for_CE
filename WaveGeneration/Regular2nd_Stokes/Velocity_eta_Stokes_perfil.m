%This script compares the SPH measured velocities with theoretical ones for
%Regulas waves

%%
%set the input parameters for SPH and CIEM

clc
clear all
close all

f_out = 10;        %sample SPH frequency (number of fps)
% nADV = 1;                  %number of ADV in SPH
tmax = 40;                %duration SPH run
%Target wave 

H=0.05; %wave height
T=1.25; %wave period
d=0.77; %water depth
L=L_lin(T,d);
kL=2*pi/L; %wave number
f=1/T; %frequency
x=6.5; %x coordinates ADV
z=[-d:d/11:0];
% z(n)=-0.15; %z(n) coordinate ADV
g=9.81; %gravity acceleration

   
%-------------------------------------------------------------------%


%% set timeline for SPH
row=f_out*tmax;

for i=1:row;
    time(i)=i/f_out;
end

%% theoretical velcoities

%1st and 2nd order wave celerity
C=g*T/(2*pi)*tanh(kL*d); 



for n=1:size(z,2)
    for i=1:row;
        % 2nd order Stokes incident
        eta_inc(i,n)=H/2*cos(2*pi*f*time(i)-kL*x)+(kL*H^2)/16*((3*((coth(kL*d))^3))-coth(kL*d))*cos(4*pi*f*time(i)-2*kL*x);   
%         eta_incSPM(i,1)=H/2*cos(2*pi*f*time(i)-kL*x)+(kL*H^2)/16*(cosh(kL*d)/(sinh(kL*d)^3)*(2+cosh(2*kL*d)))*cos(4*pi*f*time(i)-2*kL*x);   
        u_inc(i,n)=H/2*2*pi*f*cosh(kL*(d+z(n)))*(cos(2*pi*f*time(i)-kL*x))/sinh(kL*d)+3/4*(pi*H/L)^2*C*cosh(2*kL*(d+z(n)))*(cos(4*pi*f*time(i)-2*kL*x))/sinh(kL*d)^4;
        v_inc(i,n)=-H/2*2*pi*f*sinh(kL*(d+z(n)))*(sin(2*pi*f*time(i)-kL*x))/sinh(kL*d)-3/4*(pi*H/L)^2*C*sinh(2*kL*(d+z(n)))*(sin(4*pi*f*time(i)-2*kL*x))/sinh(kL*d)^4;

      
        % 1st order (Airy)
        eta(i,n)=H/2*cos(2*pi*f*time(i)-kL*x);
        u(i,n)=H/2*2*pi*f*cosh(kL*(d+z(n)))*(cos(2*pi*f*time(i)-kL*x))/sinh(kL*d);
                u_uu(i,n)=H/2*9.81*T/L*cosh(kL*(d+z(n)))*(cos(2*pi*f*time(i)-kL*x))/cosh(kL*d);
        v(i,n)=-H/2*2*pi*f*sinh(kL*(d+z(n)))*(sin(2*pi*f*time(i)-kL*x))/sinh(kL*d);
    end
    
        u_max(1,n)=max(u_inc(:,n));
        u_min(1,n)=min(u_inc(:,n));
        v_max(1,n)=max(v_inc(:,n));
        v_min(1,n)=min(v_inc(:,n));
        
end

SAVE=[z' u_max' u_min' v_max' v_min']; 
save('vel_maxmin.txt','SAVE','-ASCII'); 
%% Figure

cont=round(size(z,2)/2);

figure(1)
figsize1 = [100 100 700 600];
figure(1); clf(figure(1)); set(gcf, 'color', 'white','Position', figsize1);
for n=1:size(z,2)
    subplot(cont,2,n)
%     plot(time,u(:,size(z,2)-n+1),'b'); hold on 
    plot(time,u_inc(:,size(z,2)-n+1),'b')    
    xlabel('t [s]')
    ylabel('u [m/s]')
    xlim([0 tmax])
    ylim([(u_min(1,(size(z,2)-n+1))*1.5) (u_max(1,(size(z,2)-n+1))*1.5)])
    a=num2str(z(size(z,2)-n+1));
    name=(strcat('z=',a,' [m]'));
    title(name);
end
    
    
figure(2)
figsize2 = [700 100 700 600];
figure(2); clf(figure(2)); set(gcf, 'color', 'white','Position', figsize2);
    for n=1:size(z,2)
    subplot(cont,2,n)
%     plot(time,u(:,size(z,2)-n+1),'b'); hold on 
    plot(time,v_inc(:,size(z,2)-n+1),'m')    
    xlabel('t [s]')
    ylabel('v [m/s]')
    xlim([0 tmax])
      a=num2str(z(size(z,2)-n+1));
    name=(strcat('z=',a,' [m]'));
    title(name);
%     ylim([(v_min(1,(size(z,2)-n+1))*1.5) (v_max(1,(size(z,2)-n+1))*1.5)])
  
end

 



