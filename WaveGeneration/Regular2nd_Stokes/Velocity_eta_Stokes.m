%This script compares the SPH measured velocities with theoretical ones for
%Regulas waves

%%
%set the input parameters for SPH and CIEM

clc
clear all
close all

f_out = 40;        %sample SPH frequency (number of fps)
nADV = 1;                  %number of ADV in SPH
tmax = 30;                %duration SPH run
name='Dike004';
%Target wave 

H=0.1; %wave height
T=1.3; %wave period
d=0.266; %water depth
L=L_lin(T,d);
kL=2*pi/L; %wave number
f=1/T; %frequency
x=2; %x coordinates ADV
z=-0.15; %z coordinate ADV
g=9.81; %gravity acceleration

CR=0.0; %reflection coeffcient
   
%-------------------------------------------------------------------%


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
        % 2nd order Stokes incident
        eta_inc(i,1)=H/2*cos(2*pi*f*time(i)-kL*x)+(kL*H^2)/16*((3*((coth(kL*d))^3))-coth(kL*d))*cos(4*pi*f*time(i)-2*kL*x);   
%         eta_incSPM(i,1)=H/2*cos(2*pi*f*time(i)-kL*x)+(kL*H^2)/16*(cosh(kL*d)/(sinh(kL*d)^3)*(2+cosh(2*kL*d)))*cos(4*pi*f*time(i)-2*kL*x);   
        u_inc(i,1)=H/2*2*pi*f*cosh(kL*(d+z))*(cos(2*pi*f*time(i)-kL*x))/sinh(kL*d)+3/4*(pi*H/L)^2*C*cosh(2*kL*(d+z))*(cos(4*pi*f*time(i)-2*kL*x))/sinh(kL*d)^4;
        v_inc(i,1)=-H/2*2*pi*f*sinh(kL*(d+z))*(sin(2*pi*f*time(i)-kL*x))/sinh(kL*d)-3/4*(pi*H/L)^2*C*sinh(2*kL*(d+z))*(sin(4*pi*f*time(i)-2*kL*x))/sinh(kL*d)^4;
        
        % only second component Stokes incident for eta
        eta_inc2(i,1)=+(kL*H^2)/16*((3*((coth(kL*d))^3))-coth(kL*d))*cos(4*pi*f*time(i)-2*kL*x);   
        
        %2nd order Stokes reflected 
        eta_refl(i,1)=CR*(H/2*cos(2*pi*f*time(i)+kL*x)+(pi*H^2/(8*L))*((3*(coth(kL*d)^3)-coth(kL*d)))*cos(4*pi*f*time(i)+2*kL*x));   
        u_refl(i,1)=CR*(H/2*2*pi*f*cosh(kL*(d+z))*(cos(2*pi*f*time(i)+kL*x))/sinh(kL*d)+3/4*(pi*H/L)^2*C*cosh(2*kL*(d+z))*(cos(4*pi*f*time(i)+2*kL*x))/sinh(kL*d)^4);
        v_refl(i,1)=CR*(H/2*2*pi*f*sinh(kL*(d+z))*(sin(2*pi*f*time(i)+kL*x))/sinh(kL*d)+3/4*(pi*H/L)^2*C*sinh(2*kL*(d+z))*(sin(4*pi*f*time(i)+2*kL*x))/sinh(kL*d)^4);
    
        %2nd order Stokes total
        eta_2(i,1)=eta_inc(i,1)+eta_refl(i,1);
        u_2(i,1)=u_inc(i,1)-u_refl(i,1);
        v_2(i,1)=v_inc(i,1)-v_refl(i,1);
        
        % 1st order (Airy)
        eta(i,1)=H/2*cos(2*pi*f*time(i)-kL*x)+CR*H/2*cos(2*pi*f*time(i)+kL*x);  
        u(i,1)=H/2*2*pi*f*cosh(kL*(d+z))*(cos(2*pi*f*time(i)-kL*x))/sinh(kL*d)-CR*(H/2)*2*pi*f*cosh(kL*(d+z))*cos(2*pi*f*time(i)+kL*x)/sinh(kL*d);
        v(i,1)=-H/2*2*pi*f*sinh(kL*(d+z))*(sin(2*pi*f*time(i)-kL*x))/sinh(kL*d)-CR*(H/2)*2*pi*f*sinh(kL*(d+z))*sin(2*pi*f*time(i)+kL*x)/sinh(kL*d); 
end

%Figure


figure(1)
figsize1 = [100 100 700 600];
figure(1); clf(figure(1)); set(gcf, 'color', 'white','Position', figsize1);


    subplot(3,1,1)
    plot(time,u(:,1),'b'); hold on 
    plot(time,u_2(:,1),'r')    
    xlabel('t [s]')
    ylabel('u [m/s]')
    xlim([0 tmax])
    legend('Theoretical (Airy)','Theoretical (2nd order Stokes)' )
    ylim([-0.4 0.4])
    
    subplot(3,1,2)
    plot(time,v(:,1),'b');hold on 
    plot(time,v_2(:,1),'r')  
    xlabel('t [s]')
    ylabel('v [m/s]')
    xlim([0 tmax])
    ylim([-0.4 0.4])
    
    subplot(3,1,3)
    plot(time,eta(:,1),'b');hold on 
    plot(time,eta_2(:,1),'r')
%     plot(time,eta_incSPM(:,1),'--k')
%     plot(time,eta_inc2(:,1),'--m') 
    xlabel('t [s]')
    ylabel('\eta [m]')
    xlim([0 tmax])
    ylim([-1.5*H 1.5*H])
    grid on

%% write the output file with all the 1st and 2nd order time series
time1=time';
L={'Time','u','v','eta','u_2','v_2','eta_2'};
OUT=[time1 u v eta u_2 v_2 eta_2];
txt=sprintf('%s\t',L{:});
txt(end)='';
dlmwrite(strcat(name,'_theor.txt'),txt,'');
dlmwrite(strcat(name,'_theor.txt'),OUT,'-append','delimiter','\t','precision','%4.4f');

% save(strcat(name,'_R_time.txt'),'time','-ASCII');
% save(strcat(name,'_R_vel_1st_x.txt'),'u','-ASCII');
% save(strcat(name,'_R_vel_1st_z.txt'),'v','-ASCII');
% save(strcat(name,'_R_eta_1st.txt'),'eta','-ASCII');
% save(strcat(name,'_R_vel_2nd_x.txt'),'u_2','-ASCII');
% save(strcat(name,'_R_vel_2nd_z.txt'),'v_2','-ASCII');
% save(strcat(name,'_R_eta_2nd.txt'),'eta_2','-ASCII');
