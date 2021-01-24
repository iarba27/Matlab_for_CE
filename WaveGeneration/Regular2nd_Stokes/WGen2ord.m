%This script compares the SPH measured velocities with theoretical ones for
%Regulas waves and Calculate the time series of a PISTON-type wavemaker for
%1st and 2nd order wave generation theory (Madsen approximation, 1971)

%%
%set the input parameters for SPH and CIEM

clc
clear all
close all

%% input 
CR=1; %reflection coeffcient

H=0.15; %wave height
T=2; %wave period
d=0.66; %water depth

L=L_lin(T,d);
kL=2*pi/L; %wave number
f=1/T; %frequency

x=L/4; %x coordinates ADV
z=-0.26; %z coordinate ADV
g=9.81; %gravity acceleration
ph0=0; %initial phase

f_out = 40;        %sample SPH frequency (number of fps)
nADV = 1;                  %number of ADV in SPH
tmax = 50;                %duration SPH run
name='2ndWG';

%%
%time min and max for plots
tplotmin=0;
tplotmax=4*T;



%% for PISTON wave generation 

m1=4*(sinh(kL*d)^2)/(sinh(2*kL*d)+2*kL*d); %1st order Biesel function S/H
s0=H/m1; %piston max Stroke 

o2=(H^2)/(32*d)*(3*cosh(kL*d)/(sinh(kL*d)^3)-2/m1); %2nd order component

% check for Madsen (1971) theory Madcrit<Madlim

Madlim=8*pi*pi/3
Madcrit=H*L*L/d^3
  
%-------------------------------------------------------------------%

%% set timeline for SPH
row=f_out*tmax;

for i=1:row;
    time(i)=i/f_out;
end

%% theoretical velcoities

%1st and 2nd order wave celerity
C=g*T/(2*pi)*tanh(kL*d); 


for i=1:row;
    
        %% piston displacement 
        e1(i,1)=0.5*s0*sin(2*pi*f*time(i)+ph0); %1st term
        e2(i,1)=o2*sin(4*pi*f*time(i)+2*ph0); %2nd order term (Madsen, 1971)
        e(i,1)=e1(i,1)+e2(i,1); %total signal 
        
    
        %% incident
        % 2nd order Stokes incident
        eta_inc(i,1)=H/2*cos(2*pi*f*time(i)-kL*x+ph0)+(kL*H^2)/16*((3*((coth(kL*d))^3))-coth(kL*d))*cos(4*pi*f*time(i)-2*kL*x+2*ph0);   
        u_inc(i,1)=H/2*2*pi*f*cosh(kL*(d+z))*(cos(2*pi*f*time(i)-kL*x+ph0))/sinh(kL*d)+3/4*(pi*H/L)^2*C*cosh(2*kL*(d+z))*(cos(4*pi*f*time(i)-2*kL*x+2*ph0))/sinh(kL*d)^4;
        v_inc(i,1)=-H/2*2*pi*f*sinh(kL*(d+z))*(sin(2*pi*f*time(i)-kL*x+ph0))/sinh(kL*d)-3/4*(pi*H/L)^2*C*sinh(2*kL*(d+z))*(sin(4*pi*f*time(i)-2*kL*x+2*ph0))/sinh(kL*d)^4;
        
        % only second component Stokes incident for eta
        eta_inc2(i,1)=+(kL*H^2)/16*((3*((coth(kL*d))^3))-coth(kL*d))*cos(4*pi*f*time(i)-2*kL*x+2*ph0); 
                
        %1st order Airy reflected 
        eta_refl1(i,1)=CR*(H/2*cos(2*pi*f*time(i)+kL*x));   
        u_refl1(i,1)=CR*(H/2*2*pi*f*cosh(kL*(d+z))*(cos(2*pi*f*time(i)+kL*x))/sinh(kL*d));
        v_refl1(i,1)=CR*(H/2*2*pi*f*sinh(kL*(d+z))*(sin(2*pi*f*time(i)+kL*x))/sinh(kL*d));
        
        %2nd order Stokes reflected 
        eta_refl(i,1)=CR*(H/2*cos(2*pi*f*time(i)+kL*x)+(pi*H^2/(8*L))*((3*(coth(kL*d)^3)-coth(kL*d)))*cos(4*pi*f*time(i)+2*kL*x));   
        u_refl(i,1)=CR*(H/2*2*pi*f*cosh(kL*(d+z))*(cos(2*pi*f*time(i)+kL*x))/sinh(kL*d)+3/4*(pi*H/L)^2*C*cosh(2*kL*(d+z))*(cos(4*pi*f*time(i)+2*kL*x))/sinh(kL*d)^4);
        v_refl(i,1)=CR*(H/2*2*pi*f*sinh(kL*(d+z))*(sin(2*pi*f*time(i)+kL*x))/sinh(kL*d)+3/4*(pi*H/L)^2*C*sinh(2*kL*(d+z))*(sin(4*pi*f*time(i)+2*kL*x))/sinh(kL*d)^4);
    
        %2nd order Stokes total
        eta_2(i,1)=eta_inc(i,1)+eta_refl(i,1);
        u_2(i,1)=u_inc(i,1)-u_refl(i,1);
        v_2(i,1)=v_inc(i,1)-v_refl(i,1);
        
               % 1st order (Airy)
        eta(i,1)=H/2*cos(2*pi*f*time(i)-kL*x+ph0)+eta_refl1(i,1);
        u(i,1)=H/2*2*pi*f*cosh(kL*(d+z))*(cos(2*pi*f*time(i)-kL*x+ph0))/sinh(kL*d)-u_refl1(i,1);
        v(i,1)=-H/2*2*pi*f*sinh(kL*(d+z))*(sin(2*pi*f*time(i)-kL*x+ph0))/sinh(kL*d)-v_refl1(i,1);
end

umax=max(u_inc);
vmax=max(v_inc);

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
    ylim([-2*umax 2*umax])
    grid on
    
    subplot(2,2,2)
    plot(time,v(:,1),'b');hold on 
    plot(time,v_2(:,1),'r')  
    xlabel('t [s]')
    ylabel('v [m/s]')
    xlim([tplotmin tplotmax])
    ylim([-2*vmax 2*vmax])
    legend('Theoretical (Airy)','Theoretical (2nd order Stokes)' )
    grid on
    
    subplot(2,2,3)
    plot(time,eta(:,1),'b');hold on 
    plot(time,eta_2(:,1),'r')
    xlabel('t [s]')
    ylabel('\eta [m]')
    xlim([tplotmin tplotmax])
    ylim([-1*H*(1.1+CR) 1*H*(1.1+CR)])
    grid on
    legend('Theoretical (Airy)','Theoretical (2nd order Stokes)' )
    
    subplot(2,2,4)
    plot(time,e1(:,1),'b');      hold on
    plot(time,e(:,1),'r');
    plot(time,e2(:,1),'--k');  
    xlabel('t [s]')
    ylabel('e(t) [m]')
    xlim([tplotmin tplotmax])
    ylim([-0.7*s0 0.7*s0])
    grid on
    legend('1st-order solution','2nd-order solution','2nd-order term')

%% write the output file with all the 1st and 2nd order time series
time1=time';
S={'Time','u','v','eta','u_2','v_2','eta_2','e1','e2','e'};
OUT=[time1 u v eta u_2 v_2 eta_2 e1 e2 e];
txt=sprintf('%s\t',S{:});
txt(end)='';
if x==L/4
    name2='_node';
elseif x==L/2
    name2='_antinode';
else
    name2='';
end
dlmwrite(strcat(name,'_CR',num2str(CR),name2,'_theor.txt'),txt,'');
dlmwrite(strcat(name,'_CR',num2str(CR),name2,'_theor.txt'),OUT,'-append','delimiter','\t','precision','%4.4f');


