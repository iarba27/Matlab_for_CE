//###############################################################################################################
//-Calculo de velocidad de fluido en x,z.
//###############################################################################################################
const double f=1./WavePeriod;      //-Wave frequency.
const double kl=TWOPI/WaveLength;  //-Wave number.
const double sinhkld=sinh(kl*Depth);
const double ce=Gravity*WavePeriod/TWOPI*tanh(kl*Depth); //-1st and 2nd order wave celerity.
//-Matlab viejo.
vel.x= WaveHeight/2.*TWOPI*f*cosh(kl*(Depth+z))*(cos(TWOPI*f*time-kl*x+InitialPhase))/sinhkld;
vel.z=-WaveHeight/2.*TWOPI*f*sinh(kl*(Depth+z))*(sin(TWOPI*f*time-kl*x+InitialPhase))/sinhkld;
if(order2)vel.x+=3./4.*(PI*WaveHeight/WaveLength)*(PI*WaveHeight/WaveLength)*ce*cosh(2.*kl*(Depth+z))*(cos(4.*PI*f*time-2.*kl*x+2.*InitialPhase))/(sinhkld*sinhkld*sinhkld*sinhkld);
if(order2)vel.z-=3./4.*(PI*WaveHeight/WaveLength)*(PI*WaveHeight/WaveLength)*ce*sinh(2.*kl*(Depth+z))*(sin(4.*PI*f*time-2.*kl*x+2.*InitialPhase))/(sinhkld*sinhkld*sinhkld*sinhkld);
//-Matlab nuevo.
//vel.x= WaveHeight/2.*TWOPI*f*cosh(kl*(Depth+z))*(sin(TWOPI*f*time-kl*x+InitialPhase))/sinhkld;
//vel.z=-WaveHeight/2.*TWOPI*f*sinh(kl*(Depth+z))*(cos(TWOPI*f*time-kl*x+InitialPhase))/sinhkld;
//if(order2)vel.x+=3./4.*(PI*WaveHeight/WaveLength)*(PI*WaveHeight/WaveLength)*ce*cosh(2.*kl*(Depth+z))*(cos(4.*PI*f*time-2.*kl*x+2.*InitialPhase+PI))/(sinhkld*sinhkld*sinhkld*sinhkld);
//if(order2)vel.z-=3./4.*(PI*WaveHeight/WaveLength)*(PI*WaveHeight/WaveLength)*ce*sinh(2.*kl*(Depth+z))*(sin(4.*PI*f*time-2.*kl*x+2.*InitialPhase+PI))/(sinhkld*sinhkld*sinhkld*sinhkld);


//Codigo en ficheros de matlab:

//old   u_inc(i,1)=H/2*2*pi*f*cosh(kL*(d+z))*(cos(2*pi*f*time(i)-kL*x+ph0))/sinh(kL*d)+3/4*(pi*H/L)^2*C*cosh(2*kL*(d+z))*(cos(4*pi*f*time(i)-2*kL*x+2*ph0   ))/sinh(kL*d)^4;
//new   u_inc(i,1)=H/2*2*pi*f*cosh(kL*(d+z))*(sin(2*pi*f*time(i)-kL*x+ph0))/sinh(kL*d)+3/4*(pi*H/L)^2*C*cosh(2*kL*(d+z))*(cos(4*pi*f*time(i)-2*kL*x+2*ph0+pi))/sinh(kL*d)^4;
//diff----------------------------------------xxx--------------------------------------------------------------------------------------------------------xxx----------------

//old   v_inc(i,1)=-H/2*2*pi*f*sinh(kL*(d+z))*(sin(2*pi*f*time(i)-kL*x+ph0))/sinh(kL*d)-3/4*(pi*H/L)^2*C*sinh(2*kL*(d+z))*(sin(4*pi*f*time(i)-2*kL*x+2*ph0   ))/sinh(kL*d)^4;
//new   v_inc(i,1)=-H/2*2*pi*f*sinh(kL*(d+z))*(cos(2*pi*f*time(i)-kL*x+ph0))/sinh(kL*d)-3/4*(pi*H/L)^2*C*sinh(2*kL*(d+z))*(sin(4*pi*f*time(i)-2*kL*x+2*ph0+pi))/sinh(kL*d)^4;
//diff-----------------------------------------xxx--------------------------------------------------------------------------------------------------------xxx----------------

//-En el ultimo matlab es igual al primero...
//new2  u_inc(i,1)=H/2*2*pi*f*cosh(kL*(d+z))*(cos(2*pi*f*time(i)-kL*x+ph0))/sinh(kL*d)+3/4*(pi*H/L)^2*C*cosh(2*kL*(d+z))*(cos(4*pi*f*time(i)-2*kL*x+2*ph0))/sinh(kL*d)^4;
//new2  v_inc(i,1)=-H/2*2*pi*f*sinh(kL*(d+z))*(sin(2*pi*f*time(i)-kL*x+ph0))/sinh(kL*d)-3/4*(pi*H/L)^2*C*sinh(2*kL*(d+z))*(sin(4*pi*f*time(i)-2*kL*x+2*ph0))/sinh(kL*d)^4;
