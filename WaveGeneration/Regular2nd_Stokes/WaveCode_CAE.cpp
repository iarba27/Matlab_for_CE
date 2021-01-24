//###############################################################################################################
//-Calculo de posicion del piston (version previa).
//###############################################################################################################
WaveLength=CalcWaveLength(Depth,WavePeriod);
Bisel=CalcBisel(Depth,FixedDepth,WaveLength);
Stroke=WaveHeight/Bisel;
Amplitude=Stroke/2;
phase=InitialPhase+(PI+PI)*t/WavePeriod;
pos=Amplitude*sin(phase);    //-1st order term.

//==============================================================================
// Calcula Bisel value (B=H/So). 1st order Biesel function.
//==============================================================================
double JWavePaddle::CalcBisel(double depth,double fixeddepth,double wavelength){
  double k=(PI*2/wavelength); //-Wave number.
  double kh=k*depth;
  double kh0=k*fixeddepth;                                                                     
  double b=(2.*sinh(kh)*sinh(kh)-2.*sinh(kh0)*sinh(kh))/(sinh(kh)*cosh(kh)+kh);
  return(b);
}


//###############################################################################################################
//-Calculo de posicion del piston (formula nueva).
//###############################################################################################################
const double phase=InitialPhase+(PI+PI)*t/WavePeriod;
const double k=(PI*2/WaveLength); //-Wave number.
const double kh=k*Depth;
const double m1=4.*(sinh(kh)*sinh(kh))/(sinh(2.*kh)+2.*kh);       //-1st order Biesel function S/H.
const double o2=(WaveHeight*WaveHeight)/(32.*Depth)*(3.*cosh(kh)/(sinh(kh)*sinh(kh)*sinh(kh))-2./m1); //-2nd order component.
//const double e1=0.5*WaveHeight/m1*cos(phase);                   //-1st order term.
const double e1=0.5*WaveHeight/m1*sin(phase);                     //-1st order term (usando sin()).
//const double e2=o2*cos(4.*PI*t/WavePeriod+PI/2.+2.*InitialPhase); //-2nd order term (Madsen, 1971).
const double e2=o2*sin(4.*PI*t/WavePeriod+2.*InitialPhase); //-2nd order term (Madsen, 1971).                                    //CORRADO
pos=e1+e2;


//###############################################################################################################
//-Calculo de posicion del piston (lo que realmente quedo implementado).
//###############################################################################################################
//-Valores precalculados.
WaveLength=CalcWaveLength(Depth,WavePeriod);
Bisel=CalcBisel(Depth,FixedDepth,WaveLength);
Stroke=WaveHeight/Bisel;                     
Amplitude=Stroke/2;
Cte2ndOrder=CalcCte2ndOrder(Depth,FixedDepth,WaveLength,WaveHeight)                          
//-Codigo en CalcPosition(double t,bool order2)
phase=InitialPhase+(PI+PI)*t/WavePeriod;
pos=Amplitude*sin(phase);    //-1st order term.
if(order2)pos+=Cte2ndOrder*cos(4.*PI*t/WavePeriod + PI/2. + 2.*InitialPhase); //-2nd order term (Madsen, 1971).

//==============================================================================
// Calcula Cte2ndOrder value. 2nd order component.
//==============================================================================
double CalcCte2ndOrder(double depth,double fixeddepth,double wavelength,double waveheight){
  double k=(PI*2/wavelength); //-Wave number.
  double kh=k*depth;
  double bisel=CalcBisel(depth,fixeddepth,wavelength);
  double o2=(waveheight*waveheight)/(32.*depth)*(3.*cosh(kh)/(sinh(kh)*sinh(kh)*sinh(kh))-2./bisel);
  return(fixeddepth? 0.: o2);//-Formula no valida para fixeddepth!=0.
}







//###############################################################################################################
//-Calculo de elevacion (version previa).
//###############################################################################################################
const double kl=(PI+PI)/WaveLength;  //-Wave number.
const double phase=InitialPhase + (PI+PI)*t/WavePeriod - Kl*x;
ele=WaveHeight/2*cos(phase);                                                                                      


//###############################################################################################################
//-Calculo de posicion del piston (lo que quedo implementado).
//###############################################################################################################
//-Valores precalculados.
const double kl=(PI+PI)/WaveLength;  //-Wave number.
CteEle2nd=CalcCteEle2nd(Depth,FixedDepth,WaveLength,WaveHeight)                          
//-Codigo en CalcElevation(double t,double x,bool order2)
const double phase=InitialPhase + (PI+PI)*t/WavePeriod - kl*x;
ele=WaveHeight/2*cos(phase);                                                                                      //CORRADO
if(order2)ele+=CteEle2nd * cos(4.*PI*t/WavePeriod - 2.*kl*x + 2.*InitialPhase);                             //CORRADO

//==============================================================================
// Calcula CteEle2nd value. 2nd order component for elevation.
//==============================================================================
double CalcCteEle2nd(double depth,double fixeddepth,double wavelength,double waveheight){
  double kl=(PI+PI)/wavelength;  //-Wave number.
  double coth=fmath::coth(kl*depth);
  double cte=(kl*waveheight*waveheight)/16 * ((3.*coth*coth*coth)-coth);   
  return(fixeddepth? 0.: cte);//-Formula no valida para fixeddepth!=0.
}





