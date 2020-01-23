
#include "wrapper.h"

void GLoBESCalculator::SetDepth(double dep){ // Sets depth of the detector
     detDepth = dep;
     innerR = Rearth - detDepth;}

void GLoBESCalculator::SetAtmoHeight(double atmh){ // Sets atmospheric height
     atmoHeight = atmh;
     outerR = Rearth + atmoHeight; }

void GLoBESCalculator::SetEarthRadius(double erad){ // Sets radius  of the Earth (in vacuum)
     Rearth = erad;
     innerR = Rearth - detDepth;
     outerR = Rearth + atmoHeight; }

double GLoBESCalculator::GetBaselineInExperiment(int exp){ // return baseline
  return glbGetBaselineInExperiment(exp);}

void GLoBESCalculator::InitSteriles(int mode){
  /*
  Sterile mode 0: no sterile neutrinos
  Sterile mode 1: alpha mixing 
  Sterile mode 2: 4-3 and 4-2 mixing
  */
  sterileMode = mode;
  std::cout<<"Selecting sterile neutrino mode "<<sterileMode<<std::endl;
  switch (sterileMode){
    case 0: 
      {std::cout<<"No sterile neutrinos! "<<std::endl; 
      int rotation_order[][2] = {{2,3}, {1,3}, {1,2}};
      int phase_order[] =       { -1  ,   0  ,  -1  };      
      int res = snu_init_probability_engine(3, rotation_order, phase_order);  
      std::cout<<"Return code: "<<res<<std::endl;   
      break;}
    case 1: 
      {std::cout<<"The alpha mixing (mass mixing ) scenario is ON!"<<std::endl;
      int rotation_order[][2] = {{2,3}, {1,3}, {1,2}, {3,4}, {2,4}, {1,4} };
      int phase_order[] =       { -1  ,   0  ,   -1 ,  -1  ,   -1 ,   -1  };      
      int res = snu_init_probability_engine(4, rotation_order, phase_order);  
      std::cout<<"Return code: "<<res<<std::endl;    
      break; }
    case 2: 
      {std::cout<<"4-3 and 4-2 mixing is ON"<<std::endl;
      int rotation_order[][2] = {{3,4}, {2,4},{1,4},{2,3}, {1,3}, {1,2} };
      int phase_order[] =       {  -1,   -1  ,  -1 , -1,    0 ,  -1   };      
      int res = snu_init_probability_engine(4, rotation_order, phase_order);  
      std::cout<<"Return code: "<<res<<std::endl;           
      break;}
    default: sterileMode = -1;
             std::cout<<"ERROR !!! Wrong mode ! "<<sterileMode<<" is set! "<<std::endl;
             break;
  }
  if (sterileMode != -1 & sterileMode != 0){
      std::cout<<"Registering SNU module!\n";
      int n_osc_pars = 0;
      glbRegisterProbabilityEngine(6*4*4 - 4,
                               &snu_probability_matrix,
                               &snu_set_oscillation_parameters,
                               &snu_get_oscillation_parameters,
                               NULL);
  } else {
     if (sterileMode == 0){
     std::cout<<"Registering SNU module without steriles!\n";
     glbRegisterProbabilityEngine(6*3*3 - 3,
                               &snu_probability_matrix,
                               &snu_set_oscillation_parameters,
                               &snu_get_oscillation_parameters,
                               NULL);}
  }
  ApplyParameters();
  return;
}
int GLoBESCalculator::GetNumOfOscParams() {
  /* returns number of oscillation parameters (if SNU used)*/
  return glbGetNumOfOscParams();
}
void  GLoBESCalculator::PrintParametersOrder(){
  for (int i =0; i<12; i++){
    std::cout<<i<<" "<<snu_param_strings[i]<<std::endl;
  }
}
GLoBESCalculator::GLoBESCalculator(std::string init) {

  std::cout<<"Creating calculator"<<std::endl;
  
  glbInit(const_cast<char*> (init.c_str()) );
  std::cout<<"GLoBES: glbInit run successfully"<<std::endl;
  glbInitExperiment("Dummy_Experiment.glb" ,&glb_experiment_list[0], &glb_num_of_exps);
  std::cout<<"GLoBES wrapper: initialized experiment"<<std::endl; 
  detDepth = 2.0; atmoHeight = 20.0; 
  Rearth = 6371.0;
  innerR = 6371.0 - detDepth;
  outerR = 6371.0 + atmoHeight;
  PREM_layers = 0; 
  sterileMode = -1;
  SetParameters(0.0,0.0,0.0,0.0,0.0,0.0);
}

void GLoBESCalculator::SetParameters(
           double theta12, double theta13, double theta23, 
           double deltacp, double dm21, double dm31 ) {
  osc_params[0] = theta12;  osc_params[1] = theta13; 
  osc_params[2] = theta23;  osc_params[3] = deltacp;
  osc_params[4] = dm21;  osc_params[5] = dm31; 
  ApplyParameters();
  return;
}
/*
  Puts parameters in globes internal memory
*/
void GLoBESCalculator::ApplyParameters(){
 if (sterileMode == -1) {
    glb_params values = glbAllocParams();
    glbDefineParams(	values, 
	 				osc_params[0], osc_params[1], osc_params[2], osc_params[3],
					osc_params[4], osc_params[5] );
    glbSetOscillationParameters(values);
    glbFreeParams(values); }
  if (sterileMode == 0) {
    glb_params values = glbAllocParams();
    for (int i = 0; i < 6; i++){glbSetOscParams(values, osc_params[i], i);}
    for (int i = 6 ; i < 6*3*3 - 3; i++){glbSetOscParams(values, 0.0, i); }
    
    //glbDefineParams(	values, 
	// 				osc_params[0], osc_params[1], osc_params[2], osc_params[3],
	//				osc_params[4], osc_params[5] );
    glbSetOscillationParameters(values);
    glbFreeParams(values); }
  if (sterileMode ==1 | sterileMode ==2 ){
    for (int i = 0; i< 12; i++){   
        glb_params true_values = glbAllocParams();
        for (int i = 0; i < 12; i++){glbSetOscParams(true_values, osc_params[i], i);}
        for (int i = 12 ; i < 6*4*4 - 4; i++){glbSetOscParams(true_values, 0.0, i); }// Need to fill and NSI parameters (0 since all interactions are standard for us) 
        glbSetOscillationParameters(true_values);
        glbFreeParams(true_values);
    }
  }
  return;
}
/*
  Sets parameters for all neutrinos (standard + steriles)
*/
void GLoBESCalculator::SetParametersArr(boost::python::list& params){
  if (len(params )> 12 | len(params) < 6)
    {std::cout<<"Too small or too big number of parameteres: "<<len(params )<<" !!!\n";
     return;}
  for (int i = 0 ; i < len(params); i++){ osc_params[i] = boost::python::extract<double>(params[i]) ; }
  ApplyParameters();
} 

void GLoBESCalculator::PrintParameters() {
    std::cout<<"Current parameters:"<<std::endl;
    glb_params values = glbAllocParams();
    glbGetOscillationParameters(values);
    if (sterileMode == 0){
       printf("  \t th12 \t th13 \t th23 \t dcp \t dm21 \t dm31 \n\t");
       for (int k = 0; k < 6; k++){
         std::cout<<glbGetOscParams(values, k)<<"\t";
       }
       std::cout<<std::endl;
    }
    if (sterileMode == 1 | sterileMode == 2){
      std::cout<<"\t";
      for (int i = 0 ; i<12; i++) { std::cout<<snu_param_strings[i]<<"\t";}
      std::cout<<"\n\t"; 
      for (int i = 0 ; i<12; i++) { std::cout<<glbGetOscParams(values, i)<<"\t";} 
      std::cout<<std::endl;   
    }
    glbFreeParams(values);
    return;
}
/*
	Calculator for baseline
*/
double GLoBESCalculator::CalcPropDist(double costh){
  return (innerR)*(-costh) + sqrt( pow(outerR,2) - pow(innerR,2)*(1.0 - pow(costh,2)));
}
/*
	Calculator for vacuum probabilities
*/
double GLoBESCalculator::VacuumProbability(int m,int l ,int panti,double E,double costh){
  return glbVacuumProbability(m,l,panti,E,CalcPropDist(costh));
}
/*
	Calculator for probabilities with matter
*/
double GLoBESCalculator::MatterProbability(int m,int l ,int panti,double E,double costh){
  CalcDensityLayers(costh, false);
  return glbProfileProbability(0,m,l,panti,E);
}



double GLoBESCalculator::MatterProbPDG(int ini, int fin, double E ,double costh){
   int m= 0, l = 0, pa = 0;
   switch (ini) {
     case  12: m = 1 ; break;
     case -12: m = 1 ; break;
     case  14: m = 2 ; break;
     case -14: m = 2 ; break;
     case  16: m = 3 ; break;
     case -16: m = 3 ; break;
     case  18: m = 4 ; break;
     case -18: m = 4 ; break;
     default: std::cout<<"Unknown initial neutrino! "<<std::endl;
   }
   switch (fin) {
     case  12: l = 1 ; break;
     case -12: l = 1 ; break;
     case  14: l = 2 ; break;
     case -14: l = 2 ; break;
     case  16: l = 3 ; break; 
     case -16: l = 3 ; break;
     case  18: l = 4 ; break; 
     case -18: l = 4 ; break;
     default: std::cout<<"Unknown initial neutrino! "<<std::endl;
   }
   if (ini > 0){pa =  1 ;}
   if (ini < 0){pa = -1 ;}
   //std::cout<<"I : "<<ini<<"("<<m<<")\t"<<" F : "<<fin<<"("<<l<<")\t PA: "<<pa<<std::endl;
   return MatterProbability(m, l, pa,E, costh);
}

/*
  Sets layers of the Earth 
*/
void GLoBESCalculator::SetEarthModel(boost::python::list& rads, boost::python::list& dens){
  if ( len(rads) != len(dens) | len(rads) ==0 | len(dens) ==0){
    std::cout<<"Densities and lengths should be of equal length and non-zero!"<<std::endl;
    std::cout<<"Sizes: "<<len(rads)<<"(radius) \t"<<len(dens)<<"(dens)\n";
    return ;
  }
  //std::cout<<"Sizes: "<<len(rads)<<"(radius) \t"<<len(dens)<<"(dens)\n";
  PREM_layers = len(rads);
  numInnerLayers = 0; numOuterLayers = 0;
  RadPREM = boost::python::extract<double>(rads[PREM_layers-1]);
  innerR = RadPREM - detDepth; outerR = RadPREM + atmoHeight;
  PREM_rad = std::vector<double>(PREM_layers);
  PREM_dens = std::vector<double>(PREM_layers);
  for (int i =0 ; i<PREM_layers; i++){
     PREM_rad[i]  = boost::python::extract<double>(rads[i]);
     PREM_dens[i] = boost::python::extract<double>(dens[i]);
  }
  for (int i = 0; i<PREM_layers; i++){
     if (PREM_rad[i] < innerR){numInnerLayers ++;}
     else {numOuterLayers ++;}
  }
  numInnerLayers++; // layer below the detector
  numOuterLayers++; // PREM does not have an atmosphere in
  // Part below the detector
  inRads = std::vector<double>(numInnerLayers);
  inDens = std::vector<double>(numInnerLayers);
  for (int i =0 ; i< numInnerLayers-1; i++){
    inRads[i] = PREM_rad[i];
    inDens[i] = PREM_dens[i];
  }
  inRads[numInnerLayers-1]  = innerR; 
  inDens[numInnerLayers-1]  = PREM_dens[numInnerLayers-1]; 
  // part above the detector + atmosphere
  outRads = std::vector<double>(numOuterLayers);  
  outDens = std::vector<double>(numOuterLayers);  
  for (int i =0; i<numOuterLayers-1; i ++){
     outRads[i] = PREM_rad[numInnerLayers-1 + i];
     outDens[i] = PREM_dens[numInnerLayers-1 + i];
  }
  outRads[numOuterLayers-1] = outerR;
  outDens[numOuterLayers-1] = 0.0; //atmosphere
  CalcPREMCosines();
  return;
}
/*
  Calculates viewing angles of different layers inside the Earth for further estimation of numper of layers
*/
int GLoBESCalculator::CalcPREMCosines(){
  PREM_cosines = std::vector<double>(numInnerLayers);
  for (int i =0; i < numInnerLayers; i++){
    PREM_cosines[i] =  - sqrt(1.0 - pow(inRads[i]/innerR, 2) ); 
  }
  /*
  PREM_cosines = std::vector<double>(PREM_layers);
  for (int i =0; i < PREM_layers; i++){
    PREM_cosines[i] = - sqrt( pow(RadPREM,2) - pow(PREM_rad[i],2) ) / (RadPREM);
  }*/
  return 0;
}
/*
  Prints the Earth matter model
*/
void GLoBESCalculator::PrintEarthModel(){
  std::cout<<"Earth model with "<<PREM_layers<<" layers is set (R = "<<RadPREM<<" km ). \n";
  std::cout<<"Depth: "<<detDepth<<" km. \n R below detector: "<<innerR<<" km ("<<numInnerLayers<< " layers). Earth + atmosphere: "<<outerR<<" km("<<numOuterLayers<<" layers above detector).\n";
  std::cout<<"\t\t\tInner layers:\n";
  for (int i =0 ; i < numInnerLayers;i++){ printf("%i\t%f\t%f\t%f\n", i, inRads[i], inDens[i],PREM_cosines[i]); }
  std::cout<<"\t\t\tOuter layers:\n";  
  for (int i =0 ; i < numOuterLayers;i++){ printf("%i\t%f\t%f\t\n", i, outRads[i], outDens[i]); }
  return;
}
/* 
  Prints profile for the last call. 
*/
void GLoBESCalculator::PrintDensityProfile(){ 
  size_t laysTemp = 0;
  double *lenTemp, *densTemp;
  glbGetProfileDataInExperiment(0, &laysTemp, &lenTemp, &densTemp);
  std::vector<double> layers_len(lenTemp, lenTemp+laysTemp);
  std::vector<double> layers_dens(densTemp, densTemp+laysTemp);
  for (int i = 0; i < laysTemp; i++){
    std::cout<<i<<"\t"<<layers_len[i]<<"\t"<<layers_dens[i]<<std::endl;
    }
  free(lenTemp); free(densTemp);
  return;  
}

/*
	Calculates and sets matter profile along neutrino direction
*/
int GLoBESCalculator::CalcDensityLayers(double coszen, bool print = false){
  
  int act_layers = 0;
  if (coszen > 1.0 | coszen < -1.0) {
    std::cout<<"Error, wrong cos(thenth) value! "<<std::endl;
    return act_layers;}
  
  for (int i =0;i<numInnerLayers;i++){
    if (coszen < PREM_cosines[i]){
      act_layers = 2*(numInnerLayers - i)-1;
      break;
    }
  }

  size_t layers = act_layers + numOuterLayers;
  std::vector<double> lengths = std::vector<double>(layers);
  std::vector<double> densities = std::vector<double>(layers);  

  for (int i = 0 ; i < int((act_layers+1)/2); i++){
    densities[act_layers/2+2 + i] = inDens[numInnerLayers - (act_layers+1)/2 +  i ];
    densities[act_layers/2+2 - i] = inDens[numInnerLayers - (act_layers+1)/2 +  i ];
    if (i == 0){
      lengths[act_layers/2+2 + i]  = 2.0*sqrt(pow(inRads[numInnerLayers - (act_layers+1)/2 +  i ],2) - pow(innerR,2) * (1.0 - pow(coszen,2)));
    } else {
      lengths[act_layers/2+2 + i]  = sqrt(pow(inRads[numInnerLayers - (act_layers+1)/2 +  i]    ,2) - pow(innerR,2) * (1.0 - pow(coszen,2)))
                             - sqrt(pow(inRads[numInnerLayers - (act_layers+1)/2 +  i - 1],2) - pow(innerR,2) * (1.0 - pow(coszen,2)));
      lengths[act_layers/2+2 - i]  = sqrt(pow(inRads[numInnerLayers - (act_layers+1)/2 +  i],    2) - pow(innerR,2) * (1.0 - pow(coszen,2)))  
                             - sqrt(pow(inRads[numInnerLayers - (act_layers+1)/2 +  i - 1],2) - pow(innerR,2) * (1.0 - pow(coszen,2)));
    }
  }
  for (int i = 0; i < numOuterLayers; i++){
    if (i ==0){
      lengths[numOuterLayers -1 - i] = -innerR*coszen +  sqrt(  pow(outRads[i],2) - pow(innerR,2)*(1.0 - pow(coszen,2)) ) - (- innerR*coszen + std::abs( innerR*coszen ));
      densities[numOuterLayers -1 -i] = outDens[i]; }
    else{
      lengths[numOuterLayers -1 -i] = -innerR*coszen +  sqrt(  pow(outRads[i],2) - pow(innerR,2)*(1.0 - pow(coszen,2)) ) 
                             -( -innerR*coszen +  sqrt(  pow(outRads[i-1],2) - pow(innerR,2)*(1.0 - pow(coszen,2)) )); }
      densities[numOuterLayers -1 -i] = outDens[i]; 
  }

  if (print) for (int i = 0 ; i<layers; i++){std::cout<<i<<"\tdL="<<lengths[i]<<"\trho="<<densities[i]<<std::endl;}
  glbSetProfileDataInExperiment(0,layers,&lengths[0],&densities[0]); 
  return layers;
}

void GLoBESCalculator::SetManualDensities(boost::python::list& length, boost::python::list& dens){
  if ( len(length) != len(dens) | len(length) ==0 | len(dens) ==0){
    std::cout<<"Densities and lengths should be of equal length and non-zero!"<<std::endl;
    std::cout<<"Sizes: "<<len(length)<<"(radius) \t"<<len(dens)<<"(dens)\n";
    return ;
  }
  int layers = len(length);
  std::vector<double> lengths = std::vector<double>(layers);
  std::vector<double> densities = std::vector<double>(layers);  
  //std::cout<<"Got layers: "<< len(dens) << std::endl;
  for (int i =0; i< layers; i++){
     lengths[i]= boost::python::extract<double>(length[i]);
     densities[i] = boost::python::extract<double>(dens[i]);
  }
  glbSetProfileDataInExperiment(0, layers, &lengths[0], &densities[0]);
  return;
}
double GLoBESCalculator::MatterProbabilityPrevBaseline(int m,int l ,int panti,double E){
  return glbProfileProbability(0,m,l,panti,E);
}


