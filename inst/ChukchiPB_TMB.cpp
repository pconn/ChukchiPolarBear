// Polar bear TMB model
#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace Eigen;
  using namespace density;
  
  // options vec
  DATA_FACTOR( Options_vec );
  // Slot 0: compute SE on lambda? 
  // Slot 1: use latent tracks variable in DM for counts?
  
  // Data
  DATA_VECTOR( N_i );      //number of photographs examined at each sampled grid cell
  DATA_VECTOR( T_i );       	// number of photographs with polar bear tracks
  DATA_VECTOR( C_i );      // counts of 'on effort' polar bears in each cell, time visited
  DATA_VECTOR( U_i );      // counts of 'off effort' polar bears in each cell, time visited
  DATA_VECTOR( P_i );      // Proportion of survey unit that is surveyed for each observation i
  DATA_VECTOR( A_s );      // Relative area of each survey unit s (can set = 1.0 if all the same size)
  DATA_VECTOR( G_min1 );      // group size minus 1 vector
  DATA_IVECTOR( S_i ); // Site-time index for each sample
  DATA_IMATRIX( Y_s ); //indicator for which sites sampled/not sampled
  DATA_MATRIX( Xt_s );  //design matrix for fixed effects in tracks model
  DATA_MATRIX( Xc_s );  //design matrix for fixed effects in bear model
  // Russian data
  DATA_VECTOR( T_i_ru );       	// number of photos with polar bear tracks per grid cell
  DATA_VECTOR( C_i_ru );      // counts of 'on effort' polar bears in each cell, time visited
  DATA_VECTOR( P_i_ru );      // Proportion of survey unit that is surveyed for each observation i
  DATA_VECTOR(KM_i_ru );      // km surveyed in grid cell in a certain day (russia only)
  DATA_VECTOR(Distances);   // Distance data for Russian surveys
  DATA_SCALAR(g0);        //assumed g(0) for Russian surveys (proportion seen under aircraft)
  DATA_IVECTOR( S_i_ru ); // Site-time index for each sample
  DATA_SCALAR(strip_width); //width of Russian strip width for distance sampling calculations

  DATA_SCALAR(a_photo_us);  // average area of a US photo (km^2) - used in track proportion calculations
 
  //DATA_SCALAR(p);              //detection probability of polar bears  - fixed for bootstrap
  DATA_INTEGER(p_trials);
  DATA_INTEGER(p_detects);
  
  DATA_SPARSE_MATRIX(S);//Penalization block diagonal matrix 
  DATA_IVECTOR(Sdims);   //Dimensions of each smooth
  DATA_SPARSE_MATRIX(designMatrixForReport);//Design matrix for report of splines
  
  DATA_VECTOR(I_water);
  DATA_VECTOR(I_us);
  DATA_VECTOR(I_PBSG);
  DATA_VECTOR(I_regehr);
  
    

  // Parameters 
  PARAMETER_VECTOR(Beta);              // fixed effects on tracks
  PARAMETER_VECTOR(Alpha);     //fixed effects on bears
  PARAMETER(logN);           //log (polar bear abundance)
  PARAMETER(log_adjust);     // log (off effort sighting adjustment)
  PARAMETER(log_group_size_min1);   // log (group size minus 1) - for zero truncated Poisson
  PARAMETER(log_w);  // transect width for track observations in Russia
  PARAMETER(log_sigma); // log(half-normal distance sampling sd)
  PARAMETER_VECTOR(log_lambda);//Penalization parameters for splines
  PARAMETER(tracks_eff);  //linear effect of tracks (need to set 0 w/ Map to turn off)
  PARAMETER(logit_p);  //detection probablity for U.S. thermal sensors

  // derived sizes
  int n_i = N_i.size();  //number of spatio-temporally distinct samples
  int n_s = Y_s.col(0).size();
  int n_t = Y_s.row(0).size();  
  //
 // int n_i_ru = KM_i_ru.size();  //number of spatio-temporally distinct samples

  int n_i_ru = KM_i_ru.size();  //number of spatio-temporally distinct samples
  int n_st = Xt_s.col(0).size();
  int n_bt = Xt_s.row(0).size();
  int n_bc = Xc_s.row(0).size();
  int n_d = Distances.size();
  
  int n_lambda = log_lambda.size();
  
  // global stuff
   vector<Type> jnll_comp(9);
   jnll_comp.setZero();

  // Predicted densities
  vector<Type> Z_s(n_st);  //density of polar bear tracks
  vector<Type> linpredZ_s(n_st);  //density of polar bear tracks - log scale
  
  vector<Type> Ptrack_i(n_i);  //probability of a track in a US photo
  vector<Type> Lambda_ru_i(n_i_ru);   // track intensity for Rus data
  
  vector<Type> lambda=exp(log_lambda);
  linpredZ_s = Xt_s * Beta;
 
  for(int ist=0; ist<n_st; ist++){
    Z_s(ist) = exp( linpredZ_s(ist) );
  }

  for(int i=0;i<n_i;i++){
    Ptrack_i(i) = 1.0 - exp(-Z_s(S_i(i))*a_photo_us);
  }
  
  Type p = 1/(1+exp(-logit_p));

  vector<Type> Dev_tracks(n_i);
  vector<Type> Dev_tracks_ru(n_i_ru);
  Dev_tracks_ru.setZero();
  // Probability of track counts
  for(int i=0; i<n_i; i++){
    if(N_i(i)>0){
      jnll_comp(0) -= dbinom( T_i(i),N_i(i), Ptrack_i(i), true );
      if(T_i(i)==0) Dev_tracks(i) = (N_i(i)-T_i(i))*log((N_i(i)-T_i(i))/(N_i(i)-(N_i(i)*Ptrack_i(i))));
      if(N_i(i)==T_i(i)) Dev_tracks(i) = T_i(i)*log(T_i(i)/(N_i(i)*Ptrack_i(i)));
      if(T_i(i)>0 & N_i(i)!=T_i(i)) Dev_tracks(i) = T_i(i)*log(T_i(i)/(N_i(i)*Ptrack_i(i)))+(N_i(i)-T_i(i))*log((N_i(i)-T_i(i))/(N_i(i)-(N_i(i)*Ptrack_i(i))));
      
    }
    else Dev_tracks(i) = 0.0;
  }

  Type w = exp(log_w);   // effective transect width (km) for track detection in Russian survey
  for(int i=0;i<n_i_ru;i++){
    Lambda_ru_i(i) = Z_s(S_i_ru(i))*KM_i_ru(i)*w;
    jnll_comp(1) -= dpois( T_i_ru(i), Lambda_ru_i(i), true );
    if(T_i_ru(i)>0)Dev_tracks_ru(i) = T_i_ru(i) * log(T_i_ru(i)/Lambda_ru_i(i)) - (T_i_ru(i)-Lambda_ru_i(i));
    else Dev_tracks_ru(i) = Lambda_ru_i(i);  //what things work out to when T_i_ru(i)=0
  } 

  // group size likelihood
  Type group_size_min1 = exp(log_group_size_min1);
  Type n_g = G_min1.size();
  for(int i=0;i<n_g;i++){
    jnll_comp(2) -= dpois(G_min1(i),group_size_min1,true);
  }

  //distance data likelihood
  Type sigma = exp(log_sigma);
  vector<Type> P_dist(n_d);
  Type denom = 2.0*(pnorm(strip_width,Type(0.0),sigma)-0.5);
  for(int i=0;i<n_d;i++){
    P_dist(i) = dnorm(Distances(i),Type(0.0),sigma) / denom;
    jnll_comp(3) -= log(P_dist(i));  //half normal detection function
  }
  Type pp = 2.0*g0*(pnorm(strip_width,Type(0.0),sigma)-0.5)/(2.0*strip_width*dnorm(Type(0.0),Type(0.0),sigma));  
  
  //detection prior for U.S. thermal sensors
  Type p_det = p_detects;
  Type p_tr = p_trials; 
  jnll_comp(3) -= dbinom(p_det,p_tr,p,true);

  // polar bear abundance
  Type N = exp(logN);
  
  //spline prior
  int k=0;  // Counter
  for(int i=0;i<n_lambda;i++){
    int m_i = Sdims(i);
    vector<Type> beta_i = Alpha.segment(k,m_i);       // Recover betai
    SparseMatrix<Type> S_i = S.block(k,k,m_i,m_i);  // Recover Si
    jnll_comp(4) -= Type(0.5)*m_i*log_lambda(i) - 0.5*lambda(i)*GMRF(S_i).Quadform(beta_i); //note from Devin: m_i would need to be rank(S) if S_i not full rank (e.g. thin plate)
    k += m_i;
    jnll_comp(4) -= (Type(-0.95)*log_lambda(i)-0.005*exp(log_lambda(i)));  //gamma(0.05,0.005) prior like in Jagam
  }
  
  
  // assemble multinomial cell probabilities...
  vector<Type> Pi_s(n_st);
  vector<Type> Lambda_s(n_st);
  vector<Type> linpredC_s(n_st);

  linpredC_s = Xc_s * Alpha + I_water*Type(-20.0) + tracks_eff * Z_s;
  for(int ist=0;ist<n_st;ist++){
    Pi_s(ist)=exp(linpredC_s(ist));
  }
  Type cur_sum;
  for(int it=0;it<n_t;it++){
    cur_sum = 0;
    for(int is=0;is<n_s;is++){
      cur_sum = cur_sum + Pi_s(it*n_s+is)*A_s(is);
    }
    for(int is=0;is<n_s;is++){
      Pi_s(it*n_s+is)=Pi_s(it*n_s+is)*A_s(is)/cur_sum;
      Lambda_s(it*n_s+is)=N*Pi_s(it*n_s+is);
    }
  }
  
  // Probability of bear counts in US
  vector<Type> E_count_on(n_i);
  vector<Type> E_count_off(n_i);
  E_count_on = E_count_on.setZero();
  E_count_off= E_count_off.setZero();
  Type adjust = exp(log_adjust);
  vector<Type> Dev_on(n_i);
  vector<Type> Dev_off(n_i);
  for(int i=0; i<n_i; i++){
    E_count_on(i)=Lambda_s(S_i(i))*p*P_i(i);
    E_count_off(i)=E_count_on(i)*adjust;
    jnll_comp(5) -= dpois( C_i(i), E_count_on(i), true );
    jnll_comp(6) -= dpois( U_i(i), E_count_off(i), true );
    if(C_i(i)==0)Dev_on(i) = -(C_i(i)-E_count_on(i));
    else Dev_on(i) = C_i(i)*log(C_i(i)/E_count_on(i))-(C_i(i)-E_count_on(i));
    if(U_i(i)==0)Dev_off(i) = -(U_i(i)-E_count_off(i));
    else Dev_off(i) = U_i(i)*log(U_i(i)/E_count_off(i))-(U_i(i)-E_count_off(i));
  }
  
  // Probability of bear counts in Russia
  vector<Type> E_count_on_ru(n_i_ru);
  vector<Type> E_count_off_ru(n_i_ru);
  E_count_on_ru = E_count_on_ru.setZero();
  E_count_off_ru= E_count_off_ru.setZero();
  vector<Type> Dev_on_ru(n_i_ru);
  vector<Type> Dev_off_ru(n_i_ru);
  for(int i=0; i<n_i_ru; i++){
    E_count_on_ru(i)=Lambda_s(S_i_ru(i))*pp*P_i_ru(i);
    if(E_count_on_ru(i)>0)jnll_comp(7) -= dpois( C_i_ru(i), E_count_on_ru(i), true );
    if(C_i_ru(i)==0)Dev_on_ru(i) = -(C_i_ru(i)-E_count_on_ru(i));
    else Dev_on_ru(i) = C_i_ru(i)*log(C_i_ru(i)/E_count_on_ru(i))-(C_i_ru(i)-E_count_on_ru(i));
  }
  // Total objective
  Type group_size = group_size_min1+1;
  Type jnll = jnll_comp.sum();
  
  Type dev_off = 2*Dev_off.sum();
  Type dev_on = 2*Dev_on.sum();
  Type dev_off_ru = 2*Dev_off_ru.sum();
  Type dev_on_ru = 2*Dev_on_ru.sum();
  Type dev_tracks_us = 2*Dev_tracks.sum();
  Type dev_tracks_rus = 2*Dev_tracks_ru.sum();
  Type N_all = N * group_size;
  
  vector<Type> N_PBSG_t(n_t);
  vector<Type> N_Regehr_t(n_t);
  vector<Type> N_US_t(n_t);
  vector<Type> N_rus_t(n_t);
  N_PBSG_t.setZero();
  N_Regehr_t.setZero();
  N_US_t.setZero();
  N_rus_t.setZero();

  int counter=0;
  for(int it=0; it<n_t;it++){
    //std::cout<<"it "<<it<<std::endl;
    for(int is = 0; is<n_s; is++){
      N_PBSG_t(it)+= group_size*Lambda_s(it*n_s+is)*I_PBSG(is);
      N_US_t(it)+= group_size*Lambda_s(it*n_s+is)*I_us(is);
      N_rus_t(it)+= group_size*Lambda_s(it*n_s+is)*(1-I_us(is));
      N_Regehr_t(it)+= group_size*Lambda_s(it*n_s+is)*I_regehr(is);
    }
  }
  Type N_Regehr = N_Regehr_t.sum()/n_t;
  Type N_US = N_US_t.sum()/n_t;
  Type N_rus = N_rus_t.sum()/n_t;
  Type N_PBSG = N_PBSG_t.sum()/n_t;

  //Reporting
   REPORT( Z_s );
   REPORT( Beta );
   REPORT( Alpha );
   REPORT( jnll_comp );
   REPORT( jnll );
   REPORT( N );
   REPORT( Lambda_s);
   REPORT( Ptrack_i);
   REPORT( Lambda_ru_i);
   REPORT( dev_off);
   REPORT( dev_on);
   REPORT( dev_off_ru);
   REPORT( dev_on_ru);
   REPORT( dev_tracks_us);
   REPORT( dev_tracks_rus);
   REPORT( w);
   REPORT( adjust);
   REPORT(N_all);  //includes cubs / group sizes >1
   REPORT(E_count_on); 
   REPORT(E_count_off);
   REPORT(E_count_on_ru);
   REPORT(E_count_off_ru);
   REPORT(p);
   REPORT(pp);
   REPORT(log_sigma);
   REPORT(sigma);
   REPORT(group_size);
   REPORT(log_lambda);
   REPORT(tracks_eff);
   REPORT(N_Regehr_t);
   REPORT(N_US_t);
   REPORT(N_rus_t);
   REPORT(N_PBSG_t);
   
   

  // Bias correction output
  ADREPORT( N );
  ADREPORT( N_all);
  ADREPORT(Beta);
  ADREPORT(Alpha);
  vector<Type> splineForReport = designMatrixForReport*Alpha;
  ADREPORT(splineForReport);
  ADREPORT(tracks_eff);
  ADREPORT(N_Regehr);
  ADREPORT(N_US);
  ADREPORT(N_rus);
  ADREPORT(N_PBSG);
  if(Options_vec(0)==1){
    ADREPORT( Lambda_s);
  }
  return jnll;
  
 
}
