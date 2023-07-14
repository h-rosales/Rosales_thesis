
functions {
  real kf_lpdf(matrix mu, real tau2, vector alpha, vector beta, real sigma2, int time_all, int time_obs, int zrow, int murow, int J, int[] p, matrix Z,
  matrix a, matrix P, matrix[] w){
    int k = 0;    // increment Z matrix by animals observed or 1 when NULL
    int i = 1;    // increment by 1 for each time observed
    int l = 0;    // increment mu by animals observed
    //int all = 1; // increment for each time point
    matrix[J, 2] a_t = a;   // to overwrite a in KF algorithm 
    matrix[J, J] P_t[2];    // to overwrite P in KF algorithm
    real LL;      // Log likelihood
    matrix[time_obs, 2] det_F; // constant inside log likelihood for easier calculation
    matrix[J, J] Q;          // part of kf
    matrix[J, J] L[2]; // part of KF
    matrix[J, J] L_tran[2]; // increase effiency
    matrix[J, J] w_bar;
    matrix[J, J] omega[time_all]; //precision matrix
    matrix[J, J] omega_inv;  // inverse of omega
    matrix[J, J] A; // Attraction matrix
     matrix[J, J] A_tran; //transpose of attraction matrix
    
    
    P_t[1] = P;
    P_t[2] = P;
    
    for(all in 1:time_all){
      //matrix[J, J] w_mat = w[all, 1:J, 1:J];       // used to create Attraction matrix
     omega[all, :, :] = -1 * alpha[all] * w[all, :, :];      // Precision Matrix
      
      // creating omega and precision matrix
      for(m in 1:J){
        //int s = m + 1;
        omega[all, m, m] = sum(w[all, m, :]);
        w_bar[m, :] = w[all, m, :] / omega[all, m, m];
      }
      A = add_diag((beta[all] * w_bar), (1-beta[all]));
      omega_inv = inverse(omega[all, :, :]);
      Q = sigma2 * omega_inv;        // Covariance Matrix
      
      if(sum(Z[(k+1), 1:J]) > 0){
        matrix[p[i], 2] v;      // create v 1 col is latitude 2 col is longitude
        matrix[2, p[i]] v_tran; // increase efficiency
        matrix[p[i], p[i]] F[2];  // create F into array
        matrix[p[i], p[i]] F_inv[2];  // create F inverse into array
        matrix[J, p[i]] K[2]; // part of KF
        matrix[J, p[i]] Z_tran = Z[(k+1):(k+p[i]), 1:J]'; // increase efficiency
        
        //Kalman Filter

        v[:, :] = mu[(l+1):(l+p[i]), :] - (Z[(k+1):(k+p[i]), :] * a_t[:, 1:2]);
        //v[1:p[i], 2] = mu[(l+1):(l+p[i]), 2] - (Z[(k+1):(k+p[i]), 1:J] * a_t[1:J, 2]);
        v_tran[1, :] = v[:, 1]';
        v_tran[2, :] = v[:, 2]';

        
        F[1] = add_diag(quad_form(P_t[1], Z_tran), tau2);
        F[2] = add_diag(quad_form(P_t[2], Z_tran), tau2);
        F_inv[1] = inverse(F[1]);
        F_inv[2] = inverse(F[2]);
        det_F[i,1] = log_determinant(F[1]) + quad_form(F_inv[1], v[:,1]);
        det_F[i,2] = log_determinant(F[2]) + quad_form(F_inv[2], v[:,2]);
        
        K[1] = A * P_t[1] * Z_tran / F[1];
        K[2] = A * P_t[2] * Z_tran / F[2];

        L[1] = A - (K[1] * Z[(k+1):(k+p[i]), :]);
        L[2] = A - (K[2] * Z[(k+1):(k+p[i]), :]);
        L_tran[1] = L[1]';
        L_tran[2] = L[2]';

        a_t[:, 1] = (A * a_t[:, 1]) + (K[1] * v[:, 1]);
        a_t[:, 2] = (A * a_t[:, 2]) + (K[2] * v[:, 2]);

        P_t[1] = (A * P_t[1] * L_tran[1]) + Q;
        P_t[2] = (A * P_t[2] * L_tran[2]) + Q;
        

        // increments
        k += p[i];
        l += p[i];
        i += 1;
        //all += 1;
      }
      else{
        A_tran = A'; // Attraction matrix transpose;
        
        a_t[:, 1:2] = A * a_t[:, 1:2];
        //a_t[1:J, 2] = A * a_t[1:J, 2];
        P_t[1] = quad_form(P_t[1], A_tran) + Q;
        P_t[2] = quad_form(P_t[2], A_tran) + Q;
        //P_t[1] = A * P_t[1] * A_tran + Q;
        //P_t[2] = A * P_t[2] * A_tran + Q;
        k += 1;
      }
    }
    LL = -0.5 * sum(det_F);
    return LL;
  }
}

data {
  int<lower=0> time_all;          // length of all time
  int<lower=0> time_obs;          // length of time observed
  int<lower=0> zrow;              // number of rows of Z matrix
  int<lower=0> murow;             // total amount of observations x animals id
  int<lower=0> J;                 // amount of individuals
  int p[time_obs];                // animals observed at time t
  // real tau2;                      // measurment error to create matrix H
  // real sigma2;                   // measurement error to create Q
  // matrix[J, J] Q;                 // Precision matrix 
  matrix[zrow, J] Z;              // Z transformation matrix
  matrix[J, 2] a;                 // a initial values for mean location
  matrix[J, J] P;                 // initial covariance matrix 
  // real alpha;                     // alignment
  // real beta;                      // attraction
  // matrix[J, J] A;                 // attraction transformation matrix
  matrix[murow, 2] mu;            // longitude and latitude observations
  matrix[J,J] w[time_all];
  matrix[time_all, 3] x;                    // covariates
  vector[3] delta_alpha_mean;
  vector[3] delta_beta_mean;
  cov_matrix[3] delta_alpha_cov;
  cov_matrix[3] delta_beta_cov;
}

parameters{
  real<lower=0> tau2;            // measurement error to create H
  vector[3] delta_alpha;           // alignment
  vector[3] delta_beta;            // attraction
  real<lower=0> sigma2;          // scale parameter to create Q
}
transformed parameters{
  vector[time_all] alpha = inv_logit(x * delta_alpha);   //delta alpha for logit
  vector[time_all] beta = inv_logit(x * delta_beta);     //delta beta for logit
}
model {
   tau2 ~ inv_gamma(3, 0.5);
   sigma2 ~ inv_gamma(3, 8); 
   delta_alpha ~ multi_normal(delta_alpha_mean, delta_alpha_cov);
   delta_beta ~ multi_normal(delta_beta_mean, delta_beta_cov);
  target += kf_lpdf(mu | tau2, alpha, beta, sigma2, time_all, time_obs, zrow, murow, J, p, Z, a, P, w[]);
}



 

