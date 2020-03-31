functions {
  real[] 	cbf(real t,
              real[] x,
              real[] theta,
              real[] y_r,
              int[] y_i) {
    real dxdt[8];
    
    dxdt[1] = - theta[1]*theta[12]*x[1]*x[6]/(theta[17]+x[1]) - theta[2]*theta[14]*x[1]*x[7]/(theta[19]+x[1]);
    dxdt[2] = - theta[3]*theta[13]*x[2]*x[6]/(theta[18]+x[2]);
    dxdt[3] = theta[4]*theta[12]*x[1]*x[6]/(theta[17]+x[1]) + theta[5]*theta[13]*x[2]*x[6]/(theta[18]+x[2]) - theta[6]*theta[15]*x[3]*x[8]/(theta[20]+x[3]);
    dxdt[4] = theta[7]*theta[14]*x[1]*x[7]/(theta[19]+x[1]) - theta[8]*theta[16]*x[4]*x[8]/(theta[21]*x[8]+x[4]);
    dxdt[5] = theta[9]*theta[14]*x[1]*x[7]/(theta[19]+x[1]) + theta[10]*theta[15]*x[3]*x[8]/(theta[20]+x[3]) +theta[11]*theta[16]*x[4]*x[8]/(theta[21]*x[8]+x[4]);
    dxdt[6] = theta[12]*x[1]*x[6]/(theta[17]+x[1]) + theta[13]*x[2]*x[6]/(theta[18]+x[2]) - theta[22]*x[6]*x[3];
    dxdt[7] = theta[14]*x[1]*x[7]/(theta[19]+x[1]) - theta[23]*x[7]*x[4];
    dxdt[8] = theta[15]*x[3]*x[8]/(theta[20]+x[3]) + theta[16]*x[4]*x[8]/(theta[21]*x[8]+x[4]) - theta[24]*x[8]*x[5]^2;
    
    return dxdt;
  }
}
data {
  int<lower=1> T;
  real<lower=0> x[T,8];
  real t0;
  real ts[T];
  real x0[8];
}
transformed data {
  real y_r[0];
  int y_i[0];
}
parameters {
  vector<lower=0>[8] sigma;
  real<lower=0.005999985,upper=1.998029> mu1;
  real<lower=0.00399848,upper=0.404> mu2;
  real<lower=0.4559996,upper=1.998033> mu3;
  real<lower=0.7739987,upper=1.999988> mu4;
  real<lower=0,upper=0.01000002> mu5;
  real<lower=44.7998,upper=799.9959> ks1;
  real<lower=3.678926E-7,upper=0.8000041> ks2;
  real<lower=53.59992,upper=380.8> ks3;
  real<lower=38.39991,upper=132.8> ks4;
  real<lower=0.3999882,upper=799.6> ks5;
  real<lower=0.03391768,upper=0.05500034> k1;
  real<lower=0.005047708,upper=0.006317451> k2;
  real<lower=0.006030796,upper=0.007347895> k3;
  real<lower=13.99994,upper=1999.843> yc1;
  real<lower=19.99999,upper=34.00002> yc2;
  real<lower=13.99998,upper=1998.618> yc3;
  real<lower=1.997585,upper=932.0009> yc4;
  real<lower=1.999979,upper=552.0002> yc5;
  real<lower=675.9997,upper=1150> yc6;
  real<lower=9.70224,upper=10.00001> yc7;
  real<lower=0.9998504,upper=1999> yc8;
  real<lower=5.926618,upper=6.00001> yc9;
  real<lower=1.998582,upper=70.0001> yc10;
  real<lower=0.9999872,upper=1999> yc11;
}
transformed parameters {
  real x_hat[T,8];
  {
    real theta[24];
    theta[1] = yc1;
    theta[2] = yc2;
    theta[3] = yc3;
    theta[4] = yc4;
    theta[5] = yc5;
    theta[6] = yc6;
    theta[7] = yc7;
    theta[8] = yc8;
    theta[9] = yc9;
    theta[10] = yc10;
    theta[11] = yc11;
    theta[12] = mu1;
    theta[13] = mu2;
    theta[14] = mu3;
    theta[15] = mu4;
    theta[16] = mu5;
    theta[17] = ks1;
    theta[18] = ks2;
    theta[19] = ks3;
    theta[20] = ks4;
    theta[21] = ks5;
    theta[22] = k1;
    theta[23] = k2;
    theta[24] = k3;
    
    x_hat = integrate_ode_bdf(cbf, x0, t0, ts, theta, y_r, y_i,1.0E-6, 1.0E-6, 1.0E8);
  }
}
model{
  
  sigma~cauchy(0,2.5);
  for (t in 1:T)
    x[t]~normal(x_hat[t],sigma);
}