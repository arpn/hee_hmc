functions {

  // Linear interpolation method.
  real interpolate_linear(real x_1, real x_2, real y_1, real y_2, real x) {
    real dydx = (y_2-y_1)/(x_2-x_1);
    return dydx*(x-x_1) + y_1;
  }

  real eval_basis(int N, int i, real zs, real rho) {
    real ret = 0.0;
    if (i == 0) {
      return 1.0;
    }
    real z = zs*(1-rho)*(1+rho);

    if (i > 0) {
      // Polynomial basis
      ret = z^i - z^(N+1);
    }
    return ret;
  }

  real eval_basis_d1(int N, int i) {
    real ret = 0.0;
    if (i == 0) {
      return ret;
    }

    // Polynomial basis
    if (i == 1) {
      ret = 1.0;
    }

    return ret;
  }

  real eval_basis_d2(int N, int i) {
    real ret = 0.0;
    if (i == 0) {
      return ret;
    }

    // Polynomial basis
    if (i == 2) {
      ret = 2.0;
    }

    return ret;
  }

  real eval_basis_d3(int N, int i) {
    real ret = 0.0;
    if (i == 0) {
      return ret;
    }

    // Polynomial basis
    if (i == 3) {
      ret = 6.0;
    }

    return ret;
  }

  real eval_a(vector a, real zs, real rho) {
    int N = num_elements(a);
    real z = zs*(1-rho)*(1+rho);
    real ret = 1.0;

    for (i in 1:N) {
      // Polynomial basis
      ret += a[i]*(z^i-z^(N+1));
    }

    return ret;
  }

  real l_coef(real zs, real rho) {
    return 4*(1-rho)^3*(1+rho)^3*zs*inv_sqrt(1-zs^4*(1-rho)^4*(1+rho)^4)*inv_sqrt(6-15*rho^2+20*rho^4-15*rho^6+6*rho^8-rho^10);
  }

  real l_basis(int N, int i, real zs) {
    int N_grid = 50;
    if (zs > 0.99) {
      N_grid = 200;
    }
    real dx = inv(2*N_grid);

    real y_1 = l_coef(zs, 0.0)*eval_basis(N, i, zs, 0.0);
    real y_N = 0.0;
    real integral = y_1+y_N;

    for (j in 1:(N_grid-1)) {
      integral += 4*l_coef(zs, dx*(2*j-1))*eval_basis(N, i, zs, dx*(2*j-1));
      integral += 2*l_coef(zs, dx*(2*j))*eval_basis(N, i, zs, dx*(2*j));
    }
    integral += 4*l_coef(zs, dx*(2*N_grid-1))*eval_basis(N, i, zs, dx*(2*N_grid-1));
    integral *= dx/3;
    return integral;
  }

  real l_integral(vector a, real zs) {
    int N_grid = 50;
    if (zs > 0.99) {
      N_grid = 200;
    }
    real dx = inv(2*N_grid);

    real y_1 = l_coef(zs, 0.0)*eval_a(a, zs, 0.0);
    real y_N = 0.0;
    real integral = y_1+y_N;

    for (i in 1:(N_grid-1)) {
      integral += 4*l_coef(zs, dx*(2*i-1))*eval_a(a, zs, dx*(2*i-1));
      integral += 2*l_coef(zs, dx*(2*i))*eval_a(a, zs, dx*(2*i));
    }
    integral += 4*l_coef(zs, dx*(2*N_grid-1))*eval_a(a, zs, dx*(2*N_grid-1));
    integral *= dx/3;
    return integral;
  }

  real find_zs(vector a, real l_1, real l_2, real zs_1, real zs_2, real l) {
    real zs_low = zs_1;
    real zs_high = zs_2;
    real zs_mid;
    real l_low = l_1;
    real l_high = l_2;

    while (zs_high-zs_low > 1e-8) {
      zs_mid = (zs_low+zs_high)/2;
      real l_mid = l_integral(a, zs_mid);
      if ((l_mid < l && l < l_high) || (l_mid > l && l > l_high)) {
        zs_low = zs_mid;
        l_low = l_mid;
      } else {
        zs_high = zs_mid;
        l_high = l_mid;
      }
    }

    return interpolate_linear(l_low, l_high, zs_low, zs_high, l);
  }

  real S_basis(int N, int i, real zs, real rho_low, real rho_high, real S_low, real S_high) {
    real rho = (rho_low+rho_high)/2;
    real S_mid = inv_sqrt(1-zs^4*(1-rho)^4*(1+rho)^4)*inv_sqrt(6-15*rho^2+20*rho^4-15*rho^6+6*rho^8-rho^10)*eval_basis(N, i, zs, rho);
    if (i == 0) {
      S_mid -= rho;
    } else if (i == 1) {
      S_mid -= rho*zs*(1-rho)*(1+rho)*eval_basis_d1(N, i);
    } else if (i == 2) {
      S_mid -= rho*0.5*zs^2*(1-rho)^2*(1+rho)^2*eval_basis_d2(N, i);
    }
    S_mid *= 4*(1-rho)^(-3)*(1+rho)^(-3)*inv_square(zs);

    if ((rho_high-rho_low > 0.01) /*&& (fabs((S_low+S_high)/2-S_mid) > 1e-3)*/) {
      return S_basis(N, i, zs, rho_low, rho, S_low, S_mid) + S_basis(N, i, zs, rho, rho_high, S_mid, S_high);
    } else {
      return (rho_high-rho_low)/6*(S_low+4*S_mid+S_high);
    }
  }

  real S_integral(vector a, real zs) {
    int N = num_elements(a);
    real y_0 = 4*inv_sqrt(6)*inv_square(zs)*inv_sqrt((1-zs)*(1+zs)*(1+zs^2))*eval_basis(N, 0, zs, 0.0);
    real y_1 = 0.0;
    real ret = S_basis(N, 0, zs, 0.0, 1.0, y_0, y_1);

    for (i in 1:N) {
      y_0 = 4*inv_sqrt(6)*inv_square(zs)*inv_sqrt((1-zs)*(1+zs)*(1+zs^2))*eval_basis(N, i, zs, 0.0);
      y_1 = 2./3*zs*eval_basis_d3(N, i);
      ret += a[i]*S_basis(N, i, zs, 0.0, 1.0, y_0, y_1);
    }

    ret += -inv_square(zs);
    for (i in 1:N) {
      ret += -a[i]*2*eval_basis_d1(N, i)*inv(zs);
      ret += a[i]*eval_basis_d2(N, i)*log(zs);
    }

    return ret;
  }

}

data {
  int N_data;
  int N_basis;
  vector[N_data] l;
  vector[N_data] dS_dl;
  vector[N_data] sigma;
  vector[N_basis] a_prior_mean;
  real scale_prior_mean;
  int use_prior;
}

transformed data {
  matrix[N_basis, N_basis] R_inv;
  int N_grid = 100;
  vector[N_grid] l_0;
  matrix[N_grid, N_basis] l_mat;
  vector[N_grid] zs_grid;
  real zs_min = 0.0001;
  real zs_max = 0.99999;
  {
    for (i in 1:N_grid) {
      real zs = zs_min+(i-1)*(zs_max-zs_min)/(N_grid-1);
      zs_grid[i] = zs;
      l_0[i] = l_basis(N_basis, 0, zs);
      for (j in 1:N_basis) {
        l_mat[i, j] = l_basis(N_basis, j, zs);
      }
    }

    R_inv = inverse(qr_thin_R(l_mat)*inv_sqrt(N_data-1));
  }
}

parameters {
  // Metric coefficients
  vector[N_basis] b;
  // R^3/(4*G_N) =: scale
  real<lower=0> scale;
}

transformed parameters {
  vector[N_data] zs;
  vector[N_basis] a = R_inv*b;

  {
    real zs_temp;
    real S_temp;
    real a_temp;
    vector[N_data] S_min = rep_vector(positive_infinity(), N_data);
    real zs_1 = zs_min;
    real zs_2;
    vector[N_grid] l_grid = l_0+l_mat*a;
    real l_1 = l_grid[1];
    real l_2;

    // Check that a(z) is positive.
    for (i in 1:N_grid) {
      a_temp = eval_a(a, zs_grid[i], 0.0);
      if (a_temp < 0) {
        reject("a(z) is negative for a = ", a');
      }
    }

    for (i in 1:(N_grid-1)) {
      zs_2 = zs_grid[i+1];
      l_2 = l_grid[i+1];

      for (k in 1:N_data) {
        if ((l[k]-l_1)*(l[k]-l_2) < 0.0) {
          zs_temp = find_zs(a, l_1, l_2, zs_1, zs_2, l[k]);
          S_temp = S_integral(a, zs_temp);

          // Check if this is the minimum area phase.
          if (S_temp < S_min[k]) {
            S_min[k] = S_temp;
            zs[k] = zs_temp;
          }
        }
      }
      l_1 = l_2;
    }
  }
}

model {
  for (i in 1:N_data) {
    dS_dl[i] ~ normal(scale/zs[i]^3, sigma[i]);
  }

  if (use_prior) {
    a ~ normal(a_prior_mean, 5);
    scale ~ normal(scale_prior_mean, 5);
  }
}
