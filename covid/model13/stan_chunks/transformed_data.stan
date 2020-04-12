
  real x_r[3+K*K+K]; // 4 parameters + K*K contact matrix parameters + K age_dist parameters
  int x_i[1] = {K};
  real init[K*4] = rep_array(0.0, K*4); // initial values
  x_r[1] = tswitch;
  x_r[2] = p_incubation;
  x_r[3] = p_infectious;
  x_r[4:(3+K*K)] = contact;
  for(k in 1:K) {
    x_r[3 + K*K + k] = age_dist[k];
  }
