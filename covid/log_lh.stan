    log_lh = 0.0;
    for(i in 1:D) {
      log_lh += neg_binomial_2_lpmf( incidence_cases[i] | output_incidence_cases[i], output_incidence_cases[i]/phi[1]);
      log_lh += neg_binomial_2_lpmf( incidence_deaths[i] | output_incidence_deaths[i],output_incidence_deaths[i]/phi[2]);
    }
    log_lh += multinomial_lpmf(agedistr_cases | output_agedistr_cases);
    log_lh += multinomial_lpmf(agedistr_deaths | output_agedistr_deaths);
