
## Compute discrete distribution of time from onset to death ----

# Fom Linton et al
linton_mean = 20.2
linton_sd = 11.6
# Get lognormal parameters from mean and sd
get_par_lnorm = function(m,s) {
  mu = log(m) - 1/2*log((s/m)^2+1)
  sigma = sqrt(log((s/m)^2+1))
  return(list(mu=mu,sigma=sigma))
}
linton_pars = get_par_lnorm(linton_mean,linton_sd)
# Discretize
gamma = numeric(60)
for(i in 1:60) {
  gamma[i] = plnorm(i+.5,linton_pars$mu,linton_pars$sigma)-plnorm(i-.5,linton_pars$mu,linton_pars$sigma)
}
# Normalize
gamma = gamma/sum(gamma)
# plot(density(rlnorm(1000,linton_pars$mu,linton_pars$sigma)))
# lines(gamma,col="red")


### ADDED BY JT ###
tmax = day_max
###################

## Format for stan
t0 = 0
S = as.numeric(day_max-day_start)
t_data = as.numeric(day_data-day_start)
ts = t_data:S
D = length(ts)
tswitch = as.numeric(day_quarantine-day_start)
data_list_model13 = list(
  K=9,
  age_dist=age_dist,
  pop_t = pop_t,
  t0=0,
  t_data = t_data,
  tswitch = tswitch,
  S=S,
  ts = ts,
  D = D,
  incidence_cases = incidence_cases,
  incidence_deaths = incidence_deaths,
  agedistr_cases = cases_tmax,
  agedistr_deaths = mort_tmax,
  p_beta=1,
  p_eta=c(1,1),
  p_pi=c(1,999),
  p_incubation=5.95, # Bi et al
  p_infectious=2.4, # Bi et al
  p_epsilon=c(1,1),
  p_phi = 1/100,
  p_rho=c(1,1),
  p_xi=1,
  p_nu=1/5,
  p_chi=5,
  G=60,
  p_gamma=gamma,
  inference=0,
  doprint=0,
  contact=c(0.810810810810811, 0.0908559523547268, 0.372736406439194, 1.27360250772, 0.200569529052988, 0.375083342749019, 0.60252680195839, 0.0934189610338407, 0.0225225225225225, 0.0904095466183592, 2.4392523364486, 0.140093983348316, 0.706545801082683, 0.942937990573664, 0.27920963239528, 0.326366336169345, 0.196893495540358, 0.106045179398683, 0.289504965045504, 0.109348487445688, 1.76086956521739, 0.923069180041088, 0.93772012267962, 0.965186137047983, 0.274120168579709, 0.116564256844925, 0.0773400190233669, 0.91820215964926, 0.511898453884688, 0.85680985412458, 2.70542635658915, 1.41323192857305, 0.993399938008648, 0.719603621821669, 0.146103509716984, 0.07633130138862, 0.13119227828341, 0.619819944222649, 0.789700390093264, 1.28218991206025, 2.17699115044248, 1.1199461877854, 0.514253349451317, 0.496466649026704, 0.101504389707241, 0.259078294801222, 0.193808465356441, 0.858341528544101, 0.951750199084178, 1.18265232149625, 2.31730769230769, 0.977037933291252, 0.606164987575222, 0.4393566902894, 0.552747314447092, 0.300880970126328, 0.323770338070664, 0.915670466885606, 0.721247101248993, 1.29765260904839, 2.76146788990826, 0.959867553314515, 0.340125585278128, 0.121809161946423, 0.25799743320884, 0.1956843612527, 0.264241585561661, 0.989672909331423, 1.14428055461657, 1.36428769674242, 1.96363636363636, 1.0266513139522, 0.0447824174066075, 0.211894911958445, 0.197988778289041, 0.210517772531686, 0.308554588199316, 1.26474943927563, 0.737190168823191, 1.56555579008225, 2.0625)
)

# Create data and bash files ----
# symptomatic = 80% (A) (Bi et al) 
data_list_model13$p_psi=c(71,18) 

# symptomatic = 50% (B) (Diamond princess https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.10.2000180) 
#data_list_model13$p_psi=c(522,114)
