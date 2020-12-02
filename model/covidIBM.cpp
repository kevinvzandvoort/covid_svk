// covidIBM.cpp

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppGSL)]]

#include <iostream>
#include <vector>
#include <cstdlib>
#include "./randomizer.h"
#include <Rcpp.h>
#include <omp.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "./ibm.h"

int main(){}

//Function for model
Rcpp::List RunSimulation(Rcpp::List parameters, int seed){
  std::cout << "Run simulation" << std::endl;

  // load parameter values from parameters list
  Rcpp::DataFrame population_data = parameters["population"];
  int n_households = parameters["n_households"];
  int initial_inf = parameters["n_init_infected"];
  int agegroups = parameters["agegroups"];
  int times = parameters["days"];
  int tstep = parameters["tstep_day"];
  bool model_quick = parameters["model_quick"];
  int iv_stop = parameters["iv_stop"];
  double beta_hh = parameters["contact_household"];
  Rcpp::NumericMatrix beta = parameters["contact_matrix"];
  Rcpp::List interventions = parameters["interventions"];
  Rcpp::List testpars = parameters["testpars"];

  Rcpp::List lockdown = parameters["lockdown"];
  int lockdown_start = lockdown["start_value"];
  double lockdown_effect = lockdown["effect"];
  int lockdown_test_after = lockdown["test_after"];

  /*
   allocate memory for household- and individual-pointer vectors
   use pointers to minimize overhead
   */
  std::vector<Household*> households;
  households.reserve(n_households);
  std::vector<Individual*> individuals;
  int n_pop = population_data.nrows();
  individuals.reserve(n_pop);

  // use temporary scope to load data
  {
    // TODO is clone necessary here?
    Rcpp::NumericVector df_part_id = population_data["participant_id"];
    Rcpp::NumericVector df_household_id = population_data["household_id"];
    Rcpp::NumericVector df_age_group = population_data["age_group"];
    Rcpp::NumericVector df_dE = population_data["dE"];
    Rcpp::NumericVector df_dP = population_data["dP"];
    Rcpp::NumericVector df_dC = population_data["dC"];
    Rcpp::NumericVector df_dS = population_data["dS"];
    Rcpp::NumericVector df_y = population_data["y"];
    Rcpp::NumericVector df_u = population_data["u"];

    int household_id = -1;
    int h = -1;
    for(int i=0; i<n_pop; i++){
      // add household pointer to household vector
      if(df_household_id[i] != household_id){
        household_id = df_household_id[i];
        h++;
        households.emplace_back(new Household(household_id, 0));
      }

      // add individual pointer to household and individual vector
      Individual* indiv = new Individual(df_part_id[i], household_id, df_age_group[i], households[h], df_dE[i], df_dP[i], df_dC[i], df_dS[i], df_y[i], df_u[i]);
      individuals.emplace_back(indiv);
      households[h]->addMember(indiv);
    }
  }

  //this is only used if model_quick == true
  agegroups -= 1;
  int Iclasses = 3;
  int infectious_strata[agegroups][Iclasses];
  for(int a=0; a<agegroups; a++){
    for(int i=0; i<Iclasses; i++){
      infectious_strata[a][i] = 0;
    }
  }

  // vector to store infected individuals (pointers)
  std::vector<Individual*> susceptibles;
  std::vector<Individual*> preinfectious;
  std::vector<Individual*> infectious;

  // pre-allocate memory so mem addresses never change (using vectors)
  susceptibles.reserve(n_pop);
  preinfectious.reserve(n_pop);
  infectious.reserve(n_pop);

  // ensure random number generator works when using multiple threads
  std::vector<Randomizer*> randomizers;
  randomizers.reserve(omp_get_max_threads()+1);

  for(int p=0; p <= omp_get_max_threads(); p++){
    Randomizer* rndm = new Randomizer(seed+p);
    randomizers.emplace_back(rndm);
  }
  
  // random number generator used in the main thread
  Randomizer* random_generator = randomizers[omp_get_max_threads()];

  // randomly inocculate individual(s)
  while(infectious.size() < initial_inf){
    int rnd_infected = random_generator->UniformInt(0, n_pop-1);

    if(individuals[rnd_infected]->m_status==0){
      if(random_generator->Binomial(1, individuals[rnd_infected]->m_y) == 1){
        individuals[rnd_infected]->m_status=Individual::infected_clinical;
      } else {
        individuals[rnd_infected]->m_status=Individual::infected_preclinical;
      }
      individuals[rnd_infected]->m_infected_at = 0.0;
      individuals[rnd_infected]->m_infectious_at = 0.0;
      individuals[rnd_infected]->m_infected_by = 0;
      infectious.emplace_back(individuals[rnd_infected]);
    }
  }
  
  // array to track status of interventions
  int iv_vals[interventions.size()][7];
  for(int i=0; i < interventions.size(); i++){
    iv_vals[i][0] = -1;
  }
  int iv = 0;

  // loop through all timesteps
  for(int t=1; t<=times*tstep; t++){
    double t_day = t/tstep;
    
    // track number in susceptible and infected compartment
    // calculated at each timestap to allow sharing of susceptible vector across threads
    susceptibles.clear();
    susceptibles.reserve(n_pop);
    preinfectious.clear();
    preinfectious.reserve(n_pop);
    infectious.clear();
    infectious.reserve(n_pop);
    
    //track total number in each infectious strata, if model_quick is used
    if(model_quick){
      for(int a=0; a<agegroups; a++){
        for(int i=0; i<Iclasses; i++){
          infectious_strata[a][i] = 0;
        }
      }
    }
    
    //check whether lockdown will be implemented today
    bool start_of_lockdown = false;
    if(t_day == lockdown_start){
      std::cout << "Starting lockdown on day " << t_day << std::endl;
      start_of_lockdown = true;
      if(interventions.size() > 0){
        iv_vals[0][0] = t_day + lockdown_test_after;
      }
    }
    
    //loop through all individuals
    // - count number infectious and susceptible
    // - progress Status of those infected
    // - progress status of those in quarantine
    for(int i = 0; i < individuals.size(); ++i){
      Individual* member = individuals[i];
      bool change_status = member->progressStatus(random_generator, t_day);
      
      //lockdown_effect is reduction in transmissibility due to lockdown
      //effectively reduces Rt
      if(start_of_lockdown){
        member->m_lockdown_effect = lockdown_effect;
      }
      
      member->checkQuarantine();
    
      if( member->m_status == Individual::susceptible){
        susceptibles.emplace_back(member);
      } else if( member->m_status == Individual::preinfectious && !change_status) {
        preinfectious.emplace_back(member);
      } else if(
          member->m_status == Individual::infected_preclinical |
            member->m_status == Individual::infected_clinical |
            member->m_status == Individual::infected_subclinical
      ){
        infectious.emplace_back(member);
        
        //only track those not quarantining for community transmission
        if(!member->m_quarantine && model_quick){
          infectious_strata[member->m_age][member->m_status - 2] += 1;
        }
      }
    }
    
    //should intervention be started today?
    if(iv < interventions.size()){

      Rcpp::List intervention = interventions[iv];
      int iv_start_type = intervention["start_type"];
      double iv_start_value = intervention["start_value"];

      if(iv_vals[iv][0] == -1){
        //iv_start_type 0 starts intervention if prevalence level is reached
        if(iv_start_type == 0){
          double prevalence = infectious.size()/(double)n_pop;
          if(prevalence >= iv_start_value){
            iv_vals[iv][0] = t_day;
          }
        } else {
          iv_vals[iv][0] = iv_vals[iv-1][0] + iv_start_value;
        }
      }

      if(t_day == iv_vals[iv][0]){
        //implement mass testing
        std::cout << "Implementing mass testing on day " << t_day << std::endl;

        int masstest_results[2] = {0, 0};
        
        // how many are already in quarantine? Need to adjust proportion of those not in quarantine who 
        //  attend testing
        int nquarantine = 0;
        for(int i = 0; i < individuals.size(); ++i){
          nquarantine += individuals[i]->m_quarantine;
        }
        double pquarantine = nquarantine/(double)n_pop;
        
        //do individuals attend testing and quarantine?
        for(int i = 0; i < individuals.size(); ++i){
          Individual* member = individuals[i];

          bool attend = member->massTestAttend(intervention, random_generator, pquarantine);
          bool positive = false;
          if(attend){
            positive = member->massTest(intervention, testpars, random_generator, attend);
          }
          masstest_results[0] += attend;
          masstest_results[1] += positive;
        }
        
        //track summary data of mass testing
        iv_vals[iv][1] = masstest_results[0];
        iv_vals[iv][2] = masstest_results[1];
        iv_vals[iv][3] = nquarantine;
        iv_vals[iv][4] = susceptibles.size();
        iv_vals[iv][5] = preinfectious.size();
        iv_vals[iv][6] = infectious.size();
        
        std::cout << "Mass testing " << iv+1 << " completed. Attendance: " << masstest_results[0] << "; Positive: " << masstest_results[1] << std::endl;
        iv += 1;

      }
    } else if(iv == interventions.size() & t_day == (iv_vals[iv-1][0] + iv_stop)){
      std::cout << "Finished simulation early after interventions have completed" << std::endl;
      
      //TODO: write single export function so code is not duplicated
      std::cout << "Data export..." << std::endl;
      
      // Data export
      // set up vectors for R
      std::vector<int> df_part_id(n_pop);
      std::vector<int> df_household_id(n_pop);
      std::vector<int> df_age(n_pop);
      std::vector<int> df_status(n_pop);
      std::vector<int> df_infected_at(n_pop);
      std::vector<int> df_infected_by(n_pop);
      std::vector<int> df_infectious_at(n_pop);
      std::vector<int> df_recovered_at(n_pop);

      // add individual data to vectors
      for(int i=0; i<individuals.size(); i++){
        int status = individuals[i]->m_status;

        //do not export susceptible
        if(status > 0){
          df_part_id.emplace_back(individuals[i]->m_part_id);
          df_household_id.emplace_back(individuals[i]->m_household_id);
          df_age.emplace_back(individuals[i]->m_age);
          df_status.emplace_back(status);
          if(!model_quick){
            df_infected_by.emplace_back(individuals[i]->m_infected_by);
          }
          df_infected_at.emplace_back(individuals[i]->m_infected_at);
          df_infectious_at.emplace_back(individuals[i]->m_infectious_at);
          if(status == 5){
            df_recovered_at.emplace_back(individuals[i]->m_recovered_at);
          } else {
            df_recovered_at.emplace_back(-1.0);
          }
        }

        delete individuals[i];
      }

      std::cout << "Prepare data for R..." << std::endl;

      // make sure to copy
      Rcpp::DataFrame population;

      if(model_quick){
        population = Rcpp::DataFrame::create(
          Rcpp::Named("part_id") = df_part_id,
          Rcpp::Named("household_id") = df_household_id,
          Rcpp::Named("age") = df_age,
          Rcpp::Named("status") = df_status,
          //Rcpp::Named("infected_by") = df_infected_by,
          Rcpp::Named("infected_at") = df_infected_at,
          Rcpp::Named("infectious_at") = df_infectious_at,
          Rcpp::Named("recovered_at") = df_recovered_at
        );
      } else {
        population = Rcpp::DataFrame::create(
          Rcpp::Named("part_id") = df_part_id,
          Rcpp::Named("household_id") = df_household_id,
          Rcpp::Named("age") = df_age,
          Rcpp::Named("status") = df_status,
          Rcpp::Named("infected_by") = df_infected_by,
          Rcpp::Named("infected_at") = df_infected_at,
          Rcpp::Named("infectious_at") = df_infectious_at,
          Rcpp::Named("recovered_at") = df_recovered_at
        );
      }

      Rcpp::List intervention_data(interventions.size());
      for(int i=0; i < interventions.size(); i++){
        Rcpp::List ivlist = Rcpp::List::create(
          Rcpp::Named("start_t", iv_vals[i][0]),
          Rcpp::Named("attended", iv_vals[i][1]),
          Rcpp::Named("tested_positive") = iv_vals[i][2],
          Rcpp::Named("in_quarantine") = iv_vals[i][3],
          Rcpp::Named("susceptible") = iv_vals[i][4],
          Rcpp::Named("preinfectious") = iv_vals[i][5],
          Rcpp::Named("infectious") = iv_vals[i][6]
        );
        intervention_data[i] = ivlist;
      }

      Rcpp::List output = Rcpp::List::create(
        Rcpp::Named("intervention") = intervention_data,
        Rcpp::Named("population") = population
      );

      for(int h=0; h<n_households; h++){
        delete households[h];
      }

      std::cout << "Done" << std::endl;
      return(output);
    }

    // check for contact with all susceptible individuals
    if(model_quick){
      #pragma omp parallel for
      for(int s = 0; s < susceptibles.size(); s++){
        Randomizer* random_generator_thread = randomizers[omp_get_thread_num()];
        Individual* member = susceptibles[s];

        int infectious_strata_household[agegroups][Iclasses];
        int infectious_strata_household_quarantine[agegroups][Iclasses];
        for(int a=0; a<agegroups; a++){
          for(int i=0; i<Iclasses; i++){
            infectious_strata_household[a][i] = 0;
            infectious_strata_household_quarantine[a][i] = 0;
          }
        }

        //are any household members infectious?
        int hhmem_inf = 0;
        for(int h=0; h < member->m_household_ptr->members.size(); h++){
          Individual* hhmember = member->m_household_ptr->members[h];
          if(hhmember->m_status == Individual::infected_preclinical |
             hhmember->m_status == Individual::infected_clinical |
             hhmember->m_status == Individual::infected_subclinical
          ){
            infectious_strata_household[hhmember->m_age][hhmember->m_status - 2] += 1;
            if(hhmember->m_quarantine){
              infectious_strata_household_quarantine[hhmember->m_age][hhmember->m_status - 2] += 1;
            }
            hhmem_inf += 1;
          }
        }

        bool transmission = false;
        for(int a=0; a<agegroups; a++){
          if(transmission){break;}

          for(int i=0; i<Iclasses; i++){
            if(transmission){break;}

            int Nhh = infectious_strata_household[a][i];
            int Nhhq = infectious_strata_household_quarantine[a][i];
            bool transmission_hh = member->transmissionQuick(a, i, true, Nhh, beta, beta_hh, t_day, random_generator_thread);

            //only implement community transmission if member is not in quarantine
            if(!transmission_hh && !member->m_quarantine){
              int Ncom = infectious_strata[a][i];
              //make sure only those in the community, excluding any non-quarantining infectious household members, can infect
              bool transmission_com = member->transmissionQuick(a, i, false, Ncom - (Nhh-Nhhq), beta, beta_hh, t_day, random_generator_thread);
              transmission = transmission_com;
            } else {
              transmission = true;
            }
          }
        }
      }
    } else {
      //track who infects whom
      for(int i = 0; i < infectious.size(); i++){
        Individual* member = infectious[i];

        #pragma omp parallel for
        for(int s = 0; s < susceptibles.size(); s++){
          Randomizer* random_generator_thread = randomizers[omp_get_thread_num()];
          Individual* contact = susceptibles[s];

          bool household_member = contact->m_household_ptr == member->m_household_ptr;
          if(household_member | (!household_member & !member->m_quarantine)){
            bool transmission = contact->transmission( member, household_member, beta, beta_hh, t_day, random_generator_thread );
            
            // move contact to preinfectious vector if transmission occurs
            if(transmission){
              preinfectious.emplace_back(contact);
            }
          }
        }
      }
    }
  }

  std::cout << "Finished simulation" << std::endl;
  std::cout << "Data export..." << std::endl;
  
  // Data export
  // set up vectors for R
  std::vector<int> df_part_id(n_pop);
  std::vector<int> df_household_id(n_pop);
  std::vector<int> df_age(n_pop);
  std::vector<int> df_status(n_pop);
  std::vector<int> df_infected_at(n_pop);
  std::vector<int> df_infected_by(n_pop);
  std::vector<int> df_infectious_at(n_pop);
  std::vector<int> df_recovered_at(n_pop);

  // add individual data to vectors
  for(int i=0; i<individuals.size(); i++){
    int status = individuals[i]->m_status;

    //do not export susceptible
    if(status > 0){
      df_part_id.emplace_back(individuals[i]->m_part_id);
      df_household_id.emplace_back(individuals[i]->m_household_id);
      df_age.emplace_back(individuals[i]->m_age);
      df_status.emplace_back(status);
      if(!model_quick){
        df_infected_by.emplace_back(individuals[i]->m_infected_by);
      }
      df_infected_at.emplace_back(individuals[i]->m_infected_at);
      df_infectious_at.emplace_back(individuals[i]->m_infectious_at);
      if(status == 5){
        df_recovered_at.emplace_back(individuals[i]->m_recovered_at);
      } else {
        df_recovered_at.emplace_back(-1.0);
      }
    }

    delete individuals[i];
  }

  std::cout << "Prepare data for R..." << std::endl;

  // make sure to copy
  Rcpp::DataFrame population;

  if(model_quick){
    population = Rcpp::DataFrame::create(
      Rcpp::Named("part_id") = df_part_id,
      Rcpp::Named("household_id") = df_household_id,
      Rcpp::Named("age") = df_age,
      Rcpp::Named("status") = df_status,
      //Rcpp::Named("infected_by") = df_infected_by,
      Rcpp::Named("infected_at") = df_infected_at,
      Rcpp::Named("infectious_at") = df_infectious_at,
      Rcpp::Named("recovered_at") = df_recovered_at
    );
  } else {
    population = Rcpp::DataFrame::create(
      Rcpp::Named("part_id") = df_part_id,
      Rcpp::Named("household_id") = df_household_id,
      Rcpp::Named("age") = df_age,
      Rcpp::Named("status") = df_status,
      Rcpp::Named("infected_by") = df_infected_by,
      Rcpp::Named("infected_at") = df_infected_at,
      Rcpp::Named("infectious_at") = df_infectious_at,
      Rcpp::Named("recovered_at") = df_recovered_at
    );
  }

  Rcpp::List intervention_data(interventions.size());
  for(int i=0; i < interventions.size(); i++){
    Rcpp::List ivlist = Rcpp::List::create(
      Rcpp::Named("start_t", iv_vals[i][0]),
      Rcpp::Named("attended", iv_vals[i][1]),
      Rcpp::Named("tested_positive") = iv_vals[i][2],
      Rcpp::Named("in_quarantine") = iv_vals[i][3],
      Rcpp::Named("susceptible") = iv_vals[i][4],
      Rcpp::Named("preinfectious") = iv_vals[i][5],
      Rcpp::Named("infectious") = iv_vals[i][6]
    );
    intervention_data[i] = ivlist;
  }

  Rcpp::List output = Rcpp::List::create(
    Rcpp::Named("intervention") = intervention_data,
    Rcpp::Named("population") = population
  );

  for(int h=0; h<n_households; h++){
    delete households[h];
  }

  std::cout << "Done" << std::endl;
  return(output);
}

// [[Rcpp::export]]
Rcpp::List covidIBM(Rcpp::List parameters, int seed = 0)
{
  std::cout << "Start model" << std::endl;

  // Initialise parameters for this simulation
  Rcpp::List output;
  output = RunSimulation(parameters, seed);

  return output;
}
