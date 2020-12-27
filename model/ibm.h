//forward declarations
class Household;
class Individual;
Rcpp::NumericMatrix beta;
double beta_hh;
Rcpp::NumericVector u;

/*
 * Class Household
 */
class Household{
public:
  //household id
  const unsigned int m_household_id;
  //household size
  unsigned int m_household_size;
  //pointers to individuals who are household members
  std::vector<Individual*> members;
public:
  Household(int household_id, int household_size)
    : m_household_id(household_id), m_household_size(household_size)
  {
  }
  Household(const Household& other)
    : m_household_id(other.m_household_id), m_household_size(other.m_household_size)
  {
    members.reserve(m_household_size);
    for(int i=0; i<other.m_household_size; i++)
      members.emplace_back(other.members[i]);
  }
  ~Household(){
  }
  void addMember(Individual* indiv){
    members.emplace_back(indiv);
    m_household_size++;
  }
};

/*
 * Class Individual
 */
class Individual{
public:
  //transmission status: SEIIIR
  enum Status{
    susceptible = 0,
    preinfectious = 1,
    infected_subclinical = 2,
    infected_preclinical = 3,
    infected_clinical = 4,
    recovered = 5
  };
  const unsigned int m_part_id;
  const unsigned int m_household_id;
  const unsigned int m_age;
  Status m_status = susceptible;
  double m_y;
  double m_u;
  double m_duration_E;
  double m_duration_Is;
  double m_duration_Ip;
  double m_duration_Ic;
  unsigned int m_infected_by;
  double m_infected_at;
  double m_infectious_at;
  double m_recovered_at;
  double m_lockdown_effect;
  bool m_quarantine;
  bool m_hh_quarantine;
  int m_quarantine_duration;
  Household* m_household_ptr;
  
public:
  //require IDs and age for constructor
  Individual(int part_id, int household_id, int age, Household* household_ptr, double dE, double dP, double dC, double dS, double y, double u)
    : m_part_id(part_id), m_household_id(household_id), m_age(age), m_status(susceptible), m_household_ptr(household_ptr), m_duration_E(dE), m_duration_Ip(dP), m_duration_Ic(dC), m_duration_Is(dS), m_y(y), m_u(u), m_quarantine(false), m_hh_quarantine(false), m_quarantine_duration(0), m_infected_by(-1), m_infected_at(-1), m_infectious_at(-1), m_recovered_at(-1), m_lockdown_effect(1.0)
  {
  }

  Individual(const Individual& other)
    : m_part_id(other.m_part_id), m_household_id(other.m_household_id), m_age(other.m_age), m_status(other.m_status), m_infected_by(other.m_infected_by), m_infected_at(other.m_infected_at), m_quarantine(other.m_quarantine), m_quarantine_duration(other.m_quarantine)
  {
  }
  
  ~Individual()
  {
  }
  
  //Progress to next infectious status
  bool progressStatus(Randomizer* random_generator, double time){
    bool change_status = false;
    switch(m_status){
    case preinfectious:
      if(m_duration_E == 0){
        if(random_generator->Binomial(1, m_y) == 1){ m_status = infected_preclinical; } else { m_status = infected_subclinical; }
        change_status = true;
      } else {
        m_duration_E -= 1;
      }
      break;
    case infected_subclinical:
      if(m_duration_Is == 0){ m_status = recovered; change_status = true; m_recovered_at = time; } else { m_duration_Is -= 1; }
      break;
    case infected_preclinical:
      if(m_duration_Ip == 0){ m_status = infected_clinical; change_status = true; } else { m_duration_Ip -= 1; }
      break;
    case infected_clinical:
      if(m_duration_Ic == 0){ m_status = recovered; change_status = true; m_recovered_at = time; } else { m_duration_Ic -= 1; }
      break;
    default:
      change_status = false;
    }
    if(change_status & (m_status > 1) & (m_status < 4)){
      m_infectious_at = time;
    }
    
    return change_status;
  }
  
  //transmission used in the normal implementation of the model
  //tracks who infects whom
  bool transmission(Individual* infected_contact, bool household_member, Rcpp::NumericMatrix &beta, double &beta_hh, double time, Randomizer* random_generator){
    if(m_status != susceptible){
      return false;
    } else if(!household_member & m_quarantine){
      return false;
    } else {
      // subclinical infections are 50% less infectious
      double f = 1.0;
      if(infected_contact->m_status == infected_subclinical){
        f = 0.5;
      }
      //if hh_index is 0, contacts are household members, and homogeneous mixing is used
      double lambda = m_u * f * (household_member ? beta_hh : ( beta(infected_contact->m_age, m_age)*m_lockdown_effect ));
      
      bool is_infected = random_generator->Binomial(1, lambda);
      if(is_infected){
        m_status = preinfectious;
        m_infected_at = time;
        m_infected_by = infected_contact->m_part_id;
      }
      
      return is_infected;
    }
  }
  
  //transmission used in the quick implementation of the model
  //does not track who infects whom
  bool transmissionQuick(int contacts_age, int I_class, bool hhmember, int N, Rcpp::NumericMatrix &beta, double &beta_hh, double time, Randomizer* random_generator){
    if(m_status != susceptible){
      return false;
    } else if(N == 0){
      return false;
    } else {
      // assume subclinical infections are 50% less infectious
      double f = 1.0;
      if(I_class == 0){
        f = 0.5;
      }
      
      double lambda = m_u * f * (hhmember ? beta_hh : ( beta(contacts_age, m_age)*m_lockdown_effect ));
      bool is_infected = random_generator->Binomial(1, 1-std::pow((1-lambda),N));
      if(is_infected){
        m_status = preinfectious;
        m_infected_at = time;
      }
      return is_infected;
    }
  }
  
  //progress whether individual is (still) in quarantine
  void checkQuarantine(){
    if(m_hh_quarantine & !m_quarantine){
      m_quarantine = true;
      m_hh_quarantine = false;
    }
    if(m_quarantine){
      m_quarantine_duration -= 1;
      if(m_quarantine_duration == 0){
        m_quarantine = false;
      }
    }
  }
  
  //function to start quarantine
  void startQuarantine(bool quarantine, Rcpp::List intervention, bool household_member_infected){
    int duration = intervention["duration"];
    Rcpp::NumericVector ages_eligible = intervention["age_eligible"];
    
    //only apply if individual is not already in quarantine
    if(quarantine && !m_quarantine){
      bool eligible = false;
      for(int a=0; a<ages_eligible.size(); a++){
        if(ages_eligible[a] == m_age){
          eligible = true; break;
        }
      }
      
      if(eligible){
        m_quarantine = true;
        m_quarantine_duration = duration;
      }
      
      if(household_member_infected){
        m_hh_quarantine = true;
        m_quarantine_duration = duration;
      }
    }
  }
  
  //function to check whether individual attends mass testing
  bool massTestAttend(Rcpp::List intervention, Randomizer* random_generator, double pquarantine){
    double coverage = intervention["coverage"];
    double compliance_notest = intervention["compliance_notest"];
    
    //adjust coverage to coverage in those not in quarantine
    bool attend = random_generator->Binomial(1, std::min<double>((coverage/(1-pquarantine)), 1.0));
    
    //assume those in quarantine are not retested
    if(m_quarantine){
      attend = false;
    } else if(!attend) {
      startQuarantine(
        random_generator->Binomial(1, compliance_notest),
        intervention,
        false
      );
    }
    return(attend);
  }
  
  //function to implement testing for individual
  bool massTest(Rcpp::List intervention, Rcpp::List testpars, Randomizer* random_generator, bool attend){
    
    double compliance_positive = intervention["compliance_positive"];
    double compliance_household = intervention["compliance_household"];
    double test_sensitivity_preinfectious = testpars["test_sensitivity_preinfectious"];
    double test_sensitivity_subclinical = testpars["test_sensitivity_subclinical"];
    double test_sensitivity_preclinical = testpars["test_sensitivity_preclinical"];
    double test_sensitivity_clinical = testpars["test_sensitivity_clinical"];
    double test_specificity = testpars["test_specificity"];

    bool test_positive;

    if(attend){
      double positiv_prob;
      
      switch(m_status){
      case preinfectious:
        positiv_prob = test_sensitivity_preinfectious; break;
      case infected_subclinical:
        positiv_prob = test_sensitivity_subclinical; break;
      case infected_preclinical:
        positiv_prob = test_sensitivity_preclinical; break;
      case infected_clinical:
        positiv_prob = test_sensitivity_clinical; break;
      default:
        positiv_prob = 1-test_specificity;
      }
      test_positive = random_generator->Binomial(1, positiv_prob);
      startQuarantine(
        random_generator->Binomial(1, compliance_positive*test_positive),
        intervention,
        false
      );
      if(test_positive){
        //do household members quarantine?
        for(int h=0; h < m_household_ptr->members.size(); h++){
          m_household_ptr->members[h]->startQuarantine(
            random_generator->Binomial(1, compliance_household),
            intervention,
            true
          ); 
        }
      }
    }

    return(test_positive);
  }
};
