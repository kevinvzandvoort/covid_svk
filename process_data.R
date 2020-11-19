library(data.table)

setwd("c:/workspace/slovakia/")

#https://ec.europa.eu/eurostat/databrowser/view/CENS_01RHSIZE__custom_214489/default/table?lang=en

#temp <- tempfile(fileext=".xlsx")
temp = "./eurostat_household_composition_sk_2011.xlsx"
download.file(
  url="https://ec.europa.eu/eurostat/databrowser-backend/api/query/1.0/LIVE/xlsx/en/download/1e8ac54b-ce97-477e-a3de-c4f5ea489650",
  destfile=temp, mode="ab"
)

z_agegrps <- sapply(4:22, function(x) readxl::read_xlsx(temp, range = "C7", sheet = x, col_names=F)[[1]])
z_data <- rbindlist(
  lapply(4:22, function(x, z_agegrps) as.data.table(readxl::read_xlsx(temp, range = "B11:L13", sheet = x, col_names=T))[2,][, agegroup := z_agegrps[x-3]], z_agegrps)
)
z_data <- melt(z_data, id.vars="agegroup")[!is.na(value)]
z_data[, household_size := as.numeric(gsub("([0-9]+).*$", "\\1", variable))]
agetable <- data.table(
  agegroup = z_agegrps,
  age_from = seq(0,90,5),
  age_to = c(seq(4,89,5),120)
)
household_size_age <- merge(z_data, agetable, by="agegroup")
household_size_age <- household_size_age[, c("household_size", "age_from", "age_to", "value")]
fwrite(household_size_age, "household_size_age.csv")

#use contact data agegroups
agegroups <- data.table(
  age_from = seq(0, 75 ,5),
  age_to = c(seq(4, 74, 5), 120)
)
agegroups[, name := paste0(age_from,"-",age_to)]
agegroups[age_to == 120, name := "75+"]

popdata <- fread("unwpp_slovakia_2020.csv")
household_size_age <- fread("household_size_age.csv")
contact_home <- readRDS("contact_home.RDS")
#rows = participants; cols = contacts
dimnames(contact_home) = list("participants" = agegroups[, name], "contacts" = agegroups[, name])

for(i in 1:nrow(agegroups)){
  popdata[age >= agegroups[i, age_from] & age <= agegroups[i, age_to], agegroup := factor(agegroups[i, name], agegroups$name, agegroups$name)]
  household_size_age[age_from >= agegroups[i, age_from] & age_to <= agegroups[i, age_to], agegroup := factor(agegroups[i, name], agegroups$name, agegroups$name)]
}
popdata <- popdata[, .(N=sum(popsize)), by=agegroup]
household_size_age <- dcast(household_size_age[, .(N=sum(value)), by=c("agegroup","household_size")], agegroup~household_size, value.var="N")

#assume children do not live alone
household_size_age[agegroup %in% agegroups[age_to < 18, name], "1" := 0]
household_size_prob <- rbindlist(lapply(1:nrow(household_size_age), function(x, household_size_age) household_size_age[x, -"agegroup"]/sum(household_size_age[x, -"agegroup"]), household_size_age))

#assume contacts in home are symmetric on population level
contact_home_total <- contact_home * popdata[, N]
contact_home_symmetric <- (contact_home_total + t(contact_home_total))/2/popdata[, N]
contact_home_prob <- contact_home_symmetric/rowSums(contact_home_symmetric)

target_N <- 40000

households <- matrix(0,nrow=16,ncol=target_N)
h <- 0

target_popdist <- popdata[, N]/sum(popdata[, N]) * target_N
households_rsums <- rowSums(households)

tictoc::tic()
while(sum(rowSums(households)) < target_N){
  h <- h+1
  
  #sample agegroup
  #check prob!!
  agegrp <- sample(1:nrow(popdata), 1, prob = pmax((target_popdist - households_rsums)/sum(target_popdist - households_rsums), 0))
  
  #sample household size
  hhsize <- sample(1:6, 1, prob=household_size_prob[agegrp,])
  
  #sample other household member agegroups
  members <- matrix(0,nrow=16,ncol=1)
  #ensure at least 1 adult
  while(sum(members[4:16, 1]) == 0){
    #probabilities of household member age based on at home contact matrix
    members <- rmultinom(
      1:16,
      hhsize-1,
      contact_home_prob[agegrp,] #* pmax((target_popdist - households_rsums)/sum((target_popdist - households_rsums)),0)
    )
    members[agegrp,1] <- members[agegrp,1]+1  
  }
  
  households[,h] <- members
  households_rsums <- rowSums(households)
  if(h %% 500 == 0) print(sprintf("households: %s; popsize: %s/%s",h, sum(households_rsums), target_N))
}
tictoc::toc()

households <- households[,1:h]

contact_home_symmetric

total <- matrix(0,16,16)
for(i in 1:16){
  positions <- which(households[i,] > 0)
  for(j in 1:16){
    total[i,j] <- sum(sapply(positions, function(x,i,j) return((households[j,x] - (i==j))*households[i,x]), i, j))
  }
}

plot(target_popdist, type="l")
points(households_rsums)

image(total/households_rsums)
image(contact_home_symmetric)

rowSums(total/households_rsums)
rowSums(contact_home_symmetric)

j=16
z = table(unlist(sapply(which(households[j,] >= 1), function(x) rep(sum(households[,x]),households[j,x]) )))
plot(unlist(household_size_prob[j,]), type="l")
if(length(z)==5){
  points(c(0, z/sum(z)))
}else{
  points(z/sum(z))
}

