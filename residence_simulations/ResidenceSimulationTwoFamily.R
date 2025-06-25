suppressMessages(library(dplyr))
library(pedtools)
library(ribd)

args <- commandArgs(trailingOnly = TRUE)
matrilocal_proportion <- as.numeric(args[1]) # Proportion of matrilocal buildings, e.g., 90 means 90%
num_generations <- as.numeric(args[2]) # Number of of generations, e.g., 2,3,4
num_sims <- as.numeric(args[3])  # Number of simulations


num_buildings <- 100 # Number of buildings in the village

# Load haplogroup pools from files (one haplogroup per line)

Pool_Y <- readLines("Pool_Y.txt")
Pool_mt <- readLines("Pool_mt.txt")


results_df <- data.frame(
  sim = integer(),
  AuXThetaDiffWB = numeric(),
  YHaploHomozygosityWB = numeric(),
  mtHaploHomozygosityWB = numeric(),
  ResidenceType = character(),
  stringsAsFactors = FALSE
)

# Function to generate offspring
generate_offspring <- function(building, generation, male_adult, female_adult) {
  offspring_list <- list()
  
  # Adult offspring (1 male, 1 female)
  male_offspring <- data.frame(
    ID = paste(building, generation, 1, sep = "_"),
    FID = male_adult$ID,
    MID = female_adult$ID,
    SEX = 1,
    YDNA = male_adult$YDNA,
    mtDNA = female_adult$mtDNA,
    BUILDING = building,
    GENERATION = generation,
    AGE = "A"
  )
  
  female_offspring <- data.frame(
    ID = paste(building, generation, 2, sep = "_"),
    FID = male_adult$ID,
    MID = female_adult$ID,
    SEX = 2,
    YDNA = "-",
    mtDNA = female_adult$mtDNA,
    BUILDING = building,
    GENERATION = generation,
    AGE = "A"
  )
  
  offspring_list[[length(offspring_list) + 1]] <- male_offspring
  offspring_list[[length(offspring_list) + 1]] <- female_offspring
  
  # Randomly determine the number of subadult offspring (between 2 and 4)
  num_subadult_pairs <- sample(1:2, 1)  # 1 pair (2 offspring) or 2 pairs (4 offspring)
  
  for (i in 1:num_subadult_pairs) {
    subadult_male <- data.frame(
      ID = paste(building, generation, 2 * i + 1, sep = "_"),  # Odd ID for male subadult
      FID = male_adult$ID,
      MID = female_adult$ID,
      SEX = 1,
      YDNA = male_adult$YDNA,
      mtDNA = female_adult$mtDNA,
      BUILDING = building,
      GENERATION = generation,
      AGE = "S"
    )
    
    subadult_female <- data.frame(
      ID = paste(building, generation, 2 * i + 2, sep = "_"),  # Even ID for female subadult
      FID = male_adult$ID,
      MID = female_adult$ID,
      SEX = 2,
      YDNA = "-",
      mtDNA = female_adult$mtDNA,
      BUILDING = building,
      GENERATION = generation,
      AGE = "S"
    )
    
    offspring_list[[length(offspring_list) + 1]] <- subadult_male
    offspring_list[[length(offspring_list) + 1]] <- subadult_female
  }
  
  return(offspring_list)
}

# Function to shuffle buildings based on an intermediate residence type
shuffle_buildings <- function(IndInfo, generation, matrilocal_proportion) {
  adult_males <- IndInfo[IndInfo$SEX == 1 & IndInfo$AGE == "A" & IndInfo$GENERATION == generation, ]
  adult_females <- IndInfo[IndInfo$SEX == 2 & IndInfo$AGE == "A" & IndInfo$GENERATION == generation, ]
  
  # Determine the number of matrilocal buildings based on the input proportion
  num_matrilocal_buildings <- round(matrilocal_proportion / 100 * num_buildings)
  selected_matrilocal <- sample(1:num_buildings, num_matrilocal_buildings)
  
  # Move males in matrilocal buildings
  matrilocal_males <- adult_males[adult_males$BUILDING %in% selected_matrilocal, ]
  new_male_buildings <- sample(matrilocal_males$BUILDING)
  while (any(new_male_buildings == matrilocal_males$BUILDING)) {
    new_male_buildings <- sample(matrilocal_males$BUILDING)
  }
  IndInfo[IndInfo$ID %in% matrilocal_males$ID, "BUILDING"] <- new_male_buildings
  
  # Move females in patrilocal buildings (remaining buildings)
  patrilocal_females <- adult_females[!(adult_females$BUILDING %in% selected_matrilocal), ]
  new_female_buildings <- sample(patrilocal_females$BUILDING)
  while (any(new_female_buildings == patrilocal_females$BUILDING)) {
    new_female_buildings <- sample(patrilocal_females$BUILDING)
  }
  IndInfo[IndInfo$ID %in% patrilocal_females$ID, "BUILDING"] <- new_female_buildings
  
  return(IndInfo)
}

for (sim in 1:num_sims) {
  
  IndInfo <- data.frame(
    ID = character(),
    FID = integer(),
    MID = integer(),
    SEX = integer(),
    YDNA = character(),
    mtDNA = character(),
    BUILDING = integer(),
    GENERATION = integer(),
    AGE = character(),
    stringsAsFactors = FALSE
  )
  
  # Create founders for each building
  for (building in 1:num_buildings) {
    # Male founder
    male <- data.frame(
      ID = paste(building, 0, 1, sep = "_"),
      FID = 0,
      MID = 0,
      SEX = 1,
      YDNA = sample(Pool_Y, 1),
      mtDNA = sample(Pool_mt, 1),
      BUILDING = building,
      GENERATION = 0,
      AGE = "A"
    )
    
    # Female founder
    female <- data.frame(
      ID = paste(building, 0, 2, sep = "_"),
      FID = 0,
      MID = 0,
      SEX = 2,
      YDNA = "-",
      mtDNA = sample(Pool_mt, 1),
      BUILDING = building,
      GENERATION = 0,
      AGE = "A"
    )
    
    IndInfo <- rbind(IndInfo, male, female)
  }
  
  for (generation in 1:num_generations) {
    offspring_list <- list()
    
    # Generate offspring for each building
    for (building in 1:num_buildings) {
      male_adult <- IndInfo[IndInfo$BUILDING == building & IndInfo$SEX == 1 & IndInfo$AGE == "A" & IndInfo$GENERATION == generation - 1, ]
      female_adult <- IndInfo[IndInfo$BUILDING == building & IndInfo$SEX == 2 & IndInfo$AGE == "A" & IndInfo$GENERATION == generation - 1, ]
      
      if (nrow(male_adult) > 0 && nrow(female_adult) > 0) {
        offspring <- generate_offspring(building, generation, male_adult[1, ], female_adult[1, ])
        offspring_list <- append(offspring_list, offspring)
      }
    }
    
    # Add offspring to IndInfo dataframe
    IndInfo <- rbind(IndInfo, do.call(rbind, offspring_list))
    
    # Shuffle buildings based on residence type
    IndInfo <- shuffle_buildings(IndInfo, generation, matrilocal_proportion)
  }
  
  #For Kinship Coefficients
  IndInfoFKC <- IndInfo %>% select(ID,FID,MID,SEX)
  colnames(IndInfoFKC) <- c("id","fid","mid","sex")
  
  
  IndInfoPed <- do.call(ped, list(id = IndInfoFKC$id, fid = IndInfoFKC$fid, mid = IndInfoFKC$mid, sex = IndInfoFKC$sex))
  
  IndInfoKM_Au <- kinship(IndInfoPed, Xchrom = FALSE)
  IndInfoKM_X <- kinship(IndInfoPed, Xchrom = TRUE)
  
  IndKinship_Au <- as.data.frame(as.table(IndInfoKM_Au))
  colnames(IndKinship_Au) <- c("Ind1", "Ind2", "Theta_Au")
  IndKinship_Au <- IndKinship_Au %>% mutate(Pair = paste(Ind1,Ind2,sep = "-"))
  
  IndKinship_X <- as.data.frame(as.table(IndInfoKM_X))
  colnames(IndKinship_X) <- c("Ind1", "Ind2", "Theta_X")
  IndKinship_X <- IndKinship_X %>% mutate(Pair = paste(Ind1,Ind2,sep = "-")) %>% select(Pair,Theta_X)
  
  IndKinshipPairs <- IndKinship_Au %>% left_join(IndKinship_X, by = "Pair") %>% select(-Pair)
  
  # plot(IndInfoPed, cex = 0.4, symbolsize = 0.4, margins = 1)
  
  IndInfo1 <- IndInfo
  IndInfo2 <- IndInfo
  
  colnames(IndInfo1) <- c(	
    "Ind1",
    "FID1",
    "MID1",
    "SEX1",
    "YDNA1",
    "mtDNA1",
    "BUILDING1",
    "GENERATION1",
    "AGE1")
  
  colnames(IndInfo2) <- c(	
    "Ind2",
    "FID2",
    "MID2",
    "SEX2",
    "YDNA2",
    "mtDNA2",
    "BUILDING2",
    "GENERATION2",
    "AGE2")
  
  
  AllPairResults <- IndKinshipPairs %>% left_join(IndInfo1,by = "Ind1")
  AllPairResults <- AllPairResults %>% left_join(IndInfo2,by = "Ind2")
  
  
  AllPairResults_Filt1 <- AllPairResults[AllPairResults$Ind1 != AllPairResults$Ind2, ]
  
  AllPairResults_Filt1$Pair <- apply(AllPairResults_Filt1[, c("Ind1", "Ind2")], 1, function(x) paste(sort(x), collapse = "-"))
  
  unique_pairs_data <- AllPairResults_Filt1[!duplicated(AllPairResults_Filt1$Pair), ]
  
  PairResults <- unique_pairs_data[, setdiff(names(unique_pairs_data), "Pair")]
  
  PairResults <- PairResults %>% filter(GENERATION1 != 0 & GENERATION2 != 0 )
  
  PairResults$BUILDING1 <- ifelse(PairResults$BUILDING1 %% 2 == 0, PairResults$BUILDING1 - 1, PairResults$BUILDING1)
  PairResults$BUILDING2 <- ifelse(PairResults$BUILDING2 %% 2 == 0, PairResults$BUILDING2 - 1, PairResults$BUILDING2)
  
  PairResultsWB <- PairResults %>% filter(BUILDING1 == BUILDING2)
  
  #Select 20 Buildings
  selected_buildings <- sample(unique(PairResultsWB$BUILDING1), 10)
  
  PairResultsWB <- PairResultsWB %>% filter(BUILDING1 %in% selected_buildings)
  
  # Function to sample half of the individuals, excluding founders (generation 0)
  sample_half_individuals <- function(data) {
    data %>% group_by(BUILDING1) %>%
      sample_frac(0.5)
  }
  
  PairResultsWB <- sample_half_individuals(PairResultsWB)
  
  PairResultsWB <- PairResultsWB %>% ungroup()
  
  AuXThetaDiffWB <- PairResultsWB %>% mutate(AuXThetaDiff = Theta_Au - Theta_X) %>% summarise(MeanThetaDiff = mean(AuXThetaDiff)) %>% pull(MeanThetaDiff)
  YHaploHomozygosityWB <- PairResultsWB %>% filter(SEX1 == 1 & SEX2 == 1) %>% mutate(Y_Homozygosity = if_else(YDNA1 == YDNA2,1,0)) %>% summarise(Mean_Y_Homozygosity = mean(Y_Homozygosity)) %>% pull(Mean_Y_Homozygosity)
  mtHaploHomozygosityWB <- PairResultsWB %>% mutate(mt_Homozygosity = if_else(mtDNA1 == mtDNA2,1,0)) %>% summarise(Mean_mt_Homozygosity = mean(mt_Homozygosity)) %>% pull(Mean_mt_Homozygosity)
  
  
  # Append results to dataframe
  results_df <- rbind(results_df, data.frame(
    Simulation = sim,
    AuXThetaDiffWB = AuXThetaDiffWB,
    YHaploHomozygosityWB = YHaploHomozygosityWB,
    mtHaploHomozygosityWB = mtHaploHomozygosityWB,
    MatrilocalProp = matrilocal_proportion,
    Generation = num_generations 
    
  ))
}

OutputFileName <- paste("matrilocal",matrilocal_proportion,"patrilocal",100-matrilocal_proportion,generation,sim,"twofamily.tsv",sep = ".")
write.table(results_df,OutputFileName,sep = "\t", quote = F, row.names = F)

