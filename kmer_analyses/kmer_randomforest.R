library(randomForest)
library(tidyverse)
library(reshape2)
library(readxl)

`%nin%` = Negate(`%in%`)

oxygen <- read.table("./oxygen", header = T, sep = "\t")

# randomforest ----------------------------------------------------
data <- read_xlsx("./catal_KU.xlsx", sheet = 1)

tmp <- melt(data)

tmp <- tmp %>%
  filter(organism %nin% c("Homo sapiens", "bacterium"))%>%
  mutate(short = tolower(substr(variable,0,6))) %>%
  group_by(short) %>%
  mutate(allreads = sum(value),
         abundance = value/allreads) %>%
  arrange(-abundance, .by_group = T) %>%
  distinct(short, organism, .keep_all = T) %>%
  select(short, organism, abundance)

tmp$abundance[tmp$abundance %in% "NaN"] <- 0

temp <- tmp %>%
  pivot_wider(names_from = short, values_from = abundance, values_fill = list(abundance = 0)) %>%
  as.data.frame()

row.names(temp) <- temp$organism
temp <- temp[,-1]

category <- read_xlsx("./kmer_KU.xlsx", sheet = 2)
category <- category %>%
  distinct(library_short, .keep_all = T) %>%
  select(library_short, category) %>%
  as.data.frame()

row.names(category) <- category$library_short


otu_wide <- as.data.frame(t(temp))
otu_wide$group <- category$category[match(row.names(otu_wide), row.names(category))]
otu_wide$group <- gsub(pattern = "subadult_", replacement = "", x = otu_wide$group)
otu_wide$group <- gsub(pattern = "adult_", replacement = "", x = otu_wide$group)

otu_wide$group <- factor(otu_wide$group)


colnames(otu_wide) <- gsub(pattern = " ", replacement = "_", x = colnames(otu_wide), fixed = T)
colnames(otu_wide) <- gsub(pattern = ".", replacement = "", x = colnames(otu_wide), fixed = T)
colnames(otu_wide) <- gsub(pattern = "-", replacement = "", x = colnames(otu_wide), fixed = T)
colnames(otu_wide) <- gsub(pattern = "[", replacement = "", x = colnames(otu_wide), fixed = T)
colnames(otu_wide) <- gsub(pattern = "]", replacement = "", x = colnames(otu_wide), fixed = T)


rf_model <- randomForest(group ~ ., data = otu_wide, importance = TRUE, ntree = 10000)

# Get importance scores
importance_scores <- importance(rf_model)

# Convert importance scores to a data frame
importance_df <- as.data.frame(importance_scores)
importance_df$Microbe <- rownames(importance_df)

# Sort by MeanDecreaseGini (or MeanDecreaseAccuracy)
importance_df <- importance_df %>% arrange(desc(MeanDecreaseGini))

# Plot the top 20 most important features
#top_25_importance <- importance_df %>% top_n(100, wt = MeanDecreaseGini)
top_25_importance <- importance_df %>%
  filter(MeanDecreaseGini > 0)

top_25_importance$Microbe <- gsub(x = top_25_importance$Microbe, pattern = "_", " ")

top_25_importance <- merge(top_25_importance, oxygen, by.x = "Microbe", by.y = "TAXA_NAME", all.x = T)

ggplot(top_25_importance, aes(x = reorder(Microbe, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "#201f23") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Randomforest Top 25 Important Microbes", x = "", y = "Gini Index")+
  theme_minimal()+
  theme(strip.background = element_rect(fill = "#201f23"),
        strip.text = element_text(color = "white", face = "bold", size = 11,),
        legend.position = "none",
        axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 11, face = "bold"))

ggsave("./rf_KU_giniindex.pdf", width = 8, height = 6, device = cairo_pdf)
write.table(importance_df, "./rf_KU_importance.txt", quote = F, row.names = F, sep = "\t")

table(top_25_importance$short)

top_25_importance %>%
  filter(!is.na(short)) %>%
  ggplot(aes(x = short, y = MeanDecreaseGini)) +
  geom_boxplot(width = .35, lwd = 1)+
  geom_point(aes(fill = short), shape = 21 , color = "#201f23", size = 2, stroke = 1, 
             position = position_nudge(x = .25))+
  scale_color_manual(values = c("indianred4","navy"))+
  labs(title = "",
       subtitle = "")+
  xlab("Oxygen Tolerance")+
  ylab("Gini Impurity Score")+
  showSignificance(c(1,2), .9, -0.002, first,
                   textParams = list(size = 4, fontface = "bold"), segmentParams = list(size = 1))+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#201f23"),
        strip.text = element_text(color = "white", face = "bold", size = 11,),
        legend.position = "none",
        axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 11, face = "bold"))
  
library(superb)
sign <- function(sayı){
  x <- ifelse(sayı <= 0.05 & sayı > 0.01, "*",
              ifelse(sayı <= 0.01 & sayı > 0.001 , "**",
                     ifelse(sayı <= 0.001, "***", "n.s.")))
  return(x)
}

first <- sign(t.test(MeanDecreaseGini ~ short, 
                     data = top_25_importance[top_25_importance$short %in% c("aerobe", "anaerobe"),])$p.value)



abundance <- tmp

abundance <- merge(abundance, category, by.x = "short", by.y = "library_short", all.x = T)
abundance$category <- gsub(pattern = "subadult_", replacement = "", x = abundance$category)
abundance$category <- gsub(pattern = "adult_", replacement = "", x = abundance$category)

abundance %>%
  filter(organism %in% "Kocuria rosea") %>%
  ggplot(aes(x = category, y = abundance))+
  geom_boxplot()+
  geom_point()+
  

#  Comparison ------------------------------------------------------------
category_uniq <- read_xlsx("./kmer_KU.xlsx", sheet = 3)

prediction <- data.frame(prediction = as.character(predict(rf_model)), library = names(predict(rf_model)))
comparison <- merge(prediction, category_uniq[,c("library_short","category", "human_prop")], 
                    by.x = "library", by.y = "library_short", all.x = T)

colnames(comparison)[colnames(comparison) %in% "category"] <- "actual"


comparison <- comparison %>%
  arrange(-human_prop) %>%
  mutate(classification = ifelse(prediction == actual, "TRUE", "FALSE"))

comparison$actual[comparison$actual %in% "highprop"] <- "High Proportion"
comparison$actual[comparison$actual %in% "lowprop"] <- "Low Proportion"

comparison %>%
  ggplot(aes(x = classification, y = log10(human_prop), color = classification))+
  geom_boxplot(width = .5, lwd = 1)+
  geom_point(shape = 21, fill = "white", color = "#201f23", size = 2, stroke = 1)+
  scale_color_manual(values = c("indianred4","navy"))+
  labs(title = "Randomforest group predictions",
       subtitle = "High > 0.02 | Low < 0.0005")+
  xlab("Classification")+
  ylab(expression("log"["10"] ~ "(Human Proportion)"))+
  facet_grid(.~actual, space = "free")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#201f23"),
        strip.text = element_text(color = "white", face = "bold", size = 11,),
        legend.position = "none",
        axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 11, face = "bold"))

ggsave("./rf_KU_classification.pdf", width = 7.5, height = 5, device = cairo_pdf)
write.table(comparison, "./rf_KU_classification.txt", sep = "\t", quote = F, row.names = F)



# OBB Error Rate ------------------------------------------------------------------
ntree_values <- c(50,100, 500, 1000, 2000, 5000)
oob_errors <- numeric(length(ntree_values))

# Loop over different ntree values
for (i in seq_along(ntree_values)) {
  rf_model <- randomForest(group ~ ., data = otu_wide, ntree = ntree_values[i], importance = TRUE)
  oob_errors[i] <- rf_model$err.rate[ntree_values[i], "OOB"]
}

# Plotting the OOB error rate as a function of ntree
plot(ntree_values, oob_errors, type = "b", pch = 19, col = "blue",
     xlab = "Number of Trees", ylab = "OOB Error Rate",
     main = "Effect of ntree on OOB Error Rate")


# KMER --------------------------------------------------------------------

# randomforest ----------------------------------------------------
backup <- data
data <- fread("../rel_abundance/relabundance_5M")

tmp <- melt(data, id.vars = "kmer")
colnames(tmp) <- c("kmer", "library", "value")

tmp <- tmp %>%
  mutate(short = tolower(substr(library,0,6))) %>%
  group_by(short) %>%
  mutate(allreads = sum(value),
         abundance = value/allreads) %>%
  arrange(-abundance, .by_group = T) %>%
  distinct(short, kmer, .keep_all = T) %>%
  select(short, kmer, abundance)

tmp$abundance[tmp$abundance %in% "NaN"] <- 0


temp <- tmp %>%
  pivot_wider(names_from = short, values_from = abundance, values_fill = list(abundance = 0)) %>%
  as.data.frame()

temp <- as.data.frame(data)
row.names(temp) <- temp$kmer
temp <- temp[,-1]

category <- read_xlsx("./kmer_KU.xlsx", sheet = 2)
category <- category %>%
  distinct(library_short, .keep_all = T) %>%
  select(library_short, category) %>%
  as.data.frame()

row.names(category) <- category$library_short


otu_wide <- as.data.frame(t(temp))
otu_wide$group <- category$category[match(row.names(otu_wide), row.names(category))]
otu_wide$group <- gsub(pattern = "subadult_", replacement = "", x = otu_wide$group)
otu_wide$group <- gsub(pattern = "adult_", replacement = "", x = otu_wide$group)

otu_wide$group <- factor(otu_wide$group)
options(expressions = 5e5)
rf_model <- randomForest(x = otu_wide[, colnames(otu_wide) != "group"], 
                         y = otu_wide$group,importance = TRUE, ntree = 1000)


# Get importance scores
importance_scores <- importance(rf_model)

# Convert importance scores to a data frame
importance_df <- as.data.frame(importance_scores)
importance_df$Microbe <- rownames(importance_df)

# Sort by MeanDecreaseGini (or MeanDecreaseAccuracy)
importance_df <- importance_df %>% arrange(desc(MeanDecreaseGini))

# Plot the top 20 most important features
top_25_importance <- importance_df %>% top_n(100, wt = MeanDecreaseGini)

top_25_importance$Microbe <- gsub(x = top_25_importance$Microbe, pattern = "_", " ")

p1 <- ggplot(top_25_importance, aes(x = reorder(Microbe, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "#201f23") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Randomforest: Top 25 Important K-mers", x = "", y = "Gini Index")+
  theme_minimal()+
  theme(strip.background = element_rect(fill = "#201f23"),
        strip.text = element_text(color = "white", face = "bold", size = 11,),
        legend.position = "none",
        axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 11, face = "bold"));p1

ggsave("./rf_kmer_giniindex.pdf", width = 8, height = 6, device = cairo_pdf)
ggsave("./rf_kmer_giniindex.png", width = 8, height = 6, device = "png", dpi = 300)
write.table(importance_df, "./rf_kmer_importance.txt", quote = F, row.names = F, sep = "\t")



abundance <- tmp

abundance <- merge(abundance, category, by.x = "short", by.y = "library_short", all.x = T)
abundance$category <- gsub(pattern = "subadult_", replacement = "", x = abundance$category)
abundance$category <- gsub(pattern = "adult_", replacement = "", x = abundance$category)


abundance %>%
  filter(organism %in% "Kocuria rosea") %>%
  ggplot(aes(x = category, y = abundance))+
  geom_boxplot()+
  geom_point()
  
  
#  Comparison ------------------------------------------------------------
category_uniq <- read_xlsx(path = "./kmer_KU.xlsx", sheet = 3)

prediction <- data.frame(prediction = as.character(predict(rf_model)), library = names(predict(rf_model)))
comparison <- merge(prediction, category_uniq[,c("library_short","category", "human_prop")], 
                    by.x = "library", by.y = "library_short", all.x = T)

colnames(comparison)[colnames(comparison) %in% "category"] <- "actual"



comparison <- comparison %>%
  arrange(-human_prop) %>%
  mutate(classification = ifelse(prediction == actual, "TRUE", "FALSE"))

comparison$actual[comparison$actual %in% "highprop"] <- "High Proportion"
comparison$actual[comparison$actual %in% "lowprop"] <- "Low Proportion"

p2 <- comparison %>%
  ggplot(aes(x = classification, y = log10(human_prop), color = classification))+
  geom_boxplot(width = .5, lwd = 1)+
  geom_point(shape = 21, fill = "white", color = "#201f23", size = 2, stroke = 1)+
  scale_color_manual(values = c("indianred4","navy"))+
  labs(title = "Randomforest group predictions",
       subtitle = "5 million k-mers selected randomly out of 536.476.361 k-mers",
       caption = "High > 0.01 | Low < 0.0005")+
  xlab("Classification")+
  ylab(expression("log"["10"] ~ "(Human Proportion)"))+
  facet_grid(.~actual, space = "free")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#201f23"),
        strip.text = element_text(color = "white", face = "bold", size = 11,),
        legend.position = "none",
        axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 11, face = "bold"));p2

ggsave("./rf_kmer_classification.pdf", width = 7.5, height = 5, device = cairo_pdf)
ggsave("./rf_kmer_classification.png", width = 7.5, height = 5, device = "png", dpi = 300)
write.table(comparison, "./rf_kmer_classification.txt", sep = "\t", quote = F, row.names = F)

ggarrange(p1, p2, ncol = 1)


# OBB Error Rate ------------------------------------------------------------------
ntree_values <- c(50,100, 500, 1000, 2000, 5000)
oob_errors <- numeric(length(ntree_values))

# Loop over different ntree values
for (i in seq_along(ntree_values)) {
  rf_model <- randomForest(x = otu_wide[, colnames(otu_wide) != "group"], 
                           y = otu_wide$group,importance = TRUE,
                           ntree = ntree_values[i])
  oob_errors[i] <- rf_model$err.rate[ntree_values[i], "OOB"]
}

# Plotting the OOB error rate as a function of ntree
plot(ntree_values, oob_errors, type = "b", pch = 19, col = "blue",
     xlab = "Number of Trees", ylab = "OOB Error Rate",
     main = "Effect of ntree on OOB Error Rate")