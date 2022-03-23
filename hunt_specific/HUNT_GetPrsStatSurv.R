library(survival)
Infile = ""
OutFile = ""
MappingFile = ""

Input = read.table(Infile, header = T)
Map = read.table(MappingFile)
CovarCol = c("ID", "batch", "chip", "SEX", "DATE_OF_BIRTH", "DEATH", "DATE_OF_DEATH", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")

# item = "T2D"

Res = c("Item", "Target", "Beta", "SE", "P", "OR")
EndPtList = names(Input)[16:53]
EndPtList = EndPtList[EndPtList != "BMI"]

for (item in EndPtList) {
  
  index = which(Map$V1 == item)
  if (!is.na(Map[index, 2])) {
    
    print(item)
    predictor = paste("scale(", Map[index, 2], ")", sep = "")
    ColList = append(CovarCol, c(item, paste(item, "_DATE", sep = ""), Map[index, 2], "Pain_score"))
    
    Data = Input[,ColList]
    Data$DeathAge = as.numeric(as.difftime(as.Date(Data$DATE_OF_DEATH) - as.Date(Data$DATE_OF_BIRTH), units = "days")/365.25)
    
    data = Data[Data[[item]] == 1 & !is.na(Data[[item]]),ColList]
    data$DeathAge = as.numeric(as.difftime(as.Date(data$DATE_OF_DEATH) - as.Date(data$DATE_OF_BIRTH), units = "days")/365.25)
    data$DeathElapse = as.numeric(as.difftime(as.Date(data$DATE_OF_DEATH) - as.Date(data[[paste(item, "_DATE", sep = "")]]), units = "days")/365.25)
    data$DiagnosedAge = as.numeric(as.difftime(as.Date(data[[paste(item, "_DATE", sep = "")]]) - as.Date(data$DATE_OF_BIRTH), units = "days")/365.25)
    
    dropList = c("DATE_OF_BIRTH", "DATE_OF_DEATH", item, paste(item, "_DATE", sep = ""))
    
    RegressIn = Data[, !(names(Data) %in% dropList)]
    Target = "DeathInPopulation"
    form = paste("Surv(DeathAge, DEATH) ~ ", predictor, " +SEX+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+batch+chip", sep = "")
    m = coxph(as.formula(form), data = RegressIn, na.action=na.exclude)
    b = summary(m)$coefficient[predictor,"coef"]
    se = summary(m)$coefficient[predictor,"se(coef)"]
    p = summary(m)$coefficient[predictor,"Pr(>|z|)"]
    or = exp(b)
    res = c(item, Target, b, se, p, or)
    Res = rbind(Res, res)
    
    RegressIn = data[, !(names(data) %in% dropList)]
    if (length(unique(RegressIn$SEX[!is.na(RegressIn$SEX)])) == 2) {
      Target = "SurvivalTimeAfterDiagnosed(w.Age)"
      form = paste("Surv(DeathElapse, DEATH) ~ scale(", Map[index, 2], ") +SEX+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+batch+chip+DiagnosedAge", sep = "")
      m = coxph(as.formula(form), data = RegressIn, na.action=na.exclude)
      b = summary(m)$coefficient[predictor,"coef"]
      se = summary(m)$coefficient[predictor,"se(coef)"]
      p = summary(m)$coefficient[predictor,"Pr(>|z|)"]
      or = exp(b)
      res = c(item, Target, b, se, p, or)
      Res = rbind(Res, res)
      
      Target = "SurvivalTimeAfterDiagnosed(noAge)"
      form = paste("Surv(DeathElapse, DEATH) ~ scale(", Map[index, 2], ") +SEX+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+batch+chip", sep = "")
      m = coxph(as.formula(form), data = RegressIn, na.action=na.exclude)
      b = summary(m)$coefficient[predictor,"coef"]
      se = summary(m)$coefficient[predictor,"se(coef)"]
      p = summary(m)$coefficient[predictor,"Pr(>|z|)"]
      or = exp(b)
      res = c(item, Target, b, se, p, or)
      Res = rbind(Res, res)
    }
    
    else {
      Target = "SurvivalTimeAfterDiagnosed(w.Age)"
      form = paste("Surv(DeathElapse, DEATH) ~ scale(", Map[index, 2], ") +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+batch+chip+DiagnosedAge", sep = "")
      m = coxph(as.formula(form), data = RegressIn, na.action=na.exclude)
      b = summary(m)$coefficient[predictor,"coef"]
      se = summary(m)$coefficient[predictor,"se(coef)"]
      p = summary(m)$coefficient[predictor,"Pr(>|z|)"]
      or = exp(b)
      res = c(item, Target, b, se, p, or)
      Res = rbind(Res, res)
      
      Target = "SurvivalTimeAfterDiagnosed(noAge)"
      form = paste("Surv(DeathElapse, DEATH) ~ scale(", Map[index, 2], ") +PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+batch+chip", sep = "")
      m = coxph(as.formula(form), data = RegressIn, na.action=na.exclude)
      b = summary(m)$coefficient[predictor,"coef"]
      se = summary(m)$coefficient[predictor,"se(coef)"]
      p = summary(m)$coefficient[predictor,"Pr(>|z|)"]
      or = exp(b)
      res = c(item, Target, b, se, p, or)
      Res = rbind(Res, res)
    }
  }
}

write.table(Res, paste(OutFile, sep=""), row.names = F, quote = F, sep = "\t", col.names = F)
