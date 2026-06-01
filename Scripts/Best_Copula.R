


Best_Copula<-function(Con_var_1,Con_var_2){

  
  POT_Hmo  <- Con_var_1
  POT_RF <- Con_var_2


  #Axes limits for plots

  #Table of copula codes in Vine copula package
  copula_table<-data.frame(c(seq(0,40,1)[-(c(11,12,15,21,22,25,31,32,35)+1)],104,114,124,134,204,214,224,234),
                          c("Ind.","Gaussian", "t-copula", "Clayton", "Gumbel","Frank","Joe","BB1","BB6","BB7","BB8","Sur. Clayton","Sur. Gumbel","Sur. Joe",
                            "Sur. BB1","Sur. BB6","Sur. BB7","Sur. BB8","Rot. Clayton","Rot. Gumbel","Rot. Joe","Rot. BB1", "Rot. BB6",
                             "Rot. BB7","Rot. BB8","Rot. Clayton","Rot. Gumbel","Rot. Joe","Rot. BB1","Rot. BB6","Rot. BB7","Rot. BB8",
                            "Tawn type 1","Rot. Tawn type 1","Rot. Tawn type 1","Rot. Tawn type 1","Tawn type 2","Rot. Tawn type 2","Rot. Tawn type 2","Rot. Tawn type 2"))
  colnames(copula_table)<-c("Number","Family")





################### Conditional on Var1_Hmo ##########################
  
  correlation_Var1_Value<-numeric(1)
  correlation_Var1_Test<-numeric(1)
  correlation_Var1_N<-numeric(1)
  copula_Var1_Family<-numeric(1)
  copula_Var1_Family_Name<-numeric(1)

  Var1_df<-array(0,dim=c(length(POT_Hmo$Hmo),2))

  Var1_df[,1]<-POT_Hmo$Hmo
  Var1_df[,2]<-POT_Hmo$RF


  if(length(which(is.na(Var1_df[,1])==TRUE | is.na(Var1_df[,2])==TRUE))>0){
    z<-unique(which(is.na(Var1_df[,1])==TRUE | is.na(Var1_df[,2])==TRUE))
    Var1_df<-Var1_df[-z,]
    Var1_Var1_x<-Var1_Var1_x[-z]
  }
  correlation_Var1_Value<-cor(pobs(Var1_df[,1]), pobs(Var1_df[,2]),method="kendall")
  correlation_Var1_Test<-cor.test(pobs(Var1_df[,1]), pobs(Var1_df[,2]))$p.value
  correlation_Var1_N<-nrow(Var1_df)
  copula_Var1_Family<-BiCopSelect(pobs(Var1_df[,1]), pobs(Var1_df[,2]), familyset = NA, selectioncrit = "AIC",
                                  indeptest = FALSE, level = 0.05, weights = NA, rotations = TRUE,
                                  se = FALSE, presel = TRUE, method = "mle")$family
  copula_Var1_Family_Name<-as.character(copula_table$Family[which(copula_table$Number==copula_Var1_Family)])



################### Conditional on var2 RF ############################


  Var2_df<-array(0,dim=c(length(POT_RF$RF),2))
  
  Var2_df[,1]<-POT_RF$RF
  Var2_df[,2]<-POT_RF$Hmo


  if(length(which(is.na(Var2_df[,1])==TRUE | is.na(Var2_df[,2])==TRUE))>0){
    z<-unique(which(is.na(Var2_df[,1])==TRUE | is.na(Var2_df[,2])==TRUE))
    Var2_df<-Var2_df[-z,]
    Var2_Var1_x<-Var2_Var1_x[-z]
}
  correlation_Var2_Value<-cor(pobs(Var2_df[,1]), pobs(Var2_df[,2]),method="kendall")
  correlation_Var2_Test<-cor.test(pobs(Var2_df[,1]), pobs(Var2_df[,2]))$p.value
  correlation_Var2_N<-nrow(Var2_df)
  copula_Var2_Family<-BiCopSelect(pobs(Var2_df[,1]), pobs(Var2_df[,2]), familyset = NA, selectioncrit = "AIC",
                                  indeptest = FALSE, level = 0.05, weights = NA, rotations = TRUE,
                                  se = FALSE, presel = TRUE, method = "mle")$family
  copula_Var2_Family_Name<-as.character(copula_table$Family[which(copula_table$Number==copula_Var2_Family)])

  res<-list("copula_Var1_Family"=copula_Var1_Family, "copula_Var1_Family_Name"=copula_Var1_Family_Name,"copula_Var2_Family"=copula_Var2_Family, "copula_Var2_Family_Name"=copula_Var2_Family_Name)
  return(res)

}
