#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# This code is published as part of the paper Neyret et al. 2020 in  
# Agriculture, Ecosystem and the Environment (doi: )
# It aims at testing the effect of crop shift frequency on the biomass,
# density, richness and diversity of weeds in agricultural fields of 
# Northern Thailand.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

##### Load packages ####
library(readr)  
library(readxl)  
library(ggplot2)
library(BiodiversityR)
library(nlme)
library(lme4)
library(vegan)
library(car)
library(plyr)
library(lmerTest)
library(spdep)
library(MuMIn)
library(emmeans)
library(r2glmm)
library(broom)
library(reshape2)
library(indicspecies)
library(multcomp)


set.seed(101)

#### Import data ####
crop_type    = read.csv('plant_type.csv', sep = ";") 
crop_history = read.csv('crop_history.csv')
plant_com    = read.csv('plant_communities2.csv', sep = ";")

plant_com[, colnames(plant_com)[c(3, 5:139, 142:143)]] = apply(plant_com[, colnames(plant_com)[c(3, 5:139, 142:143)]],
                                                               2,
                                                               as.numeric)

# Species list
all_herbs = as.character(crop_type[crop_type$Plant.type == "Herb",]$Species_code)
all_trees_shrubs = as.character(crop_type[crop_type$Plant.type != "Herb",]$Species_code)


#### Calculate number of crop shifts and crop types in the three years preceding sampling ####
Year_season = c('X2015', 'X2016', 'X2016','X2017', 'X2017') # Year of cropping season
names(Year_season) = c('D16', 'R16', 'D17', 'R17', 'D18') # Sampling session: season (D = dry season, R = rainy season) + year

# Measure number of crops in the past three years
N_previous_crops = matrix(nrow = 20, ncol = 6) 
colnames(N_previous_crops) = c('Plot','D16', 'R16', 'D17', 'R17', 'D18') 
N_previous_crops[,'Plot'] = as.character(crop_history[,'X'])

for (i in 14:17){
  ncrops = apply(crop_history[,(i-2):i], 1, function(x){length(unique(x[!is.na(x)]))})
  nna = apply(crop_history[,(i-2):i], 1, function(x){length(x[!is.na(x)])})
  ncrops = ifelse(nna == 3, ncrops, NA)
  print(ncrops)
  N_previous_crops[,names(Year_season[Year_season==colnames(crop_history)[i]])] = ncrops
}


# Measure number of crop shifts (transitions) in the past three years
N_previous_trans = N_previous_crops
for (i in 14:17){
  ntrans = apply(crop_history[,(i-2):i], 1, function(x){
    t = 0
    for (k in 2:3){
      if (is.na(x[k-1])){t = NA}
      else {
        if (x[k] != x[k-1]){t = t+1}}}
    return(t)})
  nna = apply(crop_history[,(i-2):i], 1, function(x){length(x[!is.na(x)])})
  ntrans = ifelse(nna == 3, ntrans, NA)
  print(ntrans)
  N_previous_trans[,names(Year_season[Year_season==colnames(crop_history)[i]])] = ntrans
}

#### Plot Fig. 1 (Relative Importance index) ####
plant_com_01 = plant_com[, c(all_herbs, all_trees_shrubs)] 
plant_com_01[plant_com_01>0]=1 # make binary
rD = colSums(plant_com[, c(all_herbs, all_trees_shrubs)]) / sum(colSums(plant_com[, c(all_herbs, all_trees_shrubs)]))
rF = colSums(plant_com_01[, c(all_herbs, all_trees_shrubs)]) / nrow(plant_com_01)
RI = sort((rD + rF)/2, decreasing = T)
RI_df = data.frame(RI = RI, Species_code = names(RI))
RI_df = merge(RI_df, crop_type[, c('Species_code', 'Plant.type')])
RI_df$Species_code = factor(RI_df$Species_code, levels = rev(names(RI)))


ggRI = ggplot(RI_df[RI_df$RI > 0.01,], aes(RI, x = Species_code, fill = Plant.type)) + geom_col() + theme_bw() + 
  ylab('Relative Importance index') +
  xlab('Species') +
  coord_flip() +
  scale_fill_brewer(palette="Set2", name = 'Plant form', direction = -1
  ) + 
  scale_y_continuous(position = "right") +
 # scale_x_discrete(limits = rev(levels(RI_df[RI_df$RI > 0.01,]$Species_code))) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        #panel.grid.major.x = element_blank(),
       # axis.text.x = element_text(angle = 45,  hjust = 1, vjust = 1, size = 6),
        legend.position = "none")  +
  annotate('text', x =15, y = 0.22, label = ("+ 43 species\nwith RI < 1%")) 

ggsave(ggRI, file = 'RI.eps', height = 10, width = 4)

#### Plot Fig. 2 (crop sequences) ####
crop_history_fig = melt(data.frame('Plot' = crop_history$X,
                          'Sequence 2013-2015' =   apply(crop_history[,13:15], 1, paste, collapse = ' - '),
                          'Sequence 2014-2016' =   apply(crop_history[,14:16], 1, paste, collapse = ' - '),
                          'Sequence 2015-2017' =   apply(crop_history[,15:17], 1, paste, collapse = ' - ')), id.var = 'Plot')

# The initial dataset contained rubber tree plantations (OR) which are not used in this analysis
crop_history_fig = crop_history_fig[ crop_history_fig$value != 'OR - OR - OR',]

# Replace all codes by crop names
crop_history_fig$value = gsub('YRM', 'Maize + T',  crop_history_fig$value)
crop_history_fig$value = gsub('YRU', 'Rice + T',  crop_history_fig$value)
crop_history_fig$value = gsub('YR', 'T',  crop_history_fig$value)
crop_history_fig$value = gsub('ULRLON', 'Rice + T',  crop_history_fig$value)
crop_history_fig$value = gsub('MLON', 'Maize + T',  crop_history_fig$value)
crop_history_fig$value = gsub('LON', 'T',  crop_history_fig$value)
crop_history_fig$value = gsub('ULR', 'Rice', crop_history_fig$value)
crop_history_fig$value = gsub('Friche', 'Fallow', crop_history_fig$value)
crop_history_fig$value = gsub('NA.*', 'Unknown',   crop_history_fig$value)
crop_history_fig$value = gsub('- M -', '- Maize -',   crop_history_fig$value)
crop_history_fig$value = gsub('^M -', 'Maize -',   crop_history_fig$value)
crop_history_fig$value = gsub('- M$', '- Maize',   crop_history_fig$value)

# Classify the sequences into different sequence types
crop_history_fig$col = 'NA'
crop_history_fig[as.character(crop_history_fig$value) %in% 'Unknown',]$col = 'Unknown'
crop_history_fig[as.character(crop_history_fig$value) %in% c('Rice - Maize - Maize', 'Maize - Maize - Rice', 'Maize - Rice - Rice','Maize - Rice - Maize', 'Rice - Rice - Maize', 'Rice - Maize - Rice'),]$col = 'Rice and maize'
crop_history_fig[as.character(crop_history_fig$value) %in% c('Maize - Maize - Maize'),]$col = 'Maize only'
crop_history_fig[as.character(crop_history_fig$value) %in% c('Maize + T - Maize + T - Maize + T', 'Maize + T - Maize + T - T', 'Maize - Maize - T', 'Rice - Maize - Maize + T', 'Maize + T - T - T'),]$col = 'Trees intercropped with maize'
crop_history_fig[as.character(crop_history_fig$value) %in% c('Rice - Rice + T - Maize + T', 'Rice - Rice - T', 'Maize + T - Maize + T - Rice + T', 'Maize + T - Rice + T - T', 'Ri T-M T-T',  
                                           'Rice + T - Maize + T - T', 'Maize - Maize - Rice + T', 'Rice - Maize - Rice + T','Rice - T - T'),]$col = 'Trees intercropped with rice and maize'  
crop_history_fig[as.character(crop_history_fig$value) %in% c('Maize - Maize - Fallow', 'Rice - Rice - Fallow', 'Ri-Ri-Fallow'),]$col = 'Occurence of fallow'
crop_history_fig$value = factor(crop_history_fig$value, levels = unique(c('Unknown', 'Maize + T - Maize + T - Maize + T','Maize - Maize - Maize','Rice - Rice + T - Maize + T', 'Rice - Maize - Maize', 'Maize - Maize - Rice', 
                                                        
                                                        "Maize - Rice - Rice", 'Maize - Rice - Maize','Rice - Rice - Maize', 'Rice + T - Maize + T - T',
                                                        'Maize + T - Maize + T - Rice + T','Maize + T - Maize + T - T','Rice - Rice - T', 'Rice + T - Maize + T - T', 'Maize + T - Rice + T - T','Maize - Maize - Rice + T', 'Rice - Rice - Fallow',   
                                                        'Rice - Fallow - Fallow','Maize - Maize - Fallow' ,'Rice - Maize - Rice',  'Maize + T - T - T',
                                                        'Maize - Maize - T', 'Rice - Maize - Maize + T', 'Rice - T - T','Rice - Maize - Rice + T'
                                                        
)))

crop_history_fig$col = factor(crop_history_fig$col, levels= c('Maize only', 'Trees intercropped with maize', 'Rice and maize', 'Trees intercropped with rice and maize', 'Occurence of fallow', 'Unknown'))

sequences_labeller <- function(variable,value){
  return(sequences.names[value])
}

sequences.names  <- list(
  'Sequence.2013.2015'="Sequence 2013-2015",
  'Sequence.2014.2016'="Sequence 2014-2016",
  'Sequence.2015.2017'="Sequence 2015-2017")


crop_history_fig = crop_history_fig[as.character(crop_history_fig$value) != 'Unknown',]
crop_history_fig$variable = gsub('e\\.', 'e ',   crop_history_fig$variable )
crop_history_fig$variable = gsub('\\.', '-',   crop_history_fig$variable )

# Draw plots
crop_history_plot = ggplot(crop_history_fig, aes(x = value, fill = col)) + geom_bar() + 
  facet_wrap(~ variable, nrow = 3 ) + 
  theme_bw() +
  ylab('Number of fields') + xlab('Land use sequence') +
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.position = 'bottom', text = element_text(size=12), 
        axis.text.x = element_text(angle = 70, hjust = 1),
        plot.margin = unit(c(0.2,0.2,0.2,1), "cm")) +
  scale_fill_brewer(palette = 'Accent', direction = -1, labels = c('Maize only', 'Trees intercropped\n with maize', 
                                                                   'Rice and maize', 'Trees intercropped\n with rice and maize', 'Occurence of fallow', 'Unknown'), name = '')

ggsave(crop_history_plot, file = 'crop_history_plot.eps', width = 7)


#### Create main dataset ####
Main_data = data.frame(
  'plot' = plant_com$plot,
  'plot_id' = plant_com$plot_id,
  'Season' = plant_com$season, # Sampling session
  's' = plant_com$s, # Season (rainy/dry)
  'crop_annual' = plant_com$crop_annual,
  'crop_tree' = plant_com$crop_tree,
  'Shannon_tot' = diversity(plant_com[, c(all_herbs, all_trees_shrubs)]),
  'Shannon_herbs' = diversity(plant_com[, all_herbs]),
  'Shannon_treesshrubs' = diversity(plant_com[, all_trees_shrubs]),
  'Richness_tot' = specnumber(plant_com[, c(all_herbs, all_trees_shrubs)]),
  'Richness_herbs' = specnumber(plant_com[, all_herbs]),
  'Richness_treesshrubs' = specnumber(plant_com[, all_trees_shrubs]),
  'lAb_tot' =  log(rowSums(plant_com[, c(all_herbs, all_trees_shrubs)])),
  'sqbiomass' = sqrt(plant_com$Biomass),
  'Dominance_herbs' = diversityresult(plant_com[, all_herbs],
                                      index=c("Berger"),
                                      method=c("each site"), 
                                      sortit = FALSE, digits = 3)$Berger,
  'Dominance_tot' = diversityresult(plant_com[, c(all_herbs, all_trees_shrubs)],
                                index=c("Berger"),
                                method=c("each site"), 
                                sortit = FALSE, digits = 3)$Berger,
  'Dominance_treesshrubs' = diversityresult(plant_com[, all_trees_shrubs],
                                            index=c("Berger"),
                                            method=c("each site"), 
                                            sortit = FALSE, digits = 3)$Berger,
  'Easting' = plant_com$Easting,
  'Northing' = plant_com$Northing,
  'Ncrops' = NA,
  'Ntrans' = NA
)

# Add Ncrops and Ntrans to main dataset
for (k in 1:20){
  pl = N_previous_crops[k,1]
  print(pl)
  for (seas in colnames(N_previous_crops)[2:6]){
    print(seas)
    if(nrow(Main_data[Main_data$plot == pl & Main_data$Season == seas,])>0){
      Main_data[Main_data$plot == pl & Main_data$Season == seas & !is.na(Main_data$Season),]$Ncrops      = N_previous_crops[N_previous_crops[,'Plot'] == as.character(pl),seas]
      Main_data[Main_data$plot == pl & Main_data$Season == seas,]$Ntrans      = N_previous_trans[N_previous_trans[,'Plot'] == as.character(pl),seas]
    }}}
Main_data$Ncrops = as.factor(Main_data$Ncrops )
Main_data$Ntrans = as.factor(Main_data$Ntrans)

# Scale Easting & Northing and add some noise
# The noise is needed to run CorSpatial function later, which cannot deal with repeated measurements
Main_data$Easting_bis = scale(Main_data$Easting) + rnorm(n = nrow(Main_data), mean = 0, sd = 0.001) 
Main_data$Northing_bis = scale(Main_data$Northing) + rnorm(n = nrow(Main_data), mean = 0, sd =  0.001 )

#%%%%%%%%%%%%%%%%%%%%%#
#### Main Analysis ####
#%%%%%%%%%%%%%%%%%%%%%#

# Function for table printing
ptostar = function(x){
  s = ifelse(x > 0.1, "n.s",
             ifelse(x > 0.05,'.',
                    ifelse(x > 0.01,'*',
                           ifelse(x > 0.001,'**', '***'))))
  return(s)}

# The function results.model fits the required model with spatial autocorrelation
# and outputs model results.
# the option "maize_only" correspond to the additionnal analysis in which only plots
# with maize as the current crop are considered.

# Y = response variable
# X = explanatory variable (Ncrops or Ntrans)

results.model = function(Y, X, maize_only = 'no'){
  if (maize_only == 'yes'){
    print('Only maize plots considered')
   
     ### Create the models
   
    # If we restrict to maize, we have to remove from the analysis the levels for which we have too few points
    if (X == 'Ntrans'){
      Data = Main_data[!(is.na(Main_data$Ntrans)) & Main_data$crop_annual == 'M' & Main_data$Ntrans != 1,]
      }
    if (X == 'Ncrops'){
      Data = Main_data[!(is.na(Main_data$Ncrops)) & Main_data$crop_annual == 'M'& Main_data$Ncrops !=3,]}
    if (!(X %in% c('Ntrans','Ncrops'))){print("X should be 'Ntrans' or 'Ncrops'")}
    
    
    Data$Y = Data[, Y]
    Data$X = factor(Data[, X])
    Data = Data[!is.na(Data$Y),]
    mod_lme = lme( Y ~ X*s   + crop_tree, random =  ~ 1|plot, Data,
                   corr = corSpatial(form = ~Easting_bis + Northing_bis, type ="exponential", nugget = F), method = "ML")
    mod_lme$call$data = Data
    
    mod_lmer = lmer( Y ~ X*s   + crop_tree + (1|plot), Data)
    
    par(mfrow = c(1,2))
    plot(mod_lme$fitted, mod_lme$residuals, main = Y)
    qqnorm(mod_lme$residuals)
    qqline(mod_lme$residuals)
    
    ### Print the results
    M = data.frame(matrix(ncol = 8, nrow = 2))
    colnames(M) = c('X', 's', 'X*s', #'Spatial',
                    'Yr', 'Est (d season)','Est (r season)', 'R2_LR')
    r2 = r2beta(mod_lmer, method = 'kr')
    r = round(r2$Rsq, 2)*100
    anov = Anova(mod_lme)
    print(anov)
    
    an = anov$`Pr(>Chisq)`
    names(an) = rownames(anov)
    names(r) = r2$Effect
    M[1,1] =  paste(r['X'], ptostar(an['X']) )
    M[1,2] =  paste(r['s'], ptostar(an['s']) )
    M[1,3] =  paste(r['X:s'], ptostar(an['X:s']) )
    M[1,4] =  paste(r['crop_tree'], ptostar(an['crop_tree']) )
    EM = emmeans(mod_lme,  ~ X | s, lmer.df = "satterthwaite")
    print(cld(EM, details = TRUE))
    EM2 = cld(emmeans(mod_lme,  ~ s, lmer.df = "satterthwaite"), details = TRUE)
    print(EM2)
    estimates = paste(round(tidy(EM)$estimate, 1), ' (', round(tidy(EM)$conf.low, 1),' - ', round(tidy(EM)$conf.high, 1), ')', sep = '')
    print(estimates)
    letters = gsub(' ','',cld(EM, by = 's')$.group)
    M[1:2,5] = paste(estimates, letters)[1:2]  
    M[1:2,6] = paste(estimates, letters)[3:4] 
    M[1,7] =  r2beta(mod_lmer, partial = FALSE, 'kr')$Rsq
    return(M)
 }
  
  else{
    print('All plots considered')
    
    ### Create the models
    
    Data = Main_data[!(is.na(Main_data$Ntrans)),]
    Data$Y = Data[, Y]
    Data$X = Data[, X]
    Data = Data[!is.na(Data$Y),]
    mod_lme = lme(Y ~ X*s + crop_annual +  crop_tree , random =~ 1|plot,
                  data =  Data,    corr = corSpatial(form = ~Easting_bis + Northing_bis, type ="exponential", nugget = F))
    mod_lme$call$data = Data
    
    mod_lmer = lmer(Y ~ X*s + crop_annual +  crop_tree + (1|plot),
                    Data)
    
    par(mfrow = c(1,2))
    plot(mod_lme$fitted, mod_lme$residuals, main = Y)
    qqnorm(mod_lme$residuals)
    qqline(mod_lme$residuals)
    
    ### Print the results
    
    M = data.frame(matrix(ncol = 10, nrow = 3))
    colnames(M) = c('X', 's', 'X*s', 
                    'Annual crop','Trees', 
                    'Est (d season)','Est (r season)', 'R2_LR')
    r2 = r2beta(mod_lmer, partial = 'TRUE', 'kr')
    r = round(r2$Rsq, 2)*100
    anov = Anova(mod_lme)
    an = anov$`Pr(>Chisq)`
    names(an) = rownames(anov)
    names(r) = r2$Effect
    M[2,1] =  paste(r['X'], ptostar(an['X']) )
    M[2,2] =  paste(r['s'], ptostar(an['s']) )
    M[2,3] =  paste(r['X:s'], ptostar(an['X:s']) )
    M[2,4] =  paste(r['crop_annual'], ptostar(an['crop_annual']) )
    M[2,5] =  paste(r['crop_tree'], ptostar(an['crop_tree']) )
    EM = emmeans(mod_lme,  ~ X | s, lmer.df = "satterthwaite")
    print(cld(EM, details = TRUE))
    EM2 = cld(emmeans(mod_lme,  ~ s, lmer.df = "satterthwaite"), details = TRUE)
    print(EM2)
    estimates = paste(round(tidy(EM)$estimate, 1), ' (', round(tidy(EM)$conf.low, 1),' - ', round(tidy(EM)$conf.high, 1), ')', sep = '')
    letters = gsub(' ','',cld(EM, by = 's')$.group)
    M[,6] = paste(estimates, letters)[1:3]  
    M[,7] = paste(estimates, letters)[4:6] 
    M[2,8] = r2beta(mod_lmer, partial = FALSE, 'kr')$Rsq
    return(M)}
}



Res = rbind(
  results.model('Richness_tot', X='Ntrans',  maize_only = "no"),
  results.model('Richness_herbs', X='Ntrans',  maize_only = "no"),
  results.model('Richness_treesshrubs', X='Ntrans',  maize_only = "no"),
  results.model('Dominance_tot', X='Ntrans',  maize_only = "no"),
  results.model('Dominance_herbs', X='Ntrans',  maize_only = "no"),
  results.model('Dominance_treesshrubs', X='Ntrans',  maize_only = "no"),
  results.model('Shannon_tot', X='Ntrans',  maize_only = "no"),
  results.model('Shannon_herbs', X='Ntrans',  maize_only = "no"),
  results.model('Shannon_treesshrubs', X='Ntrans',  maize_only = "no"),
  results.model('sqbiomass', X='Ntrans',  maize_only = "no"),
  results.model('lAb_tot', X='Ntrans',  maize_only = "no")
)


mod = lm(sqbiomass ~ s + crop_tree, Main_data[Main_data$crop_annual == 'M',])
emmeans(mod, ~s)
plot(mod)
qqnorm(mod$residuals)
qqline(mod$residuals)

range(Main_data[Main_data$s != "R", 'Richness_tot'])

Anova(lm( Ab_tot.1 ~ s + Crop3 + yrwithlon, Div_trans2))

### Print Fig. 4 : Richness barplot ###
D = melt(Main_data[, c('Richness_herbs', 'Richness_treesshrubs', 'Ntrans','s', 'plot_id') ])
D = melt(Main_data[ , c('Richness_herbs', 'Richness_treesshrubs', 'Ntrans','s', 'plot_id') ])

D$variable = factor(D$variable, levels = c('Richness_treesshrubs','Richness_herbs'))
Dmean = melt(tapply(D$value, list(D$variable, D$Ntrans, D$s), mean))
Dsd = melt(tapply(D$value, list(D$variable, D$Ntrans, D$s), sd))
Dmean$sd = Dsd$value

Dmean$Var3 = mapvalues(Dmean$Var3, from = c('D', 'R'), to = c('Dry season', 'Rainy season'))
Dmean$lettre = c('a', 'a','a','ab','a', 'b','a','a','a','a','a','a')
ric_barplot = ggplot(Dmean, aes(y=value, x = as.factor(Var2), fill = Var1)) + 
  geom_col(alpha = 0.5, position = "dodge") + facet_wrap(~Var3) +theme_bw() +
  geom_errorbar(aes(ymax=value + sd, ymin=value - sd, color = Var1), position = 'dodge', fatten = 1) +
  xlab('Number of crop shifts') + ylab('Species richness') +
  geom_text(aes(y = value + sd + 1, label = lettre, group = Var1,color = Var1), 
            position =position_dodge(width = 1
            )) +
  scale_fill_brewer(palette="Set2",  breaks = c('Richness_herbs', 'Richness_treesshrubs'), labels = c('Herbaceous', 'Shrubs and trees'), name = '') +
  scale_color_brewer(palette="Dark2",breaks = c('Richness_herbs', 'Richness_treesshrubs'), labels = c('Herbaceous', 'Shrubs and trees'), name = '') +
  theme(legend.position = 'bottom') +
  guides(color = FALSE)



#### Indicator species ####

indval = multipatt(plant_com[!is.na(Main_data$Ntrans), all_herbs], 
                   Main_data[!is.na(Main_data$Ntrans),]$Ntrans,  control = how(nperm=999))
summary(indval,indvalcomp=TRUE)

indval = multipatt(plant_com_shrubs_trees[plant_com$plot_id %in% Div_trans2$plot_id, ], 
                   Div_trans2$Ntrans,  control = how(nperm=999))
summary(indval,indvalcomp=TRUE)

indval = multipatt(plant_com[!is.na(Main_data$Ntrans), c(all_herbs, all_trees_shrubs)], 
                   Main_data[!is.na(Main_data$Ntrans),]$Ntrans,  control = how(nperm=999))
summary(indval)



