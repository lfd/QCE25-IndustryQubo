# Copyright 2025 OTH - Laboratory for Digitalisation (LfD)
# Written by Lukas Schmidbauer
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

library(tidyverse)
library(tikzDevice)
library(stringr)
library(scales)
library(dplyr)
library(ggpubr)
library(gridExtra)
library(rlist)
library(patchwork)
library(ggh4x)
library(Cairo)

INCH.PER.CM <- 1 / 2.54
WIDTH <- 18.1 * INCH.PER.CM

COL.WIDTH <- 8.85 * INCH.PER.CM
BASE.SIZE <- 9
FORMAT <- "tex"
theme_paper_base_spacing <- function() {
    return(theme_bw(base_size = BASE.SIZE) +
        theme(
            axis.title.x = element_text(size = BASE.SIZE),
            axis.title.y = element_text(size = BASE.SIZE),
            legend.title = element_text(size = BASE.SIZE),
            legend.position = "top",
            legend.box = "vertical",
            legend.spacing.y = unit(-0.15, "cm"),
            legend.box.spacing = unit(-0.25, "cm"),
            plot.margin = unit(c(0, 0.2, 0, 0), "cm")
        ))
}

theme_paper_base <- function() {
    return(theme_bw(base_size = BASE.SIZE) +
        theme(
            axis.title.x = element_text(size = BASE.SIZE),
            axis.title.y = element_text(size = BASE.SIZE),
            legend.title = element_text(size = BASE.SIZE),
            legend.position = "top",
            legend.box = "vertical",
            plot.margin = unit(c(0, 0, 0, 0.02), "cm")
        ))
}

theme_paper_dense_spacing <- function() {
    return(theme(
            legend.text = element_text(margin=margin(t = 0, r = -0.09, b = 0, l = -0.09, unit = "cm")),

    ))
}

guide_paper_base <- function(.shape_row=1, .col_row=1) {
    return(guides(
        shape = guide_legend(order = 1, nrow = .shape_row, byrow = TRUE),
        col = guide_legend(order = 2, nrow = .col_row, byrow = TRUE, reverse = TRUE),
        linetype = guide_legend(order = 3, nrow = 1, byrow = TRUE, theme = theme(legend.text = element_text(margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"))))
    ))
}

options(
    tikzDocumentDeclaration = "\\documentclass[conference]{IEEEtran}",
    tikzLatexPackages = c(
        getOption("tikzLatexPackages"),
        "\\usepackage{amsmath}"
    ),
    tikzSanitizeCharacters = c("%", "#"),
    tikzReplacementCharacters = c("\\%", "\\#")
)

OUT.PATH <- "./../Additional_Figures/"
TIKZ.PATH <- "./tikz/"

do.save.tikz <- function(g, out.name, .width, .height, .sanitize=FALSE) {
    tikz(str_c(TIKZ.PATH, out.name, ".tex"),
        width = .width, height = .height,
        sanitize = .sanitize
    )
    print(g)
    dev.off()
}

POINT.ALPHA <- 0.6
LINE.WIDTH <- 0.4

LFD.COLOURS <- c("black", "#E69F00", "#999999", "#009371", "#beaed4", "#ed665a", "#1f78b4", "#19c716")

LFD.RAND.COL <- "#ed665a"
LFD.SA.COL <- "#1f78b4"
LFD.SHAPES <- c(15, 16, 17, 4, 5, 8, 9, 20, 25)
BOXPLOT.WIDTH <- 1

RLR.TYPE <- "solid"
RLSA.TYPE <- "11"
RLW <- 0.35

DL.TYPE <- 1

BOXPLOT.LW <- 0.5
BOXPLOT.Med <- 0.7
BOXPLOT.OUTSH <- 20
BOXPLOT.OUTAL <- 0.7

############################################################################################
# Get data from big IBM experiment
############################################################################################

relablePenalties <- function(input){
    case_when(
        #Raw
        input == "3_7" ~ "$\\lambda_m\\!\\!=\\!\\!10^3,\\lambda_t\\!\\!=\\!\\!10^7$",
        input == "4_7" ~ "$\\lambda_m\\!\\!=\\!\\!10^4,\\lambda_t\\!\\!=\\!\\!10^7$",
        input == "5_7" ~ "$\\lambda_m\\!\\!=\\!\\!10^5,\\lambda_t\\!\\!=\\!\\!10^7$",
        input == "3_8" ~ "$\\lambda_m\\!\\!=\\!\\!10^3,\\lambda_t\\!\\!=\\!\\!10^8$", # 0
        input == "4_8" ~ "$\\lambda_m\\!\\!=\\!\\!10^4,\\lambda_t\\!\\!=\\!\\!10^8$", # 2
        input == "5_8" ~ "$\\lambda_m\\!\\!=\\!\\!10^5,\\lambda_t\\!\\!=\\!\\!10^8$",        
        input == "3_9" ~ "$\\lambda_m\\!\\!=\\!\\!10^3,\\lambda_t\\!\\!=\\!\\!10^9$", # 5
        input == "4_9" ~ "$\\lambda_m\\!\\!=\\!\\!10^4,\\lambda_t\\!\\!=\\!\\!10^9$", # 6
        input == "5_9" ~ "$\\lambda_m\\!\\!=\\!\\!10^5,\\lambda_t\\!\\!=\\!\\!10^9$", # 4
        #Scaled
        input == "0_1" ~ "$\\lambda_s\\!\\!=\\!\\!1$", # 13
        input == "0_10" ~ "$\\lambda_s\\!\\!=\\!\\!0.1$", # 10
        #Rounded
        input == "0_0" ~ "rounded", #1
        #default
        TRUE ~ "NA"
    )
}

#energy
filenames <- c("results/LRQAOA/Results_raw.csv",
               "results/LRQAOA/Results_rounded.csv",
               "results/LRQAOA/Results_scaled.csv")
Ldf <- filenames %>% 
    map_df(~read_csv(.)) %>% 
    mutate(Variant = sub("_qubo.*", "", QuboName)) %>%
    mutate(Solver = "LRQAOA")

#extended with valid
filenamesV <- c("results/LRQAOA/extended_Results_raw.csv",
                "results/LRQAOA/extended_Results_rounded.csv",
                "results/LRQAOA/extended_Results_scaled.csv")
Ldfv <- filenamesV %>% 
    map_df(~read_csv(.)) %>% 
    mutate(Variant = sub("_qubo.*", "", QuboName)) %>%
    mutate(Solver = "LRQAOA")

#extended with valid and error mitigation
filenamesVE <- c("results/LRQAOA/extended_with_mitigation_Results_raw.csv",
                 "results/LRQAOA/extended_with_mitigation_Results_rounded.csv",
                 "results/LRQAOA/extended_with_mitigation_Results_scaled.csv")
Ldfve <- filenamesVE %>% 
    map_df(~read_csv(.)) %>% 
    mutate(Variant = sub("_qubo.*", "", QuboName)) %>%
    mutate(Solver = "LRQAOA")

#Random results 
filenamesRandom <- c("results/LRQAOA/extended_Results_raw_random.csv",
                     "results/LRQAOA/extended_Results_rounded_random.csv",
                     "results/LRQAOA/extended_Results_scaled_random.csv")
Ldfrand <- filenamesRandom %>% 
    map_df(~read_csv(.)) %>%
    mutate(Variant = sub("_qubo.*", "", QuboName)) %>%
    mutate(Solver = "LRQAOA")


######################
# noiseless simulation
######################

LNdfv <- read_csv("results/LRQAOA/extended_noiseless_simulation.csv") %>% 
    mutate(Solver = "Noiseless LRQAOA")
LNdfve <- read_csv("results/LRQAOA/extended_with_mitigation_noiseless_simulation.csv") %>% 
    mutate(Solver = "Noiseless LRQAOA")

############################################################################################
# Get data from big DWave experiment
############################################################################################

#energy
filenames <- c("results/Annealing/Results_raw.csv",
                "results/Annealing/Results_rounded.csv",
                "results/Annealing/Results_scaled.csv")
Adf <- filenames %>% 
    map_df(~read_csv(.)) %>%
    mutate(Solver = "Quantum Annealing")

#extended with valid
filenamesV <- c("results/Annealing/extended_Results_raw.csv",
                 "results/Annealing/extended_Results_rounded.csv",
                 "results/Annealing/extended_Results_scaled.csv")
Adfv <- filenamesV %>% 
    map_df(~read_csv(.)) %>%
    mutate(Solver = "Quantum Annealing")

#extended with valid and error mitigation
filenamesVE <- c("results/Annealing/extended_with_mitigation_Results_raw.csv",
                 "results/Annealing/extended_with_mitigation_Results_rounded.csv",
                 "results/Annealing/extended_with_mitigation_Results_scaled.csv")
Adfve <- filenamesVE %>% 
    map_df(~read_csv(.)) %>%
    mutate(Solver = "Quantum Annealing")

#Results with Embedding
filenamesEmbedding <- c("results/Annealing/Results_raw_Embedding.csv",
                        "results/Annealing/Results_rounded_Embedding.csv",
                        "results/Annealing/Results_scaled_Embedding.csv")
AdfEmbed <- filenamesEmbedding %>%
    map_df(~read_csv(.)) %>%
    mutate(Solver = "Quantum Annealing")

#Random results 
filenamesRandom <- c("results/Annealing/extended_Results_raw_random.csv",
                     "results/Annealing/extended_Results_rounded_random.csv",
                     "results/Annealing/extended_Results_scaled_random.csv")
Adfrand <- filenamesRandom %>% 
    map_df(~read_csv(.)) %>%
    mutate(Variant = sub("_qubo.*", "", QuboName)) %>%
    mutate(Solver = "Quantum Annealing")

bins <- 100


############################################################################################
# Get data from Simulated Annealing experiment
############################################################################################

#energy
filenamesSA <- c("results/Simulated/Results_raw.csv",
               "results/Simulated/Results_rounded.csv",
               "results/Simulated/Results_scaled.csv")
Sdf <- filenamesSA %>% 
    map_df(~read_csv(.)) %>%
    mutate(Solver = "Simulated")

#extended with valid
filenamesVSA <- c("results/Simulated/extended_Results_raw.csv",
                "results/Simulated/extended_Results_rounded.csv",
                "results/Simulated/extended_Results_scaled.csv")
Sdfv <- filenamesVSA %>% 
    map_df(~read_csv(.)) %>%
    mutate(Solver = "Simulated")

#extended with valid and error mitigation
filenamesVESA <- c("results/Simulated/extended_with_mitigation_Results_raw.csv",
                 "results/Simulated/extended_with_mitigation_Results_rounded.csv",
                 "results/Simulated/extended_with_mitigation_Results_scaled.csv")
Sdfve <- filenamesVESA %>% 
    map_df(~read_csv(.)) %>%
    mutate(Solver = "Simulated")

############################################################################################
# Paper plots
############################################################################################

#Reference df's for simulated annealing (1280 steps) and random (1000 shots)
SdfveRef <- Sdfve %>% 
    mutate(Toolkits = factor(Toolkits)) %>%
    filter(Annealing_time == 1280)
SdfppRefPercentValid <- SdfveRef %>% 
    group_by(Variant, Toolkits) %>%
    mutate(NumberOfSolutions = n()) %>%
    filter(Valid==TRUE) %>%
    mutate(PercentValid = n() / NumberOfSolutions) %>%
    filter(PercentValid == max(PercentValid)) %>%
    distinct(Variant, Toolkits, PercentValid)

Ldfrand <- mutate(Ldfrand, Toolkits = factor(Toolkits)) 
LdfRandPercentValid <- Ldfrand %>% 
    group_by(Variant, Toolkits) %>%
    mutate(NumberOfSolutions = n()) %>%
    filter(Valid == TRUE) %>%
    mutate(PercentValid = n() / NumberOfSolutions) %>%
    filter(PercentValid == max(PercentValid)) %>%
    distinct(Variant, Toolkits, PercentValid)

#################
Ftolerance <- 0.01

SdfppRefNumOptSolutions <- SdfveRef %>%
    filter(Valid == TRUE) %>%
    group_by(Variant, Toolkits) %>%
    mutate(NumValidSolutions = n()) %>%
    filter(abs(Cost / OptimalCost) < 1 + Ftolerance) %>%
    mutate(NumOfOptSolutions = n()) %>%
    mutate(PercentOptOfValidSolutions = NumOfOptSolutions / NumValidSolutions) %>%
    filter(PercentOptOfValidSolutions == max(PercentOptOfValidSolutions)) %>%
    distinct(Variant, Toolkits, PercentOptOfValidSolutions)

LdfppRandNumOptSolutions <- Ldfrand %>% 
    filter(Valid == TRUE) %>%
    group_by(Variant, Toolkits) %>%
    mutate(NumValidSolutions = n()) %>%
    filter(abs(Cost / OptimalCost) < 1 + Ftolerance) %>%
    mutate(NumOfOptSolutions = n()) %>%
    mutate(PercentOptOfValidSolutions = NumOfOptSolutions / NumValidSolutions) %>%
    filter(PercentOptOfValidSolutions == max(PercentOptOfValidSolutions)) %>%
    distinct(Variant, Toolkits, PercentOptOfValidSolutions)

#################
LdfppRandLowestCost <- Ldfrand %>%
    filter(Valid == TRUE) %>%
    group_by(Variant, Toolkits) %>%
    filter(Cost == min(Cost)) %>%
    mutate(Cost_OptCost = Cost / OptimalCost) %>%
    distinct(Variant, Toolkits, Cost_OptCost)

SdfppRefLowestCost <- SdfveRef %>%
    filter(Valid == TRUE) %>%
    group_by(Variant, Toolkits) %>%
    filter(Cost == min(Cost)) %>%
    mutate(Cost_OptCost = Cost / OptimalCost) %>%
    distinct(Variant, Toolkits, Cost_OptCost)

##############
# LRQAOA
##############
p.GreyScale <- c("#D3D3D3","#8C8C8C","#555555","#2A2A2A","#000000")
p.QiskitOptLevel <- 3

############################
# Noiseless simulation
############################
p.GreyScaleNoisless <- c("#D9D9D9", "#B3B3B3", "#8C8C8C", "#666666", "#4D4D4D", "#333333", "#000000")
LNdfppPercentValid <- LNdfve %>% filter(p %in% c(1,2,5,10,20,50,100)) %>%
    mutate(P1_P2 = relablePenalties(str_c(Penalty1, Penalty2, sep="_"))) %>%
    group_by(Variant, Toolkits, P1_P2, p) %>%
    mutate(NumberOfSolutions = n()) %>%
    filter(Valid == TRUE) %>%
    mutate(PercentValid = n() / NumberOfSolutions) %>%
    group_by(Variant, Toolkits, p) %>%
    filter(PercentValid == max(PercentValid)) %>%
    distinct(Variant, Toolkits, p, PercentValid, P1_P2, Solver) %>%
    mutate(Toolkits = factor(Toolkits))

LNpplotRatioValid <- ggplot(LNdfppPercentValid, aes(x=Toolkits, y=PercentValid, colour = as.factor(p), shape=as.factor(P1_P2))) +
    geom_point(alpha=POINT.ALPHA) +
    facet_grid(rows=vars(Solver), cols=vars(Variant)) +
    scale_colour_manual("$p$", values=p.GreyScaleNoisless) +
    scale_shape_manual("Penalties", values=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25)) +
    xlab("# Toolkits") +
    ylab("# Valid solutions / all solutions") +
    theme_paper_base() +
    theme_paper_dense_spacing() +
    guide_paper_base(3) 
ggsave(plot=LNpplotRatioValid, filename="NoiselessSimRatioValid.pdf", path=OUT.PATH, width=1.5*COL.WIDTH, height=COL.WIDTH)
do.save.tikz(LNpplotRatioValid, "genNoiselessSimRatioValid", 1.15 *COL.WIDTH, COL.WIDTH, TRUE)

LNdfppNumOptSolutions <- LNdfve %>% filter(Valid == TRUE) %>%
    mutate(P1_P2 = relablePenalties(str_c(Penalty1, Penalty2, sep="_"))) %>%
    group_by(Variant, Toolkits, P1_P2, p) %>%
    mutate(NumValidSolutions = n()) %>%
    filter(abs(Cost / OptimalCost) < 1 + Ftolerance) %>%
    mutate(NumOfOptSolutions = n()) %>%
    mutate(PercentOptOfValidSolutions = NumOfOptSolutions / NumValidSolutions) %>%
    group_by(Variant, Toolkits, p) %>%
    filter(PercentOptOfValidSolutions == max(PercentOptOfValidSolutions)) %>%
    distinct(Variant, Toolkits, P1_P2, p, PercentOptOfValidSolutions, Solver)

LNpplotPercentOptValid <- ggplot(LNdfppNumOptSolutions, aes(x=as.factor(Toolkits), y=PercentOptOfValidSolutions, colour = as.factor(p), shape=as.factor(P1_P2))) +
    geom_point(alpha=POINT.ALPHA) +
    facet_grid(rows=vars(Solver), cols=vars(Variant)) +
    scale_colour_manual("$p$", values = p.GreyScaleNoisless) +
    scale_shape_manual("Penalties", values=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25)) +
    xlab("# Toolkits") +
    ylab("# Optimal (1%) / # valid solutions") +
    theme_paper_base() +
    theme_paper_dense_spacing() +
    guide_paper_base(2)
ggsave(plot=LNpplotPercentOptValid, filename="NoiselessSimulationRatioOptValid.pdf", path=OUT.PATH, width=1.5*COL.WIDTH, height=COL.WIDTH)
do.save.tikz(LNpplotPercentOptValid, "genNoiselessSimRatioOptValid", 1.15 *COL.WIDTH, COL.WIDTH, TRUE)

LNdfppLowestCost <- LNdfve %>% filter(Valid == TRUE) %>%
    mutate(P1_P2 = relablePenalties(str_c(Penalty1, Penalty2, sep="_"))) %>%
    group_by(Variant, Toolkits, p) %>%
    filter(Cost == min(Cost)) %>%
    mutate(Cost_OptCost = Cost / OptimalCost) %>%
    distinct(Variant, Toolkits, p, P1_P2, Cost_OptCost, Solver)

LNpplotLowestCost <- ggplot(LNdfppLowestCost, aes(x=as.factor(Toolkits), y=Cost_OptCost, colour=as.factor(p), shape=as.factor(P1_P2))) +
    geom_point(alpha=POINT.ALPHA) +
    facet_grid(rows=vars(Solver), cols=vars(Variant)) +
    scale_colour_manual("$p$", values=p.GreyScaleNoisless) +
    scale_shape_manual("Penalties", values=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25)) +
    xlab("# Toolkits") +
    ylab("Lowest valid cost / optimal cost") +
    theme_paper_base() +
    theme_paper_dense_spacing() +
    guide_paper_base(4)
ggsave(plot=LNpplotLowestCost, filename="NoiselessSimulationLowestCost.pdf", path=OUT.PATH, width=1.5*COL.WIDTH, height=COL.WIDTH)
do.save.tikz(LNpplotLowestCost, "genNoiselessSimLowestCost", 1.15 *COL.WIDTH, COL.WIDTH, TRUE)

######
# IBM
######
LdfppPercentValid <- Ldfve %>% filter(QiskitOptLevel == p.QiskitOptLevel) %>%
    mutate(P1_P2 = relablePenalties(str_c(Penalty1, Penalty2, sep="_"))) %>%
    group_by(Variant, Toolkits, P1_P2, p) %>%
    mutate(NumberOfSolutions = n()) %>%
    filter(Valid == TRUE) %>%
    mutate(PercentValid = n() / NumberOfSolutions) %>%
    group_by(Variant, Toolkits, p) %>%
    filter(PercentValid == max(PercentValid)) %>%
    distinct(Variant, Toolkits, p, PercentValid, P1_P2, Solver) %>%
    mutate(Toolkits = factor(Toolkits)) %>%
    group_by(Variant, Toolkits) %>%
    mutate(minPercentValid = min(PercentValid), maxPercentValid = max(PercentValid))

write.csv(LdfppPercentValid, "LdfppPercentValid.csv")

LdfppNumOptSolutions <- Ldfve %>% filter(QiskitOptLevel == p.QiskitOptLevel & Valid == TRUE) %>%
    mutate(P1_P2 = relablePenalties(str_c(Penalty1, Penalty2, sep="_"))) %>%
    group_by(Variant, Toolkits, P1_P2, p) %>%
    mutate(NumValidSolutions = n()) %>%
    filter(abs(Cost / OptimalCost) < 1 + Ftolerance) %>%
    mutate(NumOfOptSolutions = n()) %>%
    mutate(PercentOptOfValidSolutions = NumOfOptSolutions / NumValidSolutions) %>%
    group_by(Variant, Toolkits, p) %>%
    filter(PercentOptOfValidSolutions == max(PercentOptOfValidSolutions)) %>%
    distinct(Variant, Toolkits, P1_P2, p, PercentOptOfValidSolutions, Solver) %>%
    mutate(Toolkits = factor(Toolkits)) %>%
    group_by(Variant, Toolkits) %>%
    mutate(minPercentOptOfValidSolutions = min(PercentOptOfValidSolutions), maxPercentOptOfValidSolutions = max(PercentOptOfValidSolutions))

write.csv(LdfppNumOptSolutions, "LdfppNumOptSolutions.csv")

LdfppLowestCost <- Ldfve %>% filter(QiskitOptLevel == p.QiskitOptLevel & Valid == TRUE) %>%
    mutate(P1_P2 = relablePenalties(str_c(Penalty1, Penalty2, sep="_"))) %>%
    group_by(Variant, Toolkits, p) %>%
    filter(Cost == min(Cost)) %>%
    mutate(Cost_OptCost = Cost / OptimalCost) %>%
    distinct(Variant, Toolkits, p, P1_P2, Cost_OptCost, Solver) %>%
    mutate(Toolkits = factor(Toolkits)) %>%
    group_by(Variant, Toolkits) %>%
    mutate(minCostOptCost = min(Cost_OptCost), maxCostOptCost = max(Cost_OptCost))

write.csv(LdfppLowestCost, "LdfppLowestCost.csv")

#Scaling plot
LdfEmbed <- Ldfve %>% 
    filter(p == 1) %>%
    mutate(QiskitOptLevel=factor(QiskitOptLevel)) %>%
    mutate(Toolkits = recode(Toolkits, "3"="3", "9"="9", "13"="13", "16"="16", "18"="", "19"="19"))

write.csv(LdfEmbed, "LdfEmbed.csv")
    
LpplotEmbeddingScaling <- ggplot(LdfEmbed, aes(x=NumQubits, y = TCNumNonLocalGates, colour=QiskitOptLevel, group=interaction(NumQubits, QiskitOptLevel))) +
    geom_boxplot(width=4*BOXPLOT.WIDTH, outlier.shape = BOXPLOT.OUTSH, outlier.alpha = BOXPLOT.OUTAL) +
    facet_grid(rows = vars(Solver), cols = vars(Variant)) +
    scale_y_continuous(labels = label_comma(big.mark = "'")) +
    scale_x_continuous(breaks=unique(LdfEmbed$NumQubits), labels=unique(LdfEmbed$Toolkits)) +
    scale_colour_manual("Qiskit optimization level", values = LFD.COLOURS) +
    xlab("# Toolkits") +
    ylab("# Transpiled non-local gates ($p=1$)") +
    theme_paper_base_spacing() +
    guide_paper_base()
ggsave(plot=LpplotEmbeddingScaling, filename="IBMEmbeddigScaling.pdf", path=OUT.PATH, width=2*COL.WIDTH, height=COL.WIDTH)

#Influence of Annealing Time
LdfLayers <- Ldf %>% filter(QiskitOptLevel == p.QiskitOptLevel & Toolkits == 19) %>%
    mutate(P1_P2 = relablePenalties(str_c(Penalty1, Penalty2, sep="_"))) %>%
    filter(P1_P2 %in% c("$\\lambda_m=10^3,\\lambda_t=10^9$", "rounded", "$\\lambda_s=1$")) %>%
    mutate(NormalisedEnergy = Energy / max(abs(Energy)), p=factor(p))

LpplotLayers <- ggplot(LdfLayers, aes(x=p, y=NormalisedEnergy)) +
    geom_boxplot(outlier.shape = BOXPLOT.OUTSH, outlier.alpha = BOXPLOT.OUTAL) +
    ggh4x::facet_grid2(rows=vars(Solver), cols=vars(Variant), scales = "free_y", independent = "y") +
    xlab("$p$") +
    ylab("Normalised energy") +
    theme_paper_base_spacing() +
    guide_paper_base()
ggsave(plot=LpplotLayers, filename="IBMAnnealingTime.pdf", path=OUT.PATH, width=2*COL.WIDTH, height=COL.WIDTH)

##############
# Annealing
##############
AnnealingTime.GreyScale <- c("#D3D3D3","#B0B0B0","#8C8C8C","#707070","#555555","#404040","#2A2A2A","#000000")
p.Schedule <- "linear"

AdfppPercentValid <- Adfve %>% filter(Schedule == p.Schedule) %>%
    mutate(P1_P2 = relablePenalties(str_c(Penalty1, Penalty2, sep="_"))) %>%
    group_by(Variant, Toolkits, P1_P2, Annealing_time) %>% 
    mutate(NumberOfSolutions = n()) %>%
    filter(Valid == TRUE) %>%
    mutate(PercentValid = n() / NumberOfSolutions) %>%
    group_by(Variant, Toolkits, Annealing_time) %>% 
    filter(PercentValid == max(PercentValid)) %>%
    distinct(Variant, Toolkits, P1_P2, Annealing_time, PercentValid, Solver) %>%
    mutate(Toolkits = factor(Toolkits)) %>%
    group_by(Variant, Toolkits) %>%
    mutate(minPercentValid = min(PercentValid), maxPercentValid = max(PercentValid))

write.csv(AdfppPercentValid, "AdfppPercentValid.csv")

AdfppNumOptSolutions <- Adfve %>% filter(Schedule == p.Schedule & Valid == TRUE) %>%
    mutate(P1_P2 = relablePenalties(str_c(Penalty1, Penalty2, sep="_"))) %>%
    group_by(Variant, Toolkits, P1_P2, Annealing_time) %>%
    mutate(NumValidSolutions = n()) %>%
    filter(abs(Cost / OptimalCost) < 1 + Ftolerance) %>%
    mutate(NumOfOptSolutions = n()) %>%
    mutate(PercentOptOfValidSolutions = NumOfOptSolutions / NumValidSolutions) %>%
    group_by(Variant, Toolkits, Annealing_time) %>% 
    filter(PercentOptOfValidSolutions == max(PercentOptOfValidSolutions)) %>%
    distinct(Variant, Toolkits, P1_P2, Annealing_time, PercentOptOfValidSolutions, Solver) %>%
    mutate(Toolkits = factor(Toolkits)) %>%
    group_by(Variant, Toolkits) %>%
    mutate(minPercentOptOfValidSolutions = min(PercentOptOfValidSolutions), maxPercentOptOfValidSolutions = max(PercentOptOfValidSolutions))

write.csv(AdfppNumOptSolutions, "AdfppNumOptSolutions.csv")

AdfppLowestCost <- Adfve %>% filter(Schedule == p.Schedule & Valid == TRUE) %>%
    mutate(P1_P2 = relablePenalties(str_c(Penalty1, Penalty2, sep="_"))) %>%
    group_by(Variant, Toolkits, Annealing_time) %>%
    filter(Cost == min(Cost)) %>%
    mutate(Cost_OptCost = Cost / OptimalCost) %>%
    distinct(Variant, Toolkits, Annealing_time, P1_P2, Cost_OptCost, Solver) %>%
    mutate(Toolkits=factor(Toolkits)) %>%
    group_by(Variant, Toolkits) %>%
    mutate(minCostOptCost = min(Cost_OptCost), maxCostOptCost = max(Cost_OptCost))

write.csv(AdfppLowestCost, "AdfppLowestCost.csv")

#Embedding Scaling
AdfEmbed <- AdfEmbed %>%
    mutate(Toolkits = recode(Toolkits, "3"="3", "9"="9", "13"="13", "16"="16", "18"="", "19"="19"))

write.csv(AdfEmbed, "AdfEmbed.csv")

ApplotEmbeddingScaling <- ggplot(AdfEmbed, aes(x=NumQubits, y=EmbeddingSize, group=as.factor(NumQubits))) +
    geom_boxplot(width=4*BOXPLOT.WIDTH, outlier.shape = BOXPLOT.OUTSH, outlier.alpha = BOXPLOT.OUTAL) +
    facet_grid(rows=vars(Solver), cols=vars(Variant)) +
    scale_y_continuous(labels = label_comma(big.mark = "'")) +
    scale_x_continuous(breaks=unique(AdfEmbed$NumQubits), labels=unique(AdfEmbed$Toolkits)) +
    xlab("# Toolkits") +
    ylab("Embedding size") +
    theme_paper_base_spacing() +
    guide_paper_base() 
ggsave(plot=ApplotEmbeddingScaling, filename="DWaveEmbeddingSize.pdf", path=OUT.PATH, width=2*COL.WIDTH, height=COL.WIDTH)


#Influence of Annealing Time
AdfAnnealingTime <- Adf %>% filter(Schedule == p.Schedule & Toolkits == 19) %>%
    mutate(P1_P2 = relablePenalties(str_c(Penalty1, Penalty2, sep="_"))) %>%
    filter(P1_P2 %in% c("$\\lambda_m=10^3,\\lambda_t=10^9$", "rounded", "$\\lambda_s=1$")) %>%
    mutate(NormalisedEnergy = energy / max(abs(energy)), Annealing_time=factor(Annealing_time))

write.csv(AdfAnnealingTime, "AdfAnnealingTime.csv")

ApplotAnnealingTime <- ggplot(AdfAnnealingTime, aes(x=Annealing_time, y=NormalisedEnergy)) +
    geom_boxplot(outlier.shape = BOXPLOT.OUTSH, outlier.alpha = BOXPLOT.OUTAL) +
    ggh4x::facet_grid2(rows=vars(Solver), cols=vars(Variant), scales = "free_y", independent = "y") +
    xlab("Annealing time") +
    ylab("Normalised energy") +
    theme_paper_base_spacing() +
    guide_paper_base()
ggsave(plot=ApplotAnnealingTime, filename="DWaveAnnealingTime.pdf", path=OUT.PATH, width=2*COL.WIDTH, height=COL.WIDTH)


###########################
# combined IBM / Annealing
###########################

pPlot.Width <- 1
pPlot.Height <- 1.6

CompplotEmbeddingScaling <- ApplotEmbeddingScaling / LpplotEmbeddingScaling
ggsave(plot=CompplotEmbeddingScaling, filename="pEmbeddingScaling.pdf", path=OUT.PATH, width=2*COL.WIDTH, height=2*COL.WIDTH)
do.save.tikz(CompplotEmbeddingScaling, "genEmbeddingScaling", pPlot.Width*COL.WIDTH, 0.8*pPlot.Height*COL.WIDTH, TRUE)

CompplotAnnealingTimeP <- ApplotAnnealingTime / LpplotLayers
ggsave(plot=CompplotAnnealingTimeP, filename="pAnnealingTimeLayers.pdf", path=OUT.PATH, width=2*COL.WIDTH, height=2*COL.WIDTH)
do.save.tikz(CompplotAnnealingTimeP, "genAnnealingTimeLayers", pPlot.Width*COL.WIDTH, 0.8*pPlot.Height*COL.WIDTH, TRUE)

#########################################
# Combined IBM / Annealing
#########################################
RandrefName <- "Random reference"
SArefName <- "Simulated Annealing reference"
dfppPercentValid <- bind_rows(rename(LdfppPercentValid, SolverParam = p), rename(AdfppPercentValid, SolverParam = Annealing_time))
dfppPercentValid$Solver_f <- factor(dfppPercentValid$Solver, levels=c("Quantum Annealing", "LRQAOA"))
dfppPercentValid$P1_P2_f <- factor(dfppPercentValid$P1_P2, levels=relablePenalties(c("3_9", "4_9", "0_1", "0_0")))
PercentValidShapes <- c(5,6,13,1)
pplotRatioValid <- ggplot(dfppPercentValid, aes(x=Toolkits, y=PercentValid)) +
    facet_grid(rows=vars(Solver_f), cols=vars(Variant)) +
    geom_boxplot(linewidth = BOXPLOT.LW, fatten = BOXPLOT.Med, outlier.shape = BOXPLOT.OUTSH, outlier.alpha = BOXPLOT.OUTAL) +
    geom_point(data=dfppPercentValid %>% filter(maxPercentValid==PercentValid & Toolkits %in% c(13,16,18,19)), aes(y=maxPercentValid, shape=P1_P2), position=position_dodge(width=1)) +
    scale_shape_manual("Penalties", values=PercentValidShapes) +
    xlab("# Toolkits") +
    ylab("# Valid solutions / all solutions") +
    geom_segment(data=LdfRandPercentValid, aes(x=as.numeric(Toolkits)-RLW, xend = as.numeric(Toolkits)+RLW, y=PercentValid, yend=PercentValid, linetype=RandrefName), colour=LFD.RAND.COL, inherit.aes=FALSE) +
    geom_segment(data=SdfppRefPercentValid, aes(x=as.numeric(Toolkits)-RLW, xend = as.numeric(Toolkits)+RLW, y=PercentValid, yend=PercentValid, linetype=SArefName), colour=LFD.SA.COL, inherit.aes=FALSE) +
    scale_linetype_manual("",values=c(RLR.TYPE,RLSA.TYPE)) +
    theme_paper_base() +
    theme_paper_dense_spacing() +
    guide_paper_base(2) 
ggsave(plot=pplotRatioValid, filename="pplotPercentValid.pdf", path=OUT.PATH, width=2*COL.WIDTH, height=2*COL.WIDTH)
do.save.tikz(pplotRatioValid, "genPercentValid", pPlot.Width*COL.WIDTH, 0.8*pPlot.Height*COL.WIDTH, TRUE)

dfppPercentOptValid <- bind_rows(rename(LdfppNumOptSolutions, SolverParam = p), rename(AdfppNumOptSolutions, SolverParam = Annealing_time))
dfppPercentOptValid$Solver_f <- factor(dfppPercentOptValid$Solver, levels=c("Quantum Annealing", "LRQAOA"))
dfppPercentOptValid$P1_P2_f <- factor(dfppPercentOptValid$P1_P2, levels=relablePenalties(c("3_8", "4_8", "4_9", "5_9", "0_10", "0_1", "0_0")))
PercentOptValidShapes <- c(0,2,6,4,10,13,1)
ppPercentOptValid <- ggplot(dfppPercentOptValid, aes(x=Toolkits, y=PercentOptOfValidSolutions)) +
    facet_grid(rows=vars(Solver_f), cols=vars(Variant)) +
    geom_boxplot(linewidth = BOXPLOT.LW, fatten = BOXPLOT.Med, outlier.shape = BOXPLOT.OUTSH, outlier.alpha = BOXPLOT.OUTAL) +
    geom_point(data=dfppPercentOptValid %>% filter(maxPercentOptOfValidSolutions==PercentOptOfValidSolutions & Toolkits %in% c(13,16,18,19)), aes(y=maxPercentOptOfValidSolutions, shape=P1_P2_f), position=position_dodge(width=1)) +
    scale_shape_manual("Penalties", values=PercentOptValidShapes) +
    xlab("# Toolkits") +
    ylab("# Optimal (1%) / # valid solutions") +
    geom_segment(data=LdfppRandNumOptSolutions, aes(x=as.numeric(Toolkits)-RLW, xend = as.numeric(Toolkits)+RLW, y=PercentOptOfValidSolutions, yend=PercentOptOfValidSolutions, linetype=RandrefName), colour=LFD.RAND.COL, inherit.aes=FALSE) +
    geom_segment(data=SdfppRefNumOptSolutions, aes(x=as.numeric(Toolkits)-RLW, xend = as.numeric(Toolkits)+RLW, y=PercentOptOfValidSolutions, yend=PercentOptOfValidSolutions, linetype=SArefName), colour=LFD.SA.COL, inherit.aes=FALSE) +
    scale_linetype_manual("", values=c(RLR.TYPE,RLSA.TYPE)) +
    theme_paper_base() +
    theme_paper_dense_spacing() +
    guide_paper_base(3) 
ggsave(plot=ppPercentOptValid, filename="pplotPercentOptValid.pdf", path=OUT.PATH, width=2*COL.WIDTH, height=2*COL.WIDTH)
do.save.tikz(ppPercentOptValid, "genPercentOptValid", pPlot.Width*COL.WIDTH, 0.9*pPlot.Height*COL.WIDTH, TRUE)

dfppLowestCost <- bind_rows(rename(LdfppLowestCost, SolverParam = p), rename(AdfppLowestCost, SolverParam = Annealing_time))
dfppLowestCost$Solver_f <- factor(dfppLowestCost$Solver, levels=c("Quantum Annealing", "LRQAOA"))
dfppLowestCost$P1_P2_f <- factor(dfppLowestCost$P1_P2, levels=relablePenalties(c("3_8", "3_9", "4_9", "0_10", "0_1", "0_0")))
LowestCostShapes <- c(0,5,6,10,13,1)
pplotLowestCost <- ggplot(dfppLowestCost, aes(x=Toolkits, y=Cost_OptCost)) +
    facet_grid(rows=vars(Solver_f), cols=vars(Variant)) +
    geom_boxplot(linewidth = BOXPLOT.LW, fatten = BOXPLOT.Med, outlier.shape = BOXPLOT.OUTSH, outlier.alpha = BOXPLOT.OUTAL) +
    geom_point(data=dfppLowestCost %>% filter(minCostOptCost==Cost_OptCost & Toolkits %in% c(13,16,18,19)), aes(y=minCostOptCost, shape=P1_P2_f), position=position_dodge(width=1)) +
    scale_shape_manual("Penalties", values=LowestCostShapes) +
    xlab("# Toolkits") +
    ylab("Lowest valid cost / optimal cost") +
    geom_segment(data=LdfppRandLowestCost, aes(x=as.numeric(Toolkits)-RLW, xend = as.numeric(Toolkits)+RLW, y=Cost_OptCost, yend=Cost_OptCost, linetype=RandrefName), colour=LFD.RAND.COL, inherit.aes=FALSE) +
    geom_segment(data=SdfppRefLowestCost, aes(x=as.numeric(Toolkits)-RLW, xend = as.numeric(Toolkits)+RLW, y=Cost_OptCost, yend=Cost_OptCost, linetype=SArefName), colour=LFD.SA.COL, inherit.aes=FALSE) +
    scale_linetype_manual("", values=c(RLR.TYPE,RLSA.TYPE)) +
    theme_paper_base() +
    theme_paper_dense_spacing() +
    guide_paper_base(2)
ggsave(plot=pplotLowestCost, filename="pplotLowestCost.pdf", path=OUT.PATH, width=2*COL.WIDTH, height=2*COL.WIDTH)
do.save.tikz(pplotLowestCost, "genLowestCost", pPlot.Width*COL.WIDTH, 0.8*pPlot.Height*COL.WIDTH, TRUE)

############################################
# Compute Percent valid solution correlation
############################################
AdfppPercentValidCorr <- AdfppPercentValid %>%
    group_by(Variant, Toolkits) %>%
    summarise(MeanPercentValid = mean(PercentValid), .groups="drop")
write.csv(AdfppPercentValidCorr, "ValidCorrMean.csv")
Ahlines <- AdfppPercentValidCorr %>% 
    group_by(Variant) %>%
    summarise(Avg = mean(MeanPercentValid))
Shlines <- SdfppRefPercentValid %>% 
    group_by(Variant) %>%
    summarise(Avg = mean(PercentValid))
PercentValidCorrPlot <- ggplot(AdfppPercentValidCorr, aes(x=as.factor(Toolkits), y = MeanPercentValid)) +
    geom_point() +
    facet_grid(cols=vars(Variant)) +
    geom_hline(data=Ahlines, aes(yintercept = Avg)) +
    geom_segment(data=SdfppRefPercentValid, aes(x=as.numeric(Toolkits)-RLW, xend = as.numeric(Toolkits)+RLW, y=PercentValid, yend=PercentValid, linetype="Simulated Annealing"), colour=LFD.SA.COL, inherit.aes=FALSE) +
    geom_hline(data=Shlines, aes(yintercept = Avg), colour="blue", linetype = "11") +
    theme_paper_base()
ggsave(plot=PercentValidCorrPlot, filename="pPercentValidCorr.pdf", path=OUT.PATH, width=2*COL.WIDTH, height=COL.WIDTH)
AdfppPercentValidCorr <- AdfppPercentValidCorr %>%
    left_join(SdfppRefPercentValid, by = c("Variant", "Toolkits")) %>%
    group_by(Variant) %>%
    summarise(Corr = cor(MeanPercentValid, PercentValid))
write.csv(AdfppPercentValidCorr, "ValidCorr.csv")

AdfppNumOptSolutions$Toolkits <- factor(AdfppNumOptSolutions$Toolkits)
SdfppRefNumOptSolutions$Toolkits <- factor(SdfppRefNumOptSolutions$Toolkits)
AdfppPercentOptValidCorr <- AdfppNumOptSolutions %>%
    group_by(Variant, Toolkits) %>%
    summarise(MeanPercentOptValid = mean(PercentOptOfValidSolutions), .groups="drop")
write.csv(AdfppPercentOptValidCorr, "OptValidCorrMean.csv")
Ahlines <- AdfppPercentOptValidCorr %>% 
    group_by(Variant) %>%
    summarise(Avg = mean(MeanPercentOptValid))
Shlines <- SdfppRefNumOptSolutions %>% 
    group_by(Variant) %>%
    summarise(Avg = mean(PercentOptOfValidSolutions))
PercentOptValidCorrPlot <- ggplot(AdfppPercentOptValidCorr, aes(x=Toolkits, y = MeanPercentOptValid)) +
    geom_point() +
    geom_hline(data=Ahlines, aes(yintercept = Avg)) +
    geom_segment(data=SdfppRefNumOptSolutions, aes(x=as.numeric(Toolkits)-RLW, xend = as.numeric(Toolkits)+RLW, y=PercentOptOfValidSolutions, yend=PercentOptOfValidSolutions, linetype="Simulated Annealing"), colour=LFD.SA.COL, inherit.aes=FALSE) +
    geom_hline(data=Shlines, aes(yintercept = Avg), colour="blue", linetype = "11") +
    facet_grid(cols=vars(Variant)) +
    theme_paper_base()
ggsave(plot=PercentOptValidCorrPlot, filename="pPercentOptValidCorr.pdf", path=OUT.PATH, width=2*COL.WIDTH, height=COL.WIDTH)
AdfppPercentOptValidCorr <- AdfppPercentOptValidCorr %>% 
    left_join(SdfppRefNumOptSolutions, by = c("Variant", "Toolkits")) %>% 
    group_by(Variant) %>% 
    summarise(Corr = cor(MeanPercentOptValid, PercentOptOfValidSolutions))
write.csv(AdfppPercentOptValidCorr, "OptValidCorr.csv")

AdfppLowestCost$Toolkits <- factor(AdfppLowestCost$Toolkits)
SdfppRefLowestCost$Toolkits <- factor(SdfppRefLowestCost$Toolkits)
AdfppLowestCostCorr <- AdfppLowestCost %>%
    group_by(Variant, Toolkits) %>%
    summarise(MeanLowestCost = mean(Cost_OptCost), .groups="drop")
write.csv(AdfppLowestCostCorr, "LowestCostCorrMean.csv")
Ahlines <- AdfppLowestCostCorr %>% 
    group_by(Variant) %>%
    summarise(Avg = mean(MeanLowestCost))
Shlines <- SdfppRefLowestCost %>% 
    group_by(Variant) %>%
    summarise(Avg = mean(Cost_OptCost))
PercentLowestCostCorrPlot <- ggplot(AdfppLowestCostCorr, aes(x=Toolkits, y = MeanLowestCost)) +
    geom_point() +
    geom_hline(data=Ahlines, aes(yintercept = Avg)) +
    geom_segment(data=SdfppRefLowestCost, aes(x=as.numeric(Toolkits)-RLW, xend = as.numeric(Toolkits)+RLW, y=Cost_OptCost, yend=Cost_OptCost, linetype="Simulated Annealing"), colour=LFD.SA.COL, inherit.aes=FALSE) +
    geom_hline(data=Shlines, aes(yintercept = Avg), colour="blue", linetype = "11") +
    facet_grid(cols=vars(Variant)) +
    theme_paper_base()
ggsave(plot=PercentLowestCostCorrPlot, filename="pLowestCostCorr.pdf", path=OUT.PATH, width=2*COL.WIDTH, height=COL.WIDTH)
AdfppLowestCostCorr <- AdfppLowestCostCorr %>%
    left_join(SdfppRefLowestCost, by = c("Variant", "Toolkits")) %>%
    group_by(Variant) %>%
    summarise(Corr = cor(MeanLowestCost, Cost_OptCost))
write.csv(AdfppLowestCostCorr, "LowestCostCorr.csv")

############################################################################################
# All plots
############################################################################################

#relable penalties
relableP1 <- function(P1){
    case_when(
        #Raw
        P1 == "3" ~ "λ_m = 10^3",
        P1 == "4" ~ "λ_m = 10^4",
        P1 == "5" ~ "λ_m = 10^5",
        #Scaled
        P1 == "0" ~ "",
        P1 == "0" ~ "",
        #Rounded
        P1 == "0" ~ "",
        #default
        TRUE ~ "NA"
    )
}
relableP2 <- function(P2){
    case_when(
        #Raw
        P2 == "7" ~ "λ_t = 10^7",
        P2 == "8" ~ "λ_t = 10^8",
        P2 == "9" ~ "λ_t = 10^9",
        #Scaled
        P2 == "1" ~ "λ_s = 1",
        P2 == "10" ~ "λ_s = 0.1",
        #Rounded
        P2 == "0" ~ "rounded",
        #default
        TRUE ~ "NA"
    )
}

plotVecs <- list()
for (TK in c(3,9,13,16,18,19)){
    plotVecs <- list.append(plotVecs, list(TK=TK, name="raw"    ))
    plotVecs <- list.append(plotVecs, list(TK=TK, name="rounded"))
    plotVecs <- list.append(plotVecs, list(TK=TK, name="scaled" ))
}

pos <- position_jitter(h = 0)

Ldf <- Ldf   %>% mutate(Penalty1 = relableP1(Penalty1), Penalty2 = relableP2(Penalty2))
Ldfv <- Ldfv  %>% mutate(Penalty1 = relableP1(Penalty1), Penalty2 = relableP2(Penalty2))
Ldfve <- Ldfve %>% mutate(Penalty1 = relableP1(Penalty1), Penalty2 = relableP2(Penalty2))

#.dataVec: (Toolkits, name)
L_multi_plot_energy <- function(.dataVec) {
    if(nrow(filter(Ldf, Toolkits == .dataVec$TK & Variant == .dataVec$name)) > 0){
        pltResEnergyLandscape <- ggplot(filter(Ldf, Toolkits == .dataVec$TK & Variant == .dataVec$name & QiskitOptLevel == 3), aes(x=Energy, colour = as.factor(p))) +
        geom_freqpoly(position = pos, bins = bins) + 
        facet_grid(cols=vars(Penalty1), rows=vars(Penalty2)) +
        scale_colour_manual("# LR-QAOA Layers", values = LFD.COLOURS) +
        ylab("Count") +
        theme_paper_base()

        ggsave(plot = pltResEnergyLandscape, filename = paste("LRQAOAResults", .dataVec$name, "EnergyLandscape", .dataVec$TK, "TK.pdf", sep=""), path=paste(OUT.PATH, .dataVec$TK, "_Toolkits/", sep=""), device=cairo_pdf)
    }
}

#.dataVec: (Toolkits, name)
L_multi_plot_valid <- function(.dataVec) {
    if(nrow(filter(Ldfv, Toolkits == .dataVec$TK & Variant == .dataVec$name & Valid == TRUE)) > 0){
        pltResEnergyLandscape <- ggplot(filter(Ldfv, Toolkits == .dataVec$TK & Variant == .dataVec$name & Valid == TRUE & QiskitOptLevel == 3), aes(x=Cost, colour = as.factor(p))) +
        geom_freqpoly(position = pos, bins = bins) + 
        facet_grid(cols=vars(Penalty1), rows=vars(Penalty2)) +
        geom_vline(aes(xintercept=OptimalCost)) + 
        scale_colour_manual("# LR-QAOA Layers", values = LFD.COLOURS) +
        ylab("Count") +
        theme_paper_base()

        ggsave(plot = pltResEnergyLandscape, filename = paste("LRQAOAResults", .dataVec$name, "ValidCostLandscape", .dataVec$TK, "TK.pdf", sep=""), path=paste(OUT.PATH, .dataVec$TK, "_Toolkits/", sep=""), device=cairo_pdf)
    }
}

#.dataVec: (Toolkits, name)
L_multi_plot_valid_we <- function(.dataVec) {
    if(nrow(filter(Ldfve, Toolkits == .dataVec$TK & Variant == .dataVec$name & Valid == TRUE)) > 0){
        pltResEnergyLandscape <- ggplot(filter(Ldfve, Toolkits == .dataVec$TK & Variant == .dataVec$name & Valid == TRUE & QiskitOptLevel == 3), aes(x=Cost, colour = as.factor(p))) +
        geom_freqpoly(position = pos, bins = bins) + 
        facet_grid(cols=vars(Penalty1), rows=vars(Penalty2)) +
        geom_vline(aes(xintercept=OptimalCost)) + 
        scale_colour_manual("# LR-QAOA Layers", values = LFD.COLOURS) +
        ylab("Count") +
        theme_paper_base()

        ggsave(plot = pltResEnergyLandscape, filename = paste("LRQAOAResults", .dataVec$name, "ValidErrorMitigatedCostLandscape", .dataVec$TK, "TK.pdf", sep=""), path=paste(OUT.PATH, .dataVec$TK, "_Toolkits/", sep=""), device=cairo_pdf)
    }
}

for (vec in plotVecs){
    L_multi_plot_energy(vec)
    L_multi_plot_valid(vec)
    L_multi_plot_valid_we(vec)
}

relable_p <- function(orig) {
    return(paste("p=", orig, sep=""))
}

relable_OptLevel <- function(orig) {
    return(paste("Opt. level ", orig, sep=""))
}

LdfLCTC <- Ldf %>% 
    rowwise() %>%
    mutate(QiskitOptLevel = relable_OptLevel(QiskitOptLevel), p = relable_p(p))

LdfLCTC$p_f = factor(LdfLCTC$p, levels=c("p=1", "p=2", "p=5", "p=10")) 

pltResLCTCInstr <-ggplot(LdfLCTC, aes(x=LCNumInstructions, y=TCNumInstructions, colour=as.factor(Toolkits))) +
    geom_point(size=LINE.WIDTH, alpha=POINT.ALPHA) +
    facet_grid(cols=vars(p_f), rows=vars(QiskitOptLevel)) +
    scale_colour_manual("# Toolkits", values = LFD.COLOURS) +
    xlab("# Instructions in logical circuit") +
    ylab("# Instructions in transpiled circuit") +
    theme_paper_base() 

ggsave(plot = pltResLCTCInstr, filename = "LRQAOAResultsLCvsTCInstructions.pdf", path = OUT.PATH, device=cairo_pdf)

pltResLCTCDepth <-ggplot(LdfLCTC, aes(x=LCDepth, y=TCDepth, colour=as.factor(Toolkits))) +
    geom_point(size=LINE.WIDTH, alpha=POINT.ALPHA) +
    facet_grid(cols=vars(p_f), rows=vars(QiskitOptLevel)) +
    scale_colour_manual("# Toolkits", values = LFD.COLOURS) +
    xlab("Logical circuit depth") +
    ylab("Transpiled circuit depth") +
    theme_paper_base() 

ggsave(plot = pltResLCTCDepth, filename = "LRQAOAResultsLCvsTCDepth.pdf", path = OUT.PATH, device=cairo_pdf)

plotVecs <- list()
for (TK in c(3,9,13,16,18,19)){
    plotVecs <- list.append(plotVecs, list(TK=TK, name="raw"    ))
    plotVecs <- list.append(plotVecs, list(TK=TK, name="rounded"))
    plotVecs <- list.append(plotVecs, list(TK=TK, name="scaled" ))
}

Adf <- Adf   %>% mutate(Penalty1 = relableP1(Penalty1), Penalty2 = relableP2(Penalty2))
Adfv <- Adfv  %>% mutate(Penalty1 = relableP1(Penalty1), Penalty2 = relableP2(Penalty2))
Adfve <- Adfve %>% mutate(Penalty1 = relableP1(Penalty1), Penalty2 = relableP2(Penalty2))

#.dataVec: (Toolkits, name)
A_multi_plot_energy <- function(.dataVec) {
    pltResEnergyLandscape <- ggplot(filter(Adf, Toolkits == .dataVec$TK & Variant == .dataVec$name), aes(x=energy, colour = as.factor(Schedule))) +
    geom_freqpoly(position = "jitter", bins = bins) + 
    facet_grid(cols=vars(Penalty1), rows=vars(Penalty2)) +
    scale_colour_manual("Schedule type", values = LFD.COLOURS) +
    ylab("Count") +
    theme_paper_base() + 
    guide_paper_base()

    ggsave(plot = pltResEnergyLandscape, filename = paste("QuantumAnnealingResults", .dataVec$name, "EnergyLandscape", .dataVec$TK, "TK.pdf", sep=""), path=paste(OUT.PATH, .dataVec$TK, "_Toolkits/", sep=""), device=cairo_pdf)
}

#.dataVec: (Toolkits, name)
A_multi_plot_valid <- function(.dataVec){
    pltValidCostResultCarlos <- ggplot(filter(Adfv, Toolkits == .dataVec$TK & Variant == .dataVec$name & Valid == TRUE), aes(x=Cost, colour=as.factor(Annealing_time))) +
        geom_histogram(bins=bins) +
        facet_grid(cols=vars(Penalty1), rows=vars(Penalty2)) +
        geom_vline(aes(xintercept=OptimalCost)) + 
        scale_colour_manual("Annealing time", values = LFD.COLOURS) +
        scale_fill_gradient2("Valid count", low="#fcfcfc", high="#aba9a9") +
        ylab("Valid count") +
        theme_paper_base() +
        guide_paper_base()
    ggsave(plot=pltValidCostResultCarlos, filename = paste("QuantumAnnealingResults", .dataVec$name, "ValidCostLandscape", .dataVec$TK, "TK.pdf", sep=""), path=paste(OUT.PATH, .dataVec$TK, "_Toolkits/", sep=""), device=cairo_pdf)
}

#.dataVec: (Toolkits, name)
A_multi_plot_valid_we <- function(.dataVec){
    pltValidCostResultCarlos <- ggplot(filter(Adfve, Toolkits == .dataVec$TK & Variant == .dataVec$name & Valid == TRUE), aes(x=Cost, colour=as.factor(Annealing_time))) +
        geom_histogram(bins=bins) +
        facet_grid(cols=vars(Penalty1), rows=vars(Penalty2)) +
        geom_vline(aes(xintercept=OptimalCost)) + 
        scale_colour_manual("Annealing time", values = LFD.COLOURS) +
        scale_fill_gradient2("Valid, error mitigated, count", low="#fcfcfc", high="#aba9a9") +
        ylab("Valid count") +
        theme_paper_base() +
        guide_paper_base()
    ggsave(plot=pltValidCostResultCarlos, filename = paste("QuantumAnnealingResults", .dataVec$name, "ValidErrorMitigatedCostLandscape", .dataVec$TK, "TK.pdf", sep=""), path=paste(OUT.PATH, .dataVec$TK, "_Toolkits/", sep=""), device=cairo_pdf)
}

for (vec in plotVecs){
    A_multi_plot_energy(vec)
    A_multi_plot_valid(vec)
    A_multi_plot_valid_we(vec)
}


#############################################################################
# Simulated Annealing plots
#############################################################################

S_multi_plot_sa_energy <- function(.dataVec){
    pltResEnergyLandscape <- ggplot(filter(Sdf, Toolkits == .dataVec$TK & Variant == .dataVec$name), aes(x=energy)) +
        geom_freqpoly(position = "jitter", bins = bins) + 
        facet_grid(cols=vars(Penalty1), rows=vars(Penalty2)) +
        ylab("Count") +
        theme_paper_base() + 
        guide_paper_base()

    ggsave(plot = pltResEnergyLandscape, filename = paste("SimulatedAnnealingResults", .dataVec$name, "EnergyLandscape", .dataVec$TK, "TK.pdf", sep=""), path=paste(OUT.PATH, .dataVec$TK, "_Toolkits/", sep=""), device=cairo_pdf)
}

S_multi_plot_sa_valid <- function(.dataVec){
    pltValidCostResultCarlos <- ggplot(filter(Sdfv, Toolkits == .dataVec$TK & Variant == .dataVec$name & Valid == TRUE), aes(x=Cost)) +
        geom_histogram(bins=bins) +
        facet_grid(cols=vars(Penalty1), rows=vars(Penalty2)) +
        geom_vline(aes(xintercept=OptimalCost)) + 
        scale_fill_gradient2("Valid count", low="#fcfcfc", high="#aba9a9") +
        ylab("Valid count") +
        theme_paper_base() +
        guide_paper_base()
    ggsave(plot=pltValidCostResultCarlos, filename = paste("SimulatedAnnealingResults", .dataVec$name, "ValidCostLandscape", .dataVec$TK, "TK.pdf", sep=""), path=paste(OUT.PATH, .dataVec$TK, "_Toolkits/", sep=""), device=cairo_pdf)
}

S_multi_plot_sa_valid_we <- function(.dataVec){
    pltValidCostResultCarlos <- ggplot(filter(Sdfve, Toolkits == .dataVec$TK & Variant == .dataVec$name & Valid == TRUE), aes(x=Cost)) +
        geom_histogram(bins=bins) +
        facet_grid(cols=vars(Penalty1), rows=vars(Penalty2)) +
        geom_vline(aes(xintercept=OptimalCost)) + 
        scale_fill_gradient2("Valid, error mitigated, count", low="#fcfcfc", high="#aba9a9") +
        ylab("Valid count") +
        theme_paper_base() +
        guide_paper_base()
    ggsave(plot=pltValidCostResultCarlos, filename = paste("SimulatedAnnealingResults", .dataVec$name, "ValidErrorMitigatedCostLandscape", .dataVec$TK, "TK.pdf", sep=""), path=paste(OUT.PATH, .dataVec$TK, "_Toolkits/", sep=""), device=cairo_pdf)
}

for (vec in plotVecs){
    S_multi_plot_sa_energy(vec)
    S_multi_plot_sa_valid(vec)
    S_multi_plot_sa_valid_we(vec)
}