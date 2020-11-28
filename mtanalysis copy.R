#####Mt DNA Probability estimation
library(tidyverse)
library(data.table)
library(variantprobs)
library(magrittr)
library(reshape)
library(reshape2)
source("calc_minfo.R")

mtvars<-data.table(read.table("TCMA-MutationSNV.tsv", header=TRUE))
mtvars$Start<-mtvars$position
mtvars$Stop<-mtvars$position

mtvars<-setDT(mtvars)
anno<-data.table(read.csv("mtanno.csv", header=TRUE))
anno<-setDT(anno[,1:3])
setkey(anno, Start, Stop)
annotated_mt<-foverlaps(mtvars, anno, by.y=key(anno), type="any", nomatch=0L)

dt<-annotated_mt %>% dplyr::select(Gene, sample_id, cancer_type, i.Start, ref, var, var_type)

mt_v_f <- dt %>% group_by(cancer_type) %>% mutate(n_tumor=length(unique(sample_id)))

mt_v_f<-mt_v_f %>% group_by(i.Start, ref, var, Gene, cancer_type) %>% mutate(v_f = length(unique(sample_id)), n_tumor = n_tumor[1])

mt_v_f<-mt_v_f %>% ungroup() %>% setDT()

mtga_newprobs_given_cancer<-mt_v_f[,
       # Calculate Good Turing probabilities of
       # at least one new variants per gene & cancer type
       {
         GT_probs <- goodturing_probs(
           counts = v_f,
           m = n_tumor[1]
         )
         .(p_atleast_1new_v = GT_probs['atleast_1new'])
       },
       by = .(Gene, cancer_type)
       ] %>%
  dcast.data.table(
    Gene ~ cancer_type,
    value.var = "p_atleast_1new_v",
    fill = 1 - exp(- 1/(length(unique(mt_v_f$sample_id)) + 1))
  ) %>%
  magrittr::set_rownames(.$Gene) %>%
  .[, Gene := NULL] %>%
  data.matrix()

mtga_newprobs_given_cancer %>% melt() %>% ggplot()+geom_point(aes(x=Var2, y=Var1, size=value, color=Var1))+theme_bw()+theme(axis.text.x=element_text(angle=90), legend.position="none")

cancer_npatient <- mt_v_f %>% group_by(cancer_type) %>% summarize(cancer_npatient = length(unique(sample_id))) %>% ungroup()
probs<-cancer_npatient$cancer_npatient/sum(cancer_npatient$cancer_npatient)
cancer_prob <- probs %>%
  setNames(cancer_npatient$cancer_type)

nmi <- calc_minfo(
  mtga_newprobs_given_cancer,
  cancer_prob,
  binary_minfo = FALSE,
  normalize = TRUE
)

nmi

mtvars<-data.table(read.table("TCMA-MutationINDEL.tsv", header=TRUE))
mtvars$Start<-mtvars$position
mtvars$Stop<-mtvars$position

mtvars<-setDT(mtvars)
anno<-data.table(read.csv("mtanno.csv", header=TRUE))
anno<-setDT(anno[,1:3])
setkey(anno, Start, Stop)
annotated_mt<-foverlaps(mtvars, anno, by.y=key(anno), type="any", nomatch=0L)

dt<-annotated_mt %>% dplyr::select(Gene, sample_id, cancer_type, i.Start, ref, var)

mt_v_f <- dt %>% group_by(cancer_type) %>% mutate(n_tumor=length(unique(sample_id)))

mt_v_f<-mt_v_f %>% group_by(i.Start, ref, var, Gene, cancer_type) %>% mutate(v_f = length(unique(sample_id)), n_tumor = n_tumor[1])

mt_v_f<-mt_v_f %>% ungroup() %>% setDT()

mtga_newprobs_given_cancer<-mt_v_f[,
                                   # Calculate Good Turing probabilities of
                                   # at least one new variants per gene & cancer type
                                   {
                                     GT_probs <- goodturing_probs(
                                       counts = v_f,
                                       m = n_tumor[1]
                                     )
                                     .(p_atleast_1new_v = GT_probs['atleast_1new'])
                                   },
                                   by = .(Gene, cancer_type)
                                   ] %>%
  dcast.data.table(
    Gene ~ cancer_type,
    value.var = "p_atleast_1new_v",
    fill = 1 - exp(- 1/(length(unique(mt_v_f$sample_id)) + 1))
  ) %>%
  magrittr::set_rownames(.$Gene) %>%
  .[, Gene := NULL] %>%
  data.matrix()

mtga_newprobs_given_cancer %>% melt() %>% ggplot()+geom_point(aes(x=Var2, y=Var1, size=value, color=Var1))+theme_bw()+theme(axis.text.x=element_text(angle=90), legend.position="none")

cancer_npatient <- mt_v_f %>% group_by(cancer_type) %>% summarize(cancer_npatient = length(unique(sample_id))) %>% ungroup()
probs<-cancer_npatient$cancer_npatient/sum(cancer_npatient$cancer_npatient)
cancer_prob <- probs %>%
  setNames(cancer_npatient$cancer_type)

nmi <- calc_minfo(
  mtga_newprobs_given_cancer,
  cancer_prob,
  binary_minfo = FALSE,
  normalize = TRUE
)

nmi

mtvars<-data.table(read.table("TCMA-NuclearTransfer.tsv", header=TRUE))
mtvars<-mtvars[,1:5]
colnames(mtvars)<-c("sample_id", "cancer_type", "chr", "Start", "Stop")

mtvars<-setDT(mtvars)
anno<-data.table(read.csv("mtanno.csv", header=TRUE))
anno<-setDT(anno[,1:3])
setkey(anno, Start, Stop)
annotated_mt<-foverlaps(mtvars, anno, by.y=key(anno), type="any", nomatch=0L)

dt<-annotated_mt %>% dplyr::select(Gene, sample_id, cancer_type, i.Start)

mt_v_f <- dt %>% group_by(cancer_type) %>% mutate(n_tumor=length(unique(sample_id)))

mt_v_f<-mt_v_f %>% group_by(i.Start, Gene, cancer_type) %>% mutate(v_f = length(unique(sample_id)), n_tumor = n_tumor[1])

mt_v_f<-mt_v_f %>% ungroup() %>% setDT()

mtga_newprobs_given_cancer<-mt_v_f[,
                                   # Calculate Good Turing probabilities of
                                   # at least one new variants per gene & cancer type
                                   {
                                     GT_probs <- goodturing_probs(
                                       counts = v_f,
                                       m = n_tumor[1]
                                     )
                                     .(p_atleast_1new_v = GT_probs['atleast_1new'])
                                   },
                                   by = .(Gene, cancer_type)
                                   ] %>%
  dcast.data.table(
    Gene ~ cancer_type,
    value.var = "p_atleast_1new_v",
    fill = 1 - exp(- 1/(length(unique(mt_v_f$sample_id)) + 1))
  ) %>%
  magrittr::set_rownames(.$Gene) %>%
  .[, Gene := NULL] %>%
  data.matrix()

mtga_newprobs_given_cancer %>% melt() %>% ggplot()+geom_point(aes(x=Var2, y=Var1, size=value, color=Var1))+theme_bw()+theme(axis.text.x=element_text(angle=90), legend.position="none")

cancer_npatient <- mt_v_f %>% group_by(cancer_type) %>% summarize(cancer_npatient = length(unique(sample_id))) %>% ungroup()
probs<-cancer_npatient$cancer_npatient/sum(cancer_npatient$cancer_npatient)
cancer_prob <- probs %>%
  setNames(cancer_npatient$cancer_type)

nmi <- calc_minfo(
  mtga_newprobs_given_cancer,
  cancer_prob,
  binary_minfo = FALSE,
  normalize = TRUE
)

nmi

