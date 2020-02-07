library(igraph)
## remove edges between M1+H and M2+Na, M3+H and M4+K, in general when there is 
## a link between an adduct and a molecular ion

## node with outgoing transf link should not have outgoing adduct link when
## the adduct-linking vertex has a ingoing adduct link

## node with ingoing transf link should not have ingoing adduct link when the 
## adduct-linking vertex has an outgoing adduct link

## pos mode: isotopic+1, isotopic+2, isotopic+3, Na adduct, K adduct, 
## formic acid adduct
## neg mode: isotopic+1, isotopic+2, isotopic+3, formic acid, Na adduct,
## Sodium formate 


##          M1+H    M1+Na   M2+H    M2+Na   M3+H    M3+Na
## M1+H     0       1       0       0       0       1       
## M1+Na    0       0       0       0       0       0
## M2+H     1       0       0       1       0       0
## M2+Na    0       1       0       0       1       0
## M3+H     0       1       0       0       0       1
## M3+Na    0       0       0       0       0       0  

# mat_test <- matrix(c(c(0, 1, 1, 0, 0, 1), ## M1+H
#                      c(1, 0, 1, 1, 1, 1), ## M1+Na
#                      c(1, 1, 0, 1, 0, 0), ## M2+H
#                      c(0, 1, 1, 0, 1, 0), ## M2+Na
#                      c(0, 1, 0, 1, 0, 1), ## M3+H
#                      c(0, 1, 0, 0, 1, 0)), ## M3+Na
#                    

## example 1
mat_1 <- matrix(c(c(0, 1, 0, 0, 0, 1), ## M1+H
                     c(0, 0, 0, 0, 0, 0), ## M1+Na
                     c(1, 0, 0, 1, 0, 0), ## M2+H
                     c(0, 1, 0, 0, 0, 0), ## M2+Na
                     c(0, 1, 0, 0, 0, 1), ## M3+H
                     c(0, 0, 0, 0, 0, 0)), ## M3+Na
                   byrow = TRUE, ncol = 6, nrow = 6)
mat_1_transf <- matrix(c(
    c("", "adduct_Na", "", "", "", "transf_gluc_F"), ## M1+H
    c("", "", "", "", "", ""), ## M1+Na
    c("transf_gluc_T", "", "", "adduct_Na", "", ""), ## M2+H
    c("", "transf_gluc_T", "", "", "", ""), ## M2+Na
    c("", "transf_gluc_F", "", "", "", "adduct_Na"), ## M3+H
    c("", "", "", "", "", "")), ## M3+Na
    byrow = TRUE, ncol = 6, nrow = 6)
colnames(mat_1) <- c("M1+H", "M1+Na", "M2+H", "M2+Na", "M3+H", "M3+Na")
rownames(mat_1) <- colnames(mat_1)
net_mat_1 <- graph_from_adjacency_matrix(mat_1, mode = "directed")
plot(net_mat_1)

## example 2
mat_2 <- matrix(c(c(0, 1, 1, 0, 1, 0, 0, 0, 0), ## "M_1"
                      c(0, 0, 1, 1, 0, 0, 0, 1, 0), ## "M+1_1"
                      c(0, 0, 0, 0, 0, 0, 0, 0, 0), ## "M+2_1"
                      c(0, 0, 0, 0, 0, 0, 0, 0, 0), ## "M+1+45_1"
                      c(0, 0, 0, 0, 0, 0, 0, 0, 0), ## "M+152_2"
                      c(0, 1, 0, 0, 0, 0, 0, 0, 0), ## "M+1-152_3"
                      c(0, 1, 0, 0, 0, 0, 0, 0, 0), ## "M+1-152_4"
                      c(0, 0, 0, 0, 0, 0, 0, 0, 0), ##"M+1+152_5"
                      c(0, 0, 1, 0, 0, 0, 0, 0, 0)), ## "M+2-152_6"
                    byrow = TRUE, ncol = 9, nrow = 9)
mat_2_transf <- matrix(c(
    c("", "adduct_+1",  "adduct_+2",  "",          "transf_152", "", "", "",           ""), ## "M_1"
    c("", "",           "adduct_+1",  "transf_45", "",           "", "", "transf_152", ""), ## "M+1_1"
    c("", "",           "",           "",          "",           "", "", "",           ""), ## "M+2_1"
    c("", "",           "",           "",          "",           "", "", "",           ""), ## "M+1+45_1"
    c("", "",           "",           "",          "",           "", "", "",           ""), ## "M+152_2"
    c("", "transf_152", "",           "",          "",           "", "", "",           ""), ## "M+1-152_3"
    c("", "transf_152", "",           "",          "",           "", "", "",           ""), ## "M+1-152_4"
    c("", "",           "",           "",          "",           "", "", "",           ""), ##"M+1+152_5"
    c("", "",           "transf_152", "",          "",           "", "", "",           "")), ## "M+2-152_6"
    byrow = TRUE, ncol = 9, nrow = 9)
colnames(mat_2) <- c("M_1", "M+1_1", "M+2_1", "M+1+45_1", "M+152_2", "M+1-152_3", "M+1-152_4", "M+1+152_5", "M+2-152_6")
rownames(mat_2) <- colnames(mat_2)
colnames(mat_2_transf) <- rownames(mat_2_transf) <- colnames(mat_2)
net_mat_2 <- graph_from_adjacency_matrix(mat_2, mode = "directed")
plot(net_mat_2)

## example 3 
mat_3 <- matrix(c(c(0, 1, 1, 0, 1, 0, 0, 0, 0), ## "M_1"
                      c(0, 0, 1, 1, 0, 0, 0, 1, 0), ## "M+1_1"
                      c(0, 0, 0, 0, 0, 0, 0, 0, 1), ## "M+2_1"
                      c(0, 0, 0, 0, 0, 0, 0, 0, 0), ## "M+1+45_1"
                      c(0, 0, 0, 0, 0, 0, 0, 1, 1), ## "M+152_2"
                      c(0, 1, 0, 0, 0, 0, 0, 0, 0), ## "M+1-152_3"
                      c(0, 1, 0, 0, 0, 0, 0, 0, 0), ## "M+1-152_4"
                      c(0, 0, 0, 0, 0, 0, 0, 0, 1), ##"M+1+152_2"
                      c(0, 0, 0, 0, 0, 0, 0, 0, 0)), ## "M+2_2"
                    byrow = TRUE, ncol = 9, nrow = 9)
mat_3_transf <- matrix(c(
    c("", "adduct_+1",  "adduct_+2",  "",          "transf_152", "", "", "",           ""), ## "M_1"
    c("", "",           "adduct_+1",  "transf_45", "",           "", "", "transf_152", ""), ## "M+1_1"
    c("", "",           "",           "",          "",           "", "", "",           "transf_152"), ## "M+2_1"
    c("", "",           "",           "",          "",           "", "", "",           ""), ## "M+1+45_1"
    c("", "",           "",           "",          "",           "", "", "adduct_+1", "adduct_+2"), ## "M+152_2"
    c("", "transf_152", "",           "",          "",           "", "", "",           ""), ## "M+1-152_3"
    c("", "transf_152", "",           "",          "",           "", "", "",           ""), ## "M+1-152_4"
    c("", "",           "",           "",          "",           "", "", "",           "adduct_+1"), ##"M+1+152_2"
    c("", "",           "", "",          "",           "", "", "",           "")), ## "M+2+152_2"
    byrow = TRUE, ncol = 9, nrow = 9)
colnames(mat_3) <- c("M_1", "M+1_1", "M+2_1", "M+1+45_1", "M+152_2", "M+1-152_3", "M+1-152_4", "M+1+152_2", "M+2+152_2")
rownames(mat_3) <- colnames(mat_3)
colnames(mat_3_transf) <- rownames(mat_3_transf) <- colnames(mat_3)
net_mat_3 <- graph_from_adjacency_matrix(mat_3, mode = "directed")
plot(net_mat_3)

## example 4
mat_4 <- matrix(c(c(0, 1, 1, 1, 1, 0, 0, 0, 0), ## "M_1"
                  c(0, 0, 1, 0, 0, 0, 0, 1, 0), ## "M+1_1"
                  c(0, 0, 0, 0, 0, 0, 0, 0, 1), ## "M+2_1"
                  c(0, 0, 0, 0, 0, 0, 0, 0, 0), ## "M+45_5"
                  c(0, 0, 0, 0, 0, 0, 0, 1, 1), ## "M+152_2"
                  c(1, 0, 0, 0, 0, 0, 0, 0, 0), ## "M+1-152_3"
                  c(0, 1, 0, 0, 0, 0, 0, 0, 0), ## "M+1-152_4"
                  c(0, 0, 0, 0, 0, 0, 0, 0, 1), ##"M+1+152_2"
                  c(0, 0, 0, 0, 0, 0, 0, 0, 0)), ## "M+2_2"
                byrow = TRUE, ncol = 9, nrow = 9)
mat_4_transf <- matrix(c(
    c("", "adduct_+1",  "adduct_+2",  "transf_45",          "transf_152", "", "", "",           ""), ## "M_1"
    c("", "",           "adduct_+1",  "", "",           "", "", "transf_152", ""), ## "M+1_1"
    c("", "",           "",           "",          "",           "", "", "",           "transf_152"), ## "M+2_1"
    c("", "",           "",           "",          "",           "", "", "",           ""), ## "M+45_5"
    c("", "",           "",           "",          "",           "", "", "adduct_+1", "adduct_+2"), ## "M+152_2"
    c("transf_152", "", "",           "",          "",           "", "", "",           ""), ## "M-152_3"
    c("", "transf_152", "",           "",          "",           "", "", "",           ""), ## "M+1-152_4"
    c("", "",           "",           "",          "",           "", "", "",           "adduct_+1"), ##"M+1+152_2"
    c("", "",           "", "",          "",           "", "", "",           "")), ## "M+2+152_2"
    byrow = TRUE, ncol = 9, nrow = 9)
colnames(mat_4) <- c("M_1", "M+1_1", "M+2_1", "M+45_5", "M+152_2", "M+1-152_3", "M+1-152_4", "M+1+152_2", "M+2+152_2")
rownames(mat_4) <- colnames(mat_4)
colnames(mat_4_transf) <- rownames(mat_4_transf) <- colnames(mat_4)
net_mat_4 <- graph_from_adjacency_matrix(mat_4, mode = "directed")
plot(net_mat_4)

## example 5
mat_5 <- matrix(c(
    c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0), ## "M_1"
    c(0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0), ## "M+1_1"
    c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0), ## "M+2_1"
    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0), ## "M+45_5"
    c(0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0), ## "M+152_2"
    c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), ## "M+1-152_3"
    c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), ## "M+1-152_4"
    c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0), ##"M+1+152_2"
    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), ## "M+2_2"
    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), ## "M+2+45_6"
    c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)), ## "M-152-45_7"            
                byrow = TRUE, ncol = 11, nrow = 11)
mat_5_transf <- matrix(c(
    c("", "adduct_+1",  "adduct_+2",  "transf_45",          "transf_152", "", "", "",           "", "", ""), ## "M_1"
    c("", "",           "adduct_+1",  "", "",           "", "", "transf_152", "", "", ""), ## "M+1_1"
    c("", "",           "",           "",          "",           "", "", "",           "transf_152", "transf_45", ""), ## "M+2_1"
    c("", "",           "",           "",          "",           "", "", "",           "", "adduct_+2", ""), ## "M+45_5"
    c("", "",           "",           "",          "",           "", "", "adduct_+1", "adduct_+2", "", ""), ## "M+152_2"
    c("transf_152", "", "",           "",          "",           "", "", "",           "", "", ""), ## "M-152_3"
    c("", "transf_152", "",           "",          "",           "", "", "",           "", "", ""), ## "M+1-152_4"
    c("", "",           "",           "",          "",           "", "", "",           "adduct_+1", "", ""), ##"M+1+152_2"
    c("", "",           "", "",          "",           "", "", "",           "", "", ""),## "M+2+152_2"
    c("", "",           "", "",          "",           "", "", "",           "", "", ""), ## "M+2+45_6"
    c("", "",           "", "",          "",           "transf_45", "", "",           "", "", "")), ## "M-152-45_7"
    byrow = TRUE, ncol = 11, nrow = 11)
colnames(mat_5) <- c("M_1", "M+1_1", "M+2_1", "M+1+45_5", "M+152_2", "M-152_3", "M+1-152_4", "M+1+152_2", "M+2+152_2", "M+2+45_6", "M-152-45_7")
rownames(mat_5) <- colnames(mat_5)
colnames(mat_5_transf) <- rownames(mat_5_transf) <- colnames(mat_5)
net_mat_5 <- graph_from_adjacency_matrix(mat_5, mode = "directed")
plot(net_mat_5)

## node with outgoing transf link should not have outgoing adduct link when
## the adduct-linking vertex has an ingoing adduct link of same type

## do for all components
## iterate through rows

removeFalseLinksAdducts <- function(mat_l) {
    
    mat_num <- mat_l[[1]]
    mat_char <- mat_l[[2]]
    
    for (i in seq_len(nrow(mat_num))) {
        
        rows_i <- mat_char[i, ]
        adduct_i <- rows_i[grep(rows_i, pattern = "adduct_")]
        
        ## if any vertex has outgoing adduct link, check all cols where rows_i
        ## has transf if ingoing adduct link exists for adduct-linking feature
        if (length(adduct_i) > 0) { 
            
            inds_row <- grep(rows_i, pattern = "transf_")
            
            for (j in inds_row) {
                cols_j <- mat_char[, j]
                adduct_j <- cols_j[grep(cols_j, pattern = "adduct")]
                if (length(adduct_j) > 0)
                    if (any(adduct_j %in% adduct_i)) {
                        mat_num[i, j] <- 0 ## remove in mat_test[i, ]
                    }
            }
        }
        
        ## new 
        cols_i <- mat_char[, i]
        adduct_i <- cols_i[grep(cols_i, pattern = "adduct_")]
        
        ## if any vertex has ingoing adduct link, check all rows where cols_i
        ## has transf if outgoing adduct link exists for adduct-linking feature
        if (length(adduct_i) > 0) { 
            
            inds_col <- grep(cols_i, pattern = "transf_")
            
            for (j in inds_col) {
                rows_j <- mat_char[j, ]
                adduct_j <- rows_j[grep(rows_j, pattern = "adduct")]
                if (length(adduct_j) > 0)
                    if (any(adduct_j %in% adduct_i)) {
                        mat_num[j, i] <- 0 ## remove in mat_test[i, ]
                    }
            }
        }
        
    }
    
    ## set mat_char to "" where mat_num == 0
    mat_char[which(mat_num == 0)] <- ""
    
    return(list(mat_num, mat_char))
}


#' @name removeFalseLinksCircular
#' @title Delete in- and outgoing links based on their relation to the molecular
#' ion
#' @author Thomas Naake <thomasnaake@@googlemail.com>
#' @description The function `removeFalseLinksCircular` deletes ingoing 
#' and outgoing links of isotopes/adducts to M if there is no circular relation from 
## M to M+isot/adduct to M+transf to M+iso/adduct+transformation.
#' Example:
#' 1) M and M+1 forming an "adduct_isotope+1", suppose there is a "transf_+162"
#' from M+1 to another mass feature M2, but there is no link between M and M2-1,
#' then remove the link M+1 and M2.
#' 2) M and M+1 forming an "adduct_isotope+1", suppose there is a "transf_+162"
#' from M+1 to another mass feature M, AND there is a link between M and M2-1,
#' then keep the the link M+1 and M2 (support that mass feature is a metabolite)
#' @details 
#' @example 
#' 
removeFalseLinksCircular <- function(mat_l) {
    
    mat_num <- mat_l[[1]]
    mat_char <- mat_l[[2]]
    
    ## find relations with adducts: 
    ## Molecular ion (MI): suppose this is the vertex that has no ingoing 
    ## adduct_ link, i.e. the column of M does not have a "adduct_"
    M <- apply(mat_char, 2, grep, pattern = "adduct_")
    
    ## iterate through M where M is not integer(0). These are adducts of
    ## the rows i 
    for (i in seq_len(length(M))) {
        
        if (length(M[[i]])) {
            ## if length(M[[i]]) != 0, M[[i]] is not the MI:
            ## check if the MI has transf_
            transf_M_og_ind <- grep(mat_char[i, ], pattern = "transf_")
            transf_M_ig_ind <- grep(mat_char[, i], pattern = "transf_")
            
            if (length(transf_M_og_ind)) {
                
                ## check if the MI has outgoing transf_ links and return the 
                ## indices
                tmp <- matrix(mat_char[M[[i]], ], byrow = FALSE,
                    ncol = ncol(mat_char))
                transf_MI_og_ind <- apply(tmp, 2, grepl, pattern = "transf_")
                transf_MI_og_ind <- which(transf_MI_og_ind, arr.ind = TRUE)
                
                if (length(transf_MI_og_ind)) {
                    ## iterate thorough all transf_M_og_ind (outgoing links 
                    ## from M) 
                    for(j in 1:length(transf_M_og_ind)) {
                        ind_j <- transf_M_og_ind[j]
                        ## if there is only one MI do
                        if (!is.matrix(transf_MI_og_ind)) {
                            remove <- TRUE
                            for (k in 1:length(transf_MI_og_ind)) {
                                ## get adduct of M
                                adduct_ij <- mat_char[transf_MI_og_ind[k], ind_j]
                                if (adduct_ij == mat_char[M[[i]], i]) {
                                    remove <- c(remove, FALSE)
                                }
                            }
                            ## remove link transf_M_og_ind[k] at row i if there
                            ## is one matching adduct
                            if (all(remove)) {mat_num[i, ind_j] <- 0}
                        ## if there are several MI do
                        } else {
                            remove <- TRUE
                            for(k in 1:nrow(transf_MI_og_ind)) {
                                ## get adduct of M
                                adduct_ij <- mat_char[transf_MI_og_ind[k, "row"], ind_j]
                                if (adduct_ij == mat_char[i, M[[i]][transf_MI_og_ind[k, "row"]]]) {
                                    remove <- c(remove, FALSE)   
                                }
                            }
                            ## remove link transf_M_og_ind[k] at row i if there
                            ## is one matching adduct
                            if (all(remove)) {mat_num[i, ind_j] <- 0}
                        }
                    }
                } else {
                    ## remove all links by default when there is no outgoing
                    ## transf_ from the MI
                    mat_num[i, transf_M_ig_ind] <- 0
                }
            }
            ## check only ingoing links now
            if (length(transf_M_ig_ind)) {

                ## check if the MI has ingoing transf_ links and return the
                ## indices
                tmp <- matrix(mat_char[, M[[i]]], byrow = TRUE,
                ncol = ncol(mat_char))
                transf_MI_ig_ind <- apply(tmp, 2, grepl, pattern = "transf_")
                transf_MI_ig_ind <- which(transf_MI_ig_ind, arr.ind = TRUE)

                if (length(transf_MI_ig_ind)) {
                    ## iterate thorough all transf_M_ig_ind (ingoing links from
                    ## M)
                    for(j in 1:length(transf_M_ig_ind)) {
                    
                        ind_j <- transf_M_ig_ind[j]
                        ## if there is only one MI do
                        if (!is.matrix(transf_MI_ig_ind)) {
                            remove <- TRUE
                            for (k in 1:length(transf_MI_ig_ind)) {
                                ## get adduct of M
                                adduct_ij <- mat_char[transf_MI_ig_ind[k], ind_j]
                                
                                if (adduct_ij == mat_char[M[[i]], i]) {
                                    remove <- c(remove, FALSE)
                                }
                            }
                            ## remove link transf_M_ig_ind[j] at row ind_j if 
                            ## there is one matching adduct
                            if (all(remove)) {mat_num[ind_j, i] <- 0}
                        ## if there are several MI do
                        } else {
                            remove <- TRUE
                            for(k in 1:nrow(transf_MI_ig_ind)) {
                                ## get adduct of M
                                adduct_ij <- mat_char[transf_MI_ig_ind[k, "row"], ind_j]
                                
                                if (adduct_ij == mat_char[M[[i]][transf_MI_ig_ind[k, "row"]], i]) {
                                    remove <- c(remove, FALSE)
                                }
                            }
                            ## remove link transf_M_ig_ind[k] at row ind_j if 
                            ## there is one matching adduct
                            if (all(remove)) {mat_num[ind_j, i] <- 0}
                        }
                    }
                } else {
                    ## remove all links by default when there is no ingoing
                    ## transf_ from the MI
                    mat_num[transf_M_ig_ind, i] <- 0
                }
            }
        }
    }
    
    ## set mat_char to "" where mat_num == 0 (removes links in mat_char)
    mat_char[which(mat_num == 0)] <- ""
    
    return(list(mat_num, mat_char))
    
}

## example 1
mat_l <- list(mat_1, mat_1_transf)
mat_rem <- removeFalseLinksAdducts(list(mat_1, mat_1_transf))
mat_rem <- removeFalseLinksCircular(mat_l)
plot(graph_from_adjacency_matrix(mat_l[[1]], mode = "directed"))
plot(graph_from_adjacency_matrix(mat_rem[[1]], mode = "directed"))

# example 2 
mat_l <- list(mat_2, mat_2_transf)
mat_rem <- removeFalseLinksAdducts(mat_l)
mat_rem <- removeFalseLinksCircular(mat_rem)
plot(graph_from_adjacency_matrix(mat_l[[1]], mode = "directed"))
plot(graph_from_adjacency_matrix(mat_rem[[1]], mode = "directed")) ## should only result in link between M_1->M+152_2, M_1-> M_2, M_1_M_3, M_2->M_3

## example 3
mat_rem <- removeFalseLinksAdducts(list(mat_3, mat_3_transf))
mat_l <- list(mat_3, mat_3_transf)
mat_rem <- removeFalseLinksCircular(list(mat_3, mat_3_transf))
plot(graph_from_adjacency_matrix(mat_3), mode = "directed")
plot(graph_from_adjacency_matrix(mat_rem[[1]], mode = "directed"))

## example 4
mat_l <- list(mat_4, mat_4_transf)
mat_rem <- removeFalseLinksAdducts(mat_l)
mat_rem <- removeFalseLinksCircular(list(mat_4, mat_4_transf))
plot(graph_from_adjacency_matrix(mat_4), mode = "directed")
plot(graph_from_adjacency_matrix(mat_rem[[1]], mode = "directed"))

## example 5
mat_rem <- removeFalseLinksAdducts(list(mat_5, mat_5_transf))
mat_l <- list(mat_5, mat_5_transf)
mat_rem <- removeFalseLinksCircular(list(mat_5, mat_5_transf))
plot(graph_from_adjacency_matrix(mat_5), mode = "directed")
plot(graph_from_adjacency_matrix(mat_rem[[1]], mode = "directed"))
