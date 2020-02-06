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

## real world example 1
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

## real world example 2
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

## real world example 3 (same as mat_2 atm)
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
    ## new
    
    ## set mat_char to "" where mat_num == 0
    mat_char[which(mat_num == 0)] <- ""
    
    return(list(mat_num, mat_char))
}

## delete ingoing isotopes/adducts to M if there is no circular relation from 
## M to M+isot/adduct to M+transf to M+iso/adduct+trans, e.g. imagine
## 1) M and M+1 forming an "adduct_isotope+1", suppose there is a "transf_+162"
## from M+1 to another mass feature M2, but there is no link between M and M2-1,
## then remove the link M+1 and M2
## 2) M and M+1 forming an "adduct_isotope+1", suppose there is a "transf_+162"
## from M+1 to another mass feature M, AND there is a link between M and M2-1,
## then keep the the link M+1 and M2 (support that mass feature is a metabolite)
removeFalseLinksCircular <- function(mat_l) {
    
    mat_num <- mat_l[[1]]
    mat_char <- mat_l[[2]]
    
    ## first find relations with adducts: 
    ## (molecular ion suppose this is the vertex that has no ingoing adduct_ link)
    ## i.e. the column of M does not have a "adduct_"
    M <- apply(mat_char, 2, grep, pattern = "adduct_")
    
    ## iterate through M where M is not integer(0) --> these are adducts of
    ## the rows i 
    for (i in seq_len(length(M))) {
        
        if (length(M[[i]])) {
            ## if length(M[[i]]) != 0, M[[i]] is not the molecular ion:
            ## 1) --> check if it has transf_
            transf_M_og_ind <- grep(mat_char[i, ], pattern = "transf_")
            transf_M_og <- mat_char[i, transf_M_og_ind]
            transf_M_ig_ind <- grep(mat_char[, i], pattern = "transf_")
            transf_M_ig <- mat_char[transf_M_ig_ind, i]
            
            print(i)
            ## check only outgoing links now 
            if (length(transf_M_og_ind)) {
                
                tmp <- matrix(mat_char[M[[i]], ], byrow = FALSE, ncol = ncol(mat_char))
                
                ## 2) check if the molecular ion has transf_ 
                transf_MI_og_ind <- apply(tmp, 2, grepl, pattern = "transf_") ## outgoing transf links from molecular ions
                transf_MI_og_ind <- which(transf_MI_og_ind, arr.ind = TRUE)
                transf_MI_og <- tmp[transf_MI_og_ind]
                
                
                if (length(transf_MI_og)) {
                    ## 3) check if transf from molecular ion and M[[i]] form the same adduct than molecular ion and M[[i]]
                    ## if so, then keep the relation, otherwise remove the link from M[[i]] with the transf
                    ## for og_ind
                    if (is.matrix(transf_MI_og_ind)) {
                        transf_MI_og_ind <- transf_MI_og_ind[transf_MI_og %in% transf_M_og, ]    
                    } else {
                        transf_MI_og_ind <- transf_MI_og_ind[transf_MI_og %in% transf_M_og]
                    }
                    transf_MI_og <- transf_MI_og[transf_MI_og %in% transf_M_og]
                    ind_match <-  match(transf_MI_og, transf_M_og)
                    ind_rm <- which(is.na(ind_match))
                    ind_match <- which(!is.na(ind_match))
                    
                    ## for ind_match check if the linking feature of M has a link
                    ## with the linking feature of MI
                    if (length(ind_match) == 1) {
                        links_check <- transf_MI_og_ind[ind_match]    
                        ## linking feature of M
                        linkFeat_M <- mat_char[, transf_M_og_ind]
                        addiso_bw_linkFeat <- linkFeat_M[links_check]
                        
                    } else {
                        links_check <- transf_MI_og_ind[ind_match, ]    
                        ## linking feature of M
                        linkFeat_M <- mat_char[, transf_M_og_ind]
                        addiso_bw_linkFeat <- linkFeat_M[links_check[, "col"]]
                    }
                    
                    ## check if the link between M and MI is the same as 
                    ## the link between the linking feature of M and the linking
                    ## feature of MI (addiso_bw_linkFeat)
                    inds_link <- match(addiso_bw_linkFeat, mat_char[M[[i]], i])
                    ind_rm <- c(ind_rm, ind_match[which(is.na(inds_link))])

                    ## remove:
                    #if (length(ind_rm))
                    if (is.matrix(transf_M_og_ind)) {
                        ind_rm <- matrix(transf_M_og_ind[ind_rm, ], ncol = 2)    
                        mat_num[ind_rm] <- 0
                    } else {
                        mat_num[i, transf_M_og_ind[ind_rm]] <- 0
                    }
                } else {
                    ## if MI has no outgoing links, remove outgoing links to M
                    mat_num[i, transf_M_og_ind] <- 0
                }
                
            }
            ## check only ingoing links now
            if (length(transf_M_ig_ind)) {

                tmp <- matrix(mat_char[, M[[i]]], byrow = TRUE, ncol = ncol(mat_char)) ## was t()

                ## 2) check if the molecular ion(s) have ingoing transf_
                transf_MI_ig_ind <- apply(tmp, 2, grepl, pattern = "transf_") ## ingoing transf links from molecular ions
                transf_MI_ig_ind <- which(transf_MI_ig_ind, arr.ind = TRUE)
                transf_MI_ig <- tmp[transf_MI_ig_ind]

                if (length(transf_MI_ig)) {
                    ## 3) check, if transf from molecular ion and M[[i]] form the same adduct than molecular ion and M[[i]]
                    ## if so, then keep the relation, otherwise remove the link from M[[i]] with the transf
                    ## for og_ind
                    if (is.matrix(transf_MI_ig_ind)) {
                        transf_MI_ig_ind <- transf_MI_ig_ind[transf_MI_ig %in% transf_M_ig, ]    
                    } else {
                        transf_MI_ig_ind <- transf_MI_ig_ind[transf_MI_ig %in% transf_M_ig]
                    }
                    transf_MI_ig <- transf_MI_ig[transf_MI_ig %in% transf_M_ig]
                    ind_match <- match(transf_MI_ig, transf_M_ig)
                    ind_rm <- which(is.na(ind_match))
                    ind_match <- which(!is.na(ind_match))

                    ## for ind_match check if the linking feature of M has a link
                    ## with the linking feature of MI
                    if (length(ind_match) == 1) {
                        links_check <- transf_MI_ig_ind[ind_match]
                        ## linking features of M
                        linkFeat_M <- mat_char[transf_M_ig_ind, ]
                        addiso_bw_linkFeat <- linkFeat_M[links_check]

                    } else {
                        links_check <- transf_MI_ig_ind[ind_match, ]
                        ## linking features of M
                        linkFeat_M <- mat_char[transf_M_ig_ind, ]
                        addiso_bw_linkFeat <- linkFeat_M[links_check[, "col"]]
                    }

                    ## check if the link between M and MI is the same as the
                    ## link between the linking feature of M and the linking
                    ## feature of MI (addiso_bw_linkFeat)
                    inds_link <- match(addiso_bw_linkFeat, mat_char[M[[i]], i])
                    ind_rm <- c(ind_rm, ind_match[which(is.na(inds_link))])

                    # remove:
                    if (is.matrix(transf_M_ig_ind)) {
                        ind_rm <- matrix(transf_M_ig_ind[ind_rm, 2:1], ncol = 2)
                        mat_num[ind_rm] <- 0
                    } else {
                        mat_num[transf_M_ig_ind[ind_rm], i] <- 0
                    }
                } else {
                    ## if MI has no ingoing links, remove ingoing links to M
                    mat_num[transf_M_ig_ind, i] <- 0
                }

            }
        }
    }
    ## set mat_char to "" where mat_num == 0
    mat_char[which(mat_num == 0)] <- ""
    
    return(list(mat_num, mat_char))
   
}

removeFalseLinksCircular2 <- function(mat_l) {
    
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
                            for (k in 1:length(transf_MI_og_ind)) {
                                ## get adduct of M
                                adduct_ij <- mat_char[transf_MI_og_ind[k], ind_j]
                                if (adduct_ij != mat_char[M[[i]], i]) {
                                    ## remove link transf_M_og_ind[k] at row i
                                    mat_num[i, ind_j] <- 0
                                }   
                            }
                        ## if there are several MI do
                        } else {
                            for(k in 1:nrow(transf_MI_og_ind)) {
                                ## get adduct of M
                                adduct_ij <- mat_char[transf_MI_og_ind[k, "row"], ind_j]
                                if (adduct_ij != mat_char[i, M[[i]][transf_MI_og_ind[k, "row"]]]) {
                                    ## remove link transf_M_og_ind[k] at row i
                                    mat_num[i, ind_j] <- 0
                                }
                            } 
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
                            for (k in 1:length(transf_MI_ig_ind)) {
                                ## get adduct of M
                                adduct_ij <- mat_char[transf_MI_ig_ind[k], ind_j]
                                if (adduct_ij != mat_char[M[[i]], i]) {
                                    ## remove link transf_M_ig_ind[j] at row ind_j
                                    mat_num[ind_j, i] <- 0
                                }   
                            }
                        ## if there are several MI do
                        } else {
                            for(k in 1:nrow(transf_MI_ig_ind)) {
                                ## get adduct of M
                                adduct_ij <- mat_char[transf_MI_ig_ind[k, "row"], ind_j]
                                if (adduct_ij != mat_char[M[[i]][transf_MI_ig_ind[k, "row"]], i]) {
                                    ## remove link transf_M_ig_ind[k] at row ind_j
                                    mat_num[ind_j, i] <- 0
                                }
                            } 
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

mat_l <- list(mat_1, mat_1_transf)
mat_rem <- removeFalseLinksAdducts(list(mat_1, mat_1_transf))
mat_rem <- removeFalseLinksCircular2(mat_rem)
plot(graph_from_adjacency_matrix(mat_l[[1]], mode = "directed"))
plot(graph_from_adjacency_matrix(mat_rem[[1]], mode = "directed"))

mat_l <- list(mat_2, mat_2_transf)
mat_rem <- removeFalseLinksAdducts(mat_l)
mat_rem <- removeFalseLinksCircular2(mat_rem)
plot(graph_from_adjacency_matrix(mat_l[[1]], mode = "directed"))
plot(graph_from_adjacency_matrix(mat_rem[[1]], mode = "directed")) ## should only result in link between M_1->M+152_2, M_1-> M_2, M_1_M_3, M_2->M_3

#mat_rem <- removeFalseLinksAdducts(list(mat_3, mat_3_transf))
mat_l <- list(mat_3, mat_3_transf)
mat_rem <- removeFalseLinksCircular2(list(mat_3, mat_3_transf))
plot(graph_from_adjacency_matrix(mat_3), mode = "directed")
plot(graph_from_adjacency_matrix(mat_rem[[1]], mode = "directed"))
## raises warnings adduct_j %in% adduct_i

