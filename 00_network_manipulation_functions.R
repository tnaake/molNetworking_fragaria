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


mat_test <- matrix(c(c(0, 1, 0, 0, 0, 1), ## M1+H
                     c(0, 0, 0, 0, 0, 0), ## M1+Na
                     c(1, 0, 0, 1, 0, 0), ## M2+H
                     c(0, 1, 0, 0, 0, 0), ## M2+Na
                     c(0, 1, 0, 0, 0, 1), ## M3+H
                     c(0, 0, 0, 0, 0, 0)), ## M3+Na
                   byrow = TRUE, ncol = 6, nrow = 6)
mat_test_transformation <- matrix(c(c("", "adduct_Na", "", "", "", "transf_gluc_F"), ## M1+H
                                    c("", "", "", "", "", ""), ## M1+Na
                                    c("transf_gluc_T", "", "", "adduct_Na", "", ""), ## M2+H
                                    c("", "transf_gluc_T", "", "", "", ""), ## M2+Na
                                    c("", "transf_gluc_F", "", "", "", "adduct_Na"), ## M3+H
                                    c("", "", "", "", "", "")), ## M3+Na
                                  byrow = TRUE, ncol = 6, nrow = 6)
colnames(mat_test) <- c("M1+H", "M1+Na", "M2+H", "M2+Na", "M3+H", "M3+Na")
rownames(mat_test) <- colnames(mat_test)

## real world example 1
mat_test2 <- matrix(c(c(0, 1, 1, 0, 1, 0, 0, 0, 0), ## "M_1"
                      c(0, 0, 1, 1, 0, 0, 0, 1, 0), ## "M+1_1"
                      c(0, 0, 0, 0, 0, 0, 0, 0, 0), ## "M+2_1"
                      c(0, 0, 0, 0, 0, 0, 0, 0, 0), ## "M+1+45_1"
                      c(0, 0, 0, 0, 0, 0, 0, 0, 0), ## "M+152_2"
                      c(0, 1, 0, 0, 0, 0, 0, 0, 0), ## "M+1-152_3"
                      c(0, 1, 0, 0, 0, 0, 0, 0, 0), ## "M+1-152_4"
                      c(0, 0, 0, 0, 0, 0, 0, 0, 0), ##"M+1+152_5"
                      c(0, 0, 1, 0, 0, 0, 0, 0, 0)), ## "M+2-152_6"
                    byrow = TRUE, ncol = 9, nrow = 9)
mat_test2_transformation <- matrix(c(
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
colnames(mat_test2) <- c("M_1", "M+1_1", "M+2_1", "M+1+45_1", "M+152_2", "M+1-152_3", "M+1-152_4", "M+1+152_5", "M+2-152_6")
rownames(mat_test2) <- colnames(mat_test2)
colnames(mat_test2_transformation) <- rownames(mat_test2_transformation) <- colnames(mat_test2)
net_mat2 <- graph_from_adjacency_matrix(mat_test2, mode = "directed")
plot(net_mat2)

## real world example 3
mat_test3 <- matrix(c(c(0, 1, 1, 0, 1, 0, 0, 0, 0), ## "M_1"
                      c(0, 0, 1, 1, 0, 0, 0, 1, 0), ## "M+1_1"
                      c(0, 0, 0, 0, 0, 0, 0, 0, 0), ## "M+2_1"
                      c(0, 0, 0, 0, 0, 0, 0, 0, 0), ## "M+1+45_1"
                      c(0, 0, 0, 0, 0, 0, 0, 0, 0), ## "M+152_2"
                      c(0, 1, 0, 0, 0, 0, 0, 0, 0), ## "M+1-152_3"
                      c(0, 1, 0, 0, 0, 0, 0, 0, 0), ## "M+1-152_4"
                      c(0, 0, 0, 0, 0, 0, 0, 0, 0), ##"M+1+152_5"
                      c(0, 0, 1, 0, 0, 0, 0, 0, 0)), ## "M+2-152_6"
                    byrow = TRUE, ncol = 9, nrow = 9)
mat_test3_transformation <- matrix(c(
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
colnames(mat_test3) <- c("M_1", "M+1_1", "M+2_1", "M+1+45_1", "M+152_2", "M+1-152_3", "M+1-152_4", "M+1+152_5", "M+2-152_6")
rownames(mat_test3) <- colnames(mat_test3)
colnames(mat_test3_transformation) <- rownames(mat_test3_transformation) <- colnames(mat_test3)
net_mat3 <- graph_from_adjacency_matrix(mat_test3, mode = "directed")
plot(net_mat3)

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
removeFalsLinksCircular <- function(mat_l) {
    
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
                
                mat <- matrix(mat_char[M[[i]], ], byrow = FALSE, ncol = ncol(mat_char))
                colnames(mat) <- colnames(mat_char)
                
                ## 2) check if the molecular ion has transf_ 
                transf_MI_og_ind <- apply(mat, 2, grepl, pattern = "transf_") ## outgoing transf links from molecular ions
                transf_MI_og_ind <- which(transf_MI_og_ind, arr.ind = TRUE)
                transf_MI_og <- mat_char[M[[i]], ][transf_MI_og_ind]
                
                #if (length(transf_MI_og)) {
                ## 3) check if transf from molecular ion and M[[i]] form the same adduct than molecular ion and M[[i]]
                ## if so, then keep the relation, otherwise remove the link from M[[i]] with the transf
                ## for og_ind
                #if (length(transf_MI_og_ind)) {
                
                ind_match <- transf_M_og %in% transf_MI_og #match(transf_M_og, transf_MI_og) ## if NA, then it means that there is no shared transf_ --> remove
                ind_rm <- which(!ind_match) ## which(is.na(ind_match))
                ind_match <- which(ind_match) ##which(!is.na(ind_match))
                
                ## for ind_match check if the linking feature of M has a link
                ## with the linking feature of MI
                addiso_bw_linkFeat <- mat_char[transf_M_og_ind[ind_match], ][transf_MI_og_ind]
                
                ## mat_char[M[[i]], i] is the type of link between MI and M
                #if (addiso_bw_linkFeat %in% mat_char[M[[i]], i]) {
                ind_rm <- c(ind_rm, ind_match[!addiso_bw_linkFeat %in% mat_char[M[[i]], i]])
                
                ## remove:
                mat_num[i, transf_M_og_ind[ind_rm]] <- 0
                
            }
            ## check only ingoing links now
            if (length(transf_M_ig_ind)) {
                ## 2) check if the molecular ion(s) have transf_ 
                tmp <- t(matrix(mat_char[, M[[i]]], byrow = TRUE, ncol = ncol(mat_char)))
                
                transf_MI_ig_ind <- apply(tmp, 2, grepl, pattern = "transf_") ## ingoing transf links from molecular ions
                transf_MI_ig_ind <- which(transf_MI_ig_ind, arr.ind = TRUE)
                transf_MI_ig <- mat_char[, M[[i]]][transf_MI_ig_ind]
                
                ## 3) check, if transf from molecular ion and M[[i]] form the same adduct than molecular ion and M[[i]]
                ## if so, then keep the relation, otherwise remove the link from M[[i]] with the transf
                ## for og_ind
                ind_match <- transf_M_ig %in% transf_MI_ig ## if NA, then it means that there is no shared transf_ --> remove
                ind_rm <- which(!ind_match)
                ind_match <- which(ind_match)
                
                ## for ind_match check if the linking feature of M has a link
                ## with the linking feature of MI
                addiso_bw_linkFeat <- mat_char[, transf_M_ig_ind[ind_match]][transf_MI_ig_ind]
                
                ## mat_char[M[[i]], i] is the type of link between MI and M
                ## if addiso_bw_linkFeat %in% mat_char[M[[i]], i]
                ind_rm <- c(ind_rm, ind_match[!addiso_bw_linkFeat %in% mat_char[M[[i]], i]])
                
                ## remove:
                mat_num[transf_M_ig_ind[ind_rm], i] <- 0
            }
        }
    }
    ## set mat_char to "" where mat_num == 0
    mat_char[which(mat_num == 0)] <- ""
    
    return(list(mat_num, mat_char))
    
}


mat_rem <- removeFalseLinksAdducts(list(mat_test2, mat_test2_transformation))

mat_rem <- removeFalsLinksCircular(mat_rem)
plot(graph_from_adjacency_matrix(mat_rem[[1]], mode = "directed")) ## should only result in link between M_1 and M+152_2
mat_rem <- removeFalseLinksAdducts(list(mat_test, mat_test_transformation))
mat_rem <- removeFalsLinksCircular(mat_rem)
plot(graph_from_adjacency_matrix(mat_rem[[1]], mode = "directed"))
## raises warnings adduct_j %in% adduct_i

