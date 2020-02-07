library(igraph)


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

 

## do for all components
## iterate through rows
#' @name removeFalseLinksAdducts
#' @title
#' @description node with outgoing transf link should not have outgoing adduct link when
## the adduct-linking vertex has an ingoing adduct link of same type
#' @details 
#' Example: remove edges between M1+H and M2+Na, M3+H and M4+K, in general when there is 
## a link between an adduct and a molecular ion based on the following rule

## node with outgoing transf link should not have outgoing adduct link when
## the adduct-linking vertex has a ingoing adduct link

## node with ingoing transf link should not have ingoing adduct link when the 
## adduct-linking vertex has an outgoing adduct link
#' @author Thomas Naake, <thomasnaake@@googlemail.com>
#' @example
#' 
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
#' @author Thomas Naake <thomasnaake@@googlemail.com>
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
                    for(j in seq_len(length(transf_M_og_ind))) {
                        ind_j <- transf_M_og_ind[j]
                        ## if there is only one MI do
                        if (!is.matrix(transf_MI_og_ind)) {
                            remove <- TRUE
                            for (k in seq_len(length(transf_MI_og_ind))) {
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
                            for(k in seq_len(nrow(transf_MI_og_ind))) {
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
                    for(j in seq_len(length(transf_M_ig_ind))) {
                    
                        ind_j <- transf_M_ig_ind[j]
                        ## if there is only one MI do
                        if (!is.matrix(transf_MI_ig_ind)) {
                            remove <- TRUE
                            for (k in seq_len(length(transf_MI_ig_ind))) {
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
                            for(k in seq_len(nrow(transf_MI_ig_ind))) {
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
