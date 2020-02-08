#' @name removeFalseLinksAdducts
#'
#' @title Remove false links based on adduct information
#'
#' @description 
#' Node with outgoing `transf_` link should not have outgoing
#' adduct link when the adduct-linking vertex has an ingoing adduct link of the
#' same type. The function `removeFalseLinksAdducts` removes these links that
#' have the links of this kind, e.g. a wrong link that is by change denoted 
#' as a transformation, however the transformation is actually not a 
#' transformation (the m/z difference is the effect of the difference in 
#' isotopes + remaining mass difference). 
#'
#' @details
#' An isotope adduct in this sense above is not an adduct, even though it is 
#' codified as `adduct_isotopic`. All isotopes (`adduct_isotopic[1-3]``) will
#' be ignored by `removeFalseLinksAdducts` to avoid that true links are removed.
#' The user has to assure that isotopes have to be encoded as `adduct_isotopic`,
#' otherwise it is possible that `removeFalseLinksAdducts` will delete true 
#' links.
#'
#' The function `removeFalseLinksAdducts` will remove a link, when there is a
#' transformation link (`transf_`) between an adduct and a molecular ion,
#' e.g. `removeFalseLinksAdducts` removes edges between the mass features 
#' 1) `M1+H` and `M2+Na`, or 
#' 2) `M3+H` and `M4+K`,
#' if there is a transformation link between `M1+H` and `M2+Na`, or `M3+H` and
#' `M4+K`.
#'
#' The following rules are applied to remove the links:
#'
#' 1) A node with an outgoing `transf_` link should not have an outgoing adduct
#' link, if the adduct-linking node has a ingoing adduct link.
#' 2) A node with an ingoing `transf_` link should not have ingoing adduct link,
#' if the adduct-linking vertex has an outgoing adduct link.
#'
#' @author Thomas Naake, <thomasnaake@@googlemail.com>
#' @example
#' mat_1 <- matrix(c(
#'     c(0, 1, 0, 0, 0, 1), ## M1+H
#'     c(0, 0, 0, 0, 0, 0), ## M1+Na
#'     c(1, 0, 0, 1, 0, 0), ## M2+H
#'     c(0, 1, 0, 0, 0, 0), ## M2+Na
#'     c(0, 1, 0, 0, 0, 1), ## M3+H
#'     c(0, 0, 0, 0, 0, 0)), ## M3+Na
#'     byrow = TRUE, ncol = 6, nrow = 6)
#'
#' ## transf_gluc_T denotes a true underlying transformation (glucosylation)
#' ## transf_gluc_T denotes an incorrect underlying transformation 
#' ## (glycosylation)
#' mat_1_transf <- matrix(c(
#'     c("",              "adduct_Na",     "", "",          "", "transf_gluc_F"),
#'     c("",              "",              "", "",          "", ""), 
#'     c("transf_gluc_T", "",              "", "adduct_Na", "", ""), 
#'     c("",              "transf_gluc_T", "", "",          "", ""), 
#'     c("",              "transf_gluc_F", "",  "",         "", "adduct_Na"),
#'     c("",              "",              "", "",          "", "")), 
#'     byrow = TRUE, ncol = 6, nrow = 6)
#'     colnames(mat_1) <- c("M1+H", "M1+Na", "M2+H", "M2+Na", "M3+H", "M3+Na")
#'     rownames(mat_1) <- colnames(mat_1)
#'     rownames(mat_1_transf) <- colnames(mat_1_transf) <- colnames(mat_1)
#'
#' removeFalseLinksAdducts(list(mat_1, mat_1_transf))
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
                adduct_j <- cols_j[grep(cols_j, pattern = "adduct_")]
                adduct_j <- adduct_j[-grep(adduct_j, pattern = "duct_isotopic")]
                if (length(adduct_j) > 0)
                    if (any(adduct_j %in% adduct_i)) {
                        mat_num[i, j] <- 0 ## remove in mat_num[i, j]
                    }
            }
        }
        
        cols_i <- mat_char[, i]
        adduct_i <- cols_i[grep(cols_i, pattern = "adduct_")]
        
        ## if any vertex has ingoing adduct link, check all rows where cols_i
        ## has transf if outgoing adduct link exists for adduct-linking feature
        if (length(adduct_i) > 0) { 
            
            inds_col <- grep(cols_i, pattern = "transf_")
            
            for (j in inds_col) {
                rows_j <- mat_char[j, ]
                adduct_j <- rows_j[grep(rows_j, pattern = "adduct")]
                adduct_j <- adduct_j[-grep(adduct_j, pattern = "duct_isotopic")]
                if (length(adduct_j) > 0)
                    if (any(adduct_j %in% adduct_i)) {
                        mat_num[j, i] <- 0 ## remove in mat_num[j, i]
                    }
            }
        }
        
    }
    
    ## set mat_char to "" where mat_num == 0
    mat_char[which(mat_num == 0)] <- ""
    
    return(list(mat_num, mat_char))
}

#' @name removeFalseLinksCircular
#' @title Remove in- and outgoing links based on their relation to the molecular
#' ion
#'
#' @description
#' The function `removeFalseLinksCircular` deletes ingoing 
#' and outgoing links of isotopes/adducts to M if there is no circular 
#' relation from the molecular ion M to M+isotope/adduct to M+transf to 
#' M+iso/adduct+transformation.
#' 
#' @details
#' Examples:
#' 1) M and M+1 forming an `adduct_isotopic+1`, suppose there is a `transf_+162`
#' from M+1 to another mass feature M2, but there is no link between M and M2-1,
#' then remove the link M+1 and M2.
#'
#' 2) M and M+1 forming an `adduct_isotopic+1`, suppose there is a `transf_+162`
#' from M+1 to another mass feature M, AND there is a link between M and M2-1,
#' then keep the the link M+1 and M2 (support that mass feature is a metabolite)
#'
#' @author Thomas Naake <thomasnaake@@googlemail.com>
#'
#' @example
#' mat_2 <- matrix(c(
#'     c(0, 1, 1, 0, 1, 0, 0, 0, 0),
#'     c(0, 0, 1, 1, 0, 0, 0, 1, 0),
#'     c(0, 0, 0, 0, 0, 0, 0, 0, 0),
#'     c(0, 0, 0, 0, 0, 0, 0, 0, 0),
#'     c(0, 0, 0, 0, 0, 0, 0, 0, 0),
#'     c(0, 1, 0, 0, 0, 0, 0, 0, 0),
#'     c(0, 1, 0, 0, 0, 0, 0, 0, 0),
#'     c(0, 0, 0, 0, 0, 0, 0, 0, 0),
#'     c(0, 0, 1, 0, 0, 0, 0, 0, 0)),
#'     byrow = TRUE, ncol = 9, nrow = 9)
#'     
#' mat_2_transf <- matrix(c(
#'     c("", "adduct_isotopic+1", "adduct_isotopic+2", "",          "transf_152", "", "", "",           ""),
#'     c("", "",                  "adduct_isotopic+1", "transf_45", "",           "", "", "transf_152", ""),
#'     c("", "",                  "",                  "",          "",           "", "", "",           ""),
#'     c("", "",                  "",                  "",          "",           "", "", "",           ""),
#'     c("", "",                  "",                  "",          "",           "", "", "",           ""),
#'     c("", "transf_152",        "",                  "",          "",           "", "", "",           ""),
#'     c("", "transf_152",        "",                  "",          "",           "", "", "",           ""),
#'     c("", "",                  "",                  "",          "",           "", "", "",           ""),
#'     c("", "",                  "transf_152",        "",          "",           "", "", "",           "")),
#'     byrow = TRUE, ncol = 9, nrow = 9)
#' colnames(mat_2) <- c("M_1", "M+1_1", "M+2_1", "M+1+45_1", 
#'     "M+152_2", "M+1-152_3", "M+1-152_4", "M+1+152_5", "M+2-152_6")
#' rownames(mat_2) <- colnames(mat_2)
#' colnames(mat_2_transf) <- rownames(mat_2_transf) <- colnames(mat_2)
#' 
#' removeFalseLinksCircular(list(mat_2, mat_2_transf)) 
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
