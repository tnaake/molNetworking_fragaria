## define test examples

library(igraph)
## example 1
mat_1 <- matrix(c(
    c(0, 1, 0, 0, 0, 1), ## M1+H
    c(0, 0, 0, 0, 0, 0), ## M1+Na
    c(1, 0, 0, 1, 0, 0), ## M2+H
    c(0, 1, 0, 0, 0, 0), ## M2+Na
    c(0, 1, 0, 0, 0, 1), ## M3+H
    c(0, 0, 0, 0, 0, 0)), ## M3+Na
    byrow = TRUE, ncol = 6, nrow = 6)
mat_1_transf <- matrix(c(
    c("",              "adduct_Na",     "", "",          "", "transf_gluc_F"), ## M1+H
    c("",              "",              "", "",          "", ""), ## M1+Na
    c("transf_gluc_T", "",              "", "adduct_Na", "", ""), ## M2+H
    c("",              "transf_gluc_T", "", "",          "", ""), ## M2+Na
    c("",              "transf_gluc_F", "",              "", "", "adduct_Na"), ## M3+H
    c("",              "",              "", "",          "", "")), ## M3+Na
    byrow = TRUE, ncol = 6, nrow = 6)
colnames(mat_1) <- c("M1+H", "M1+Na", "M2+H", "M2+Na", "M3+H", "M3+Na")
rownames(mat_1) <- colnames(mat_1)
rownames(mat_1_transf) <- colnames(mat_1_transf) <- colnames(mat_1)

net_mat_1 <- graph_from_adjacency_matrix(mat_1, mode = "directed")
plot(net_mat_1)

## example 2
mat_2 <- matrix(c(
    c(0, 1, 1, 0, 1, 0, 0, 0, 0), ## "M_1"
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
colnames(mat_4) <- c("M_1", "M+1_1", "M+2_1", "M+45_5", "M+152_2", "M-152_3", "M+1-152_4", "M+1+152_2", "M+2+152_2")
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

## example 6
mat_6 <- matrix(c(
    c(0, 1, 1, 0, 1, 0, 0, 0, 0), ## "M_1"
    c(0, 0, 1, 1, 0, 0, 0, 1, 0), ## "M+1_1"
    c(0, 0, 0, 0, 0, 0, 0, 0, 0), ## "M+2_1"
    c(0, 0, 0, 0, 0, 0, 0, 0, 0), ## "M+1+45_1"
    c(0, 0, 0, 0, 0, 0, 0, 0, 0), ## "M+152_2"
    c(0, 1, 0, 0, 1, 0, 0, 0, 0), ## "M+1-152_3"
    c(0, 1, 0, 0, 0, 0, 0, 0, 0), ## "M+1-152_4"
    c(0, 0, 0, 0, 0, 0, 0, 0, 0), ##"M+1+152_5"
    c(0, 0, 1, 0, 0, 0, 0, 0, 0)), ## "M+2-152_6"
    byrow = TRUE, ncol = 9, nrow = 9)

mat_6_transf <- matrix(c(
    c("", "adduct_+1",  "adduct_+2",  "",          "transf_152", "", "", "",           ""), ## "M_1"
    c("", "",           "adduct_+1",  "transf_45", "",           "", "", "transf_152", ""), ## "M+1_1"
    c("", "",           "",           "",          "",           "", "", "",           ""), ## "M+2_1"
    c("", "",           "",           "",          "",           "", "", "",           ""), ## "M+1+45_1"
    c("", "",           "",           "",          "",           "", "", "",           ""), ## "M+152_2"
    c("", "transf_152", "",           "",          "transf_304", "", "", "",           ""), ## "M+1-152_3"
    c("", "transf_152", "",           "",          "",           "", "", "",           ""), ## "M+1-152_4"
    c("", "",           "",           "",          "",           "", "", "",           ""), ##"M+1+152_5"
    c("", "",           "transf_152", "",          "",           "", "", "",           "")), ## "M+2-152_6"
    byrow = TRUE, ncol = 9, nrow = 9)
colnames(mat_6) <- c("M_1", "M+1_1", "M+2_1", "M+1+45_1", "M+152_2", "M+1-152_3", "M+1-152_4", "M+1+152_5", "M+2-152_6")
rownames(mat_6) <- colnames(mat_6)
colnames(mat_6_transf) <- rownames(mat_6_transf) <- colnames(mat_6)
net_mat_6 <- graph_from_adjacency_matrix(mat_6, mode = "directed")
plot(net_mat_6)

## apply functions and perform tests

## example 1
mat_l <- list(mat_1, mat_1_transf)
mat_rem <- removeFalseLinksAdducts(mat_l)
mat_rem <- removeFalseLinksCircular(mat_rem)
plot(graph_from_adjacency_matrix(mat_l[[1]], mode = "directed"))
plot(graph_from_adjacency_matrix(mat_rem[[1]], mode = "directed"))

## should result in M1+H->M1+Na, M2+H->M2+Na, M2+H->M1+H, M2+Na->M1+Na, M3+H->M3+Na
mat_res_num_1 <- matrix(c(
    c(0, 1, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0),
    c(1, 0, 0, 1, 0, 0),
    c(0, 1, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 1),
    c(0, 0, 0, 0, 0, 0)),
    byrow = TRUE, ncol = 6, nrow = 6)

mat_res_char_1 <- matrix(c(
    c("",              "adduct_Na",     "",   "",          "",   ""),         
    c("",              "",              "",   "",          "",   ""),         
    c("transf_gluc_T", "",              "",   "adduct_Na", "",   ""),         
    c("",              "transf_gluc_T", "",   "",          "",   ""),       
    c("",              "",              "",   "",          "",   "adduct_Na"),
    c("",              "",              "",   "",          "",   "")),
    byrow = TRUE, ncol = 6, nrow = 6)

expect_true(all(mat_rem[[1]] == mat_res_num_1))
expect_true(all(mat_rem[[2]] == mat_res_char_1))

# example 2 
mat_l <- list(mat_2, mat_2_transf)
mat_rem <- removeFalseLinksAdducts(mat_l)
mat_rem <- removeFalseLinksCircular(mat_rem)
plot(graph_from_adjacency_matrix(mat_l[[1]], mode = "directed"))
plot(graph_from_adjacency_matrix(mat_rem[[1]], mode = "directed")) 

## should only result in link between M_1->M+152_2, M_1-> M_2, M_1_M_3, M_2->M_3, M+1-152_3->M+152_2
mat_res_num_2 <- matrix(c(
    c(0, 1, 1, 0, 1, 0, 0, 0, 0),
    c(0, 0, 1, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 1, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 0)),
    byrow = TRUE, ncol = 9, nrow = 9)

mat_res_char_2 <- matrix(c(
    c("", "adduct_+1", "adduct_+2", "", "transf_152", "", "", "", ""), 
    c("", "",          "adduct_+1", "", "",           "", "", "", ""),   
    c("", "",          "",          "", "",           "", "", "", ""),     
    c("", "",          "",          "", "",           "", "", "", ""),     
    c("", "",          "",          "", "",           "", "", "", ""),     
    c("", "",          "",          "", "transf_304", "", "", "", ""),     
    c("", "",          "",          "", "",           "", "", "", ""),     
    c("", "",          "",          "", "",           "", "", "", ""),     
    c("", "",          "",          "", "",           "", "", "", "")),
    byrow = TRUE, ncol = 9, nrow = 9)

expect_true(all(mat_rem[[1]] == mat_res_num_2))
expect_true(all(mat_rem[[2]] == mat_res_char_2))

## example 3
mat_l <- list(mat_3, mat_3_transf)
mat_rem <- removeFalseLinksAdducts(mat_l)
mat_rem <- removeFalseLinksCircular(mat_rem)
plot(graph_from_adjacency_matrix(mat_3), mode = "directed")
plot(graph_from_adjacency_matrix(mat_rem[[1]], mode = "directed"))

## should result in links M_1->M+1_1, M_1->M+2_1, M+1_1->M+2_1, 
## M+152_2->M+1+152_2, M+152_2->M+2+152_2, M+1+152_2->M+2+152_2, M_1->M+152_2, 
## M+1_1->M+1+152_2, M+2_1->M+2+152_2
mat_res_num_3 <- matrix(c(
    c(0, 1, 1, 0, 1, 0, 0, 0, 0),
    c(0, 0, 1, 0, 0, 0, 0, 1, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 1),
    c(0, 0, 0, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 1, 1),
    c(0, 0, 0, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 1),
    c(0, 0, 0, 0, 0, 0, 0, 0, 0)),
    byrow = TRUE, ncol = 9, nrow = 9)
mat_res_char_3 <- matrix(c(
    c("", "adduct_+1", "adduct_+2", "", "transf_152", "", "", "",           ""),      
    c("", "",          "adduct_+1", "", "",           "", "", "transf_152", ""),          
    c("", "",          "",          "", "",           "", "", "",           "transf_152"),
    c("", "",          "",          "", "",           "", "", "",           ""),
    c("", "",          "",          "", "",           "", "", "adduct_+1",  "adduct_+2"),
    c("", "",          "",          "", "",           "", "", "",           ""),       
    c("", "",          "",          "", "",           "", "", "",           ""),        
    c("", "",          "",          "", "",           "", "", "",           "adduct_+1"),
    c("", "",          "",          "", "",           "", "", "",           "")),
    byrow = TRUE, ncol = 9, nrow = 9)

expect_true(all(mat_rem[[1]] == mat_res_num_3))
expect_true(all(mat_rem[[2]] == mat_res_char_3))

## example 4
mat_l <- list(mat_4, mat_4_transf)
mat_rem <- removeFalseLinksAdducts(mat_l)
mat_rem <- removeFalseLinksCircular(mat_rem)
plot(graph_from_adjacency_matrix(mat_4), mode = "directed")
plot(graph_from_adjacency_matrix(mat_rem[[1]], mode = "directed"))

## should result in M_1->M+1_1, M_1->M+2_1, M+1_1->M+2_1, M+152_2->M+1+152_2, 
## M+152_2->M+2+152_2, M+1+152_2->M+2+152_2, M_1->M+152_2, M+1_1->M+1+152_2, 
## M+2_1->M+2+152_2, M_1->M+45_5, M-152_3->M_1
mat_res_num_4 <- matrix(c(
    c(0, 1, 1, 1, 1, 0, 0, 0, 0),
    c(0, 0, 1, 0, 0, 0, 0, 1, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 1),
    c(0, 0, 0, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 1, 1),
    c(1, 0, 0, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 1),
    c(0, 0, 0, 0, 0, 0, 0, 0, 0)),
    byrow = TRUE, ncol = 9, nrow = 9)

mat_res_char_4 <- matrix(c(
    c("",           "adduct_+1", "adduct_+2", "transf_45", "transf_152", "", "", "",           ""),
    c("",           "",          "adduct_+1", "",          "",           "", "", "transf_152", ""),
    c("",           "",          "",          "",          "",           "", "", "",           "transf_152"),
    c("",           "",          "",          "",          "",           "", "", "",           ""),
    c("",           "",          "",          "",          "",           "", "", "adduct_+1",  "adduct_+2"),
    c("transf_152", "",          "",          "",          "",           "", "", "",           ""),
    c("",           "",          "",          "",          "",           "", "", "",           ""),        
    c("",           "",          "",          "",          "",           "", "", "",           "adduct_+1"),
    c("",           "",          "",          "",          "",           "", "", "",           "")),
    byrow = TRUE, ncol = 9, nrow = 9)

expect_true(all(mat_rem[[1]] == mat_res_num_4))
expect_true(all(mat_rem[[2]] == mat_res_char_4))

## example 5
mat_l <- list(mat_5, mat_5_transf)
mat_rem <- removeFalseLinksAdducts(mat_l)
mat_rem <- removeFalseLinksCircular(mat_rem)
plot(graph_from_adjacency_matrix(mat_5), mode = "directed")
plot(graph_from_adjacency_matrix(mat_rem[[1]], mode = "directed"))

mat_res_num_5 <- matrix(c(
    c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0),
    c(0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
    c(0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0),
    c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)), 
    byrow = TRUE, ncol = 11, nrow = 11)

mat_res_char_5 <- matrix(c(
    c("",           "adduct_+1", "adduct_+2", "transf_45", "transf_152", "",          "", "",           "",           "",          ""),
    c("",           "",          "adduct_+1", "",          "",           "",          "", "transf_152", "",           "",          ""),        
    c("",           "",          "",          "",          "",           "",          "", "",           "transf_152", "transf_45", ""),        
    c("",           "",          "",          "",          "",           "",          "", "",           "",           "adduct_+2", ""),        
    c("",           "",          "",          "",          "",           "",          "", "adduct_+1",  "adduct_+2",  "",          ""),        
    c("transf_152", "",          "",          "",          "",           "",          "", "",           "",           "",          ""),        
    c("",           "",          "",          "",          "",           "",          "", "",           "",           "",          ""),      
    c("",           "",          "",          "",          "",           "",          "", "",           "adduct_+1",  "",          ""),        
    c("",           "",          "",          "",          "",           "",          "", "",           "",           "",          ""),      
    c("",           "",          "",          "",          "",           "",          "", "",           "",           "",          ""),      
    c("",           "",          "",          "",          "",           "transf_45", "", "",           "",           "",          "")),
    byrow = TRUE, ncol = 11, nrow = 11)

expect_true(all(mat_res_num_5 == mat_rem[[1]]))
expect_true(all(mat_res_char_5 == mat_rem[[2]]))

## example 6
mat_l <- list(mat_6, mat_6_transf)
mat_rem <- removeFalseLinksAdducts(mat_l)
mat_rem <- removeFalseLinksCircular(mat_rem)
plot(graph_from_adjacency_matrix(mat_6), mode = "directed")
plot(graph_from_adjacency_matrix(mat_rem[[1]], mode = "directed"))

## results in M1-> M1+1, M1-> M2+1, M1-> M+152, M+1-152->M+152
mat_res_num_6 <- matrix(c(
    c(0, 1, 1, 0, 1, 0, 0, 0, 0),
    c(0, 0, 1, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 1, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 0, 0, 0, 0)), 
    byrow = TRUE, ncol = 9, nrow = 9)

## transf_304 does not make chemically sense, it is just here to check the 
## correctness of the function
mat_res_char_6 <- c(
    c("", "adduct_+1", "adduct_+2", "", "transf_152", "", "", "", ""),       
    c("", "",          "adduct_+1", "", "",           "", "", "", ""),
    c("", "",          "",          "", "",           "", "", "", ""),      
    c("", "",          "",          "", "",           "", "", "", ""),      
    c("", "",          "",          "", "",           "", "", "", ""),   
    c("", "",          "",          "", "transf_304", "", "", "", ""),      
    c("", "",          "",          "", "",           "", "", "", ""),    
    c("", "",          "",          "", "",           "", "", "", ""),      
    c("", "",          "",          "", "",           "", "", "", ""),
    byrow = TRUE, ncol = 9, nrow = 9)
