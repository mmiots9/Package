mrt.par.2aov <- function(x, y, z, type.of.int = NA){

  # 1. controlli
    # x, y, z numeriche
    if (!is.numeric(x)) {
     stop("x must be a numeric vector")
    }

  # valid type of int
  if (is.na(type.of.int)){
    stop("type.of.int not specified with no default")}
  if (type.of.int != "+" && type.of.int != "*") {
    stop("type.of.int must be '+' or '*'" )
  }

  # 2. shapiro categorie
  shapy <- tapply(x, y, shapiro.test)
  shapz <- tapply(x, z, shapiro.test)

  # 3. Inizializzo vettori utili
  categy <- c()
  fatty  <- c()
  wy     <- c()
  pvaly  <- c()
  categz <- c()
  fattz  <- c()
  wz     <- c()
  pvalz  <- c()

  # 4. Creo nomi categorie ingresso
  caty <- deparse1(substitute(y))
  catz <- deparse1(substitute(z))

  # 5. creo loop per estrarre cose importanti da shapy e shapz
  # shapy
  for (i in 1:length(shapy)){
    categy[i] <- caty
    fatty[i]  <- names(shapy[i])
    wy[i]     <- shapy[[i]]$statistic
    pvaly[i]  <- shapy[[i]]$p.value
  }

  # shapz
  for (i in 1:length(shapz)){
    categz[i] <- catz
    fattz[i]  <- names(shapz[i])
    wz[i]     <- shapz[[i]]$statistic
    pvalz[i]  <- shapz[[i]]$p.value
  }

  # 6. creo dataframe
  shcat <- data.frame(categoric = c(categy, categz), fattore = c(fatty, fattz), W = c(wy, wz), pval = c(pvaly, pvalz))

  #7. Bartlett
    # categorie
    bary <- bartlett.test(x, y)
    barz <- bartlett.test(x, z)

    # estraggo valori
    barcat    <- caty
    ksq       <- bary$statistic
    bardf     <- bary$parameter
    barpv     <- bary$p.value
    barcat[2] <- catz
    ksq[2]    <- barz$statistic
    bardf[2]  <- barz$parameter
    barpv[2]  <- barz$p.value

    # se interazione
    if(type.of.int == "*"){
      bartint <- bartlett.test(x, interaction(y,z))
      barcat[3] <- "interaction"
      ksq[3]    <- bartint$statistic
      bardf[3]  <- bartint$parameter
      barpv[3]  <- bartint$p.value
    }
    else {
      barcat[3] <- NA
      ksq[3]    <- NA
      bardf[3]  <- NA
      barpv[3]  <- NA
      }

    # creo dataframe
    bartot <- data.frame(barcat, ksq, bardf, barpv)
    names(bartot) <- c("categoric", "k-squared", "df", "pval")
    bartot <- bartot[-3,]


  # 8. Shapiro se *
    if(type.of.int == "*"){

    # creo lista con tutti i risultati
    mat <- tapply(x, list(y, z), shapiro.test)
    nm1 <- outer(rownames(mat), colnames(mat), FUN=paste)
    nm2 <- setNames(c(mat), nm1)
    nom <- strsplit(names(nm2), split = "") # creo lista con fattori categorie

    # creo vettori interesse
    prima   <- c()
    seconda <- c()
    wint    <- c()
    pvalint <- c()

    # loop per riempire vettori
    for(i in 1:length(nm2)){
      prima[i]   <- nom[[i]][1]
      seconda[i] <- nom[[i]][3]
      wint[i]    <- nm2[[i]]$statistic
      pvalint[i] <- nm2[[i]]$p.value
    }

    # creo dataframe
    shint <- data.frame(prima, seconda, wint, pvalint)
    names(shint) <- c(caty, catz, "W", "pval")

    # return
    res <- list(shcat, shint, bartot)
    names(res) <- c("Shapiro categorie", "Shapiro interazione", "Bartlett")
    }
    else {
      res <- c(list(shcat, bartot))
      names(res) <- c("Shapiro categorie", "Bartlett")
      }
    return(res)
}

vicenzi_lab_pcr = function(){

  library(readxl)
  library(svGUI)
  library(svDialogs)
  library(dplyr)
  library(xlsx)
  options(digits = 15)

  # 1. selezionare il file
  file_to_use <- file.choose()
  name_end    <- strsplit(file_to_use, '/')
  name_end    <- strsplit(name_end[[1]][length(name_end[[1]])],'[.]')
  name_end    <- name_end[[1]][1]

  # 2. import table
  my_table = read_excel(file_to_use)

  # 3. trovare posizione NA prima colonna
  nas = which(is.na(my_table[,1]))

  # 4. ridimensiono la tabella
  my_table = my_table[(nas[1]+1):(nas[2]-1),]

  # 5. cambiare indici
  names(my_table) = as.character(my_table[1,])

  # 6. eliminare prima riga
  my_table = my_table[-1,]

  # 7. cambiare nomi colonne aggiungendo numeri
  col_names = names(my_table) # creare vettore coi nomi
  for(i in seq_along(col_names)) {
    col_names[i] = paste(i,'. ', col_names[i])
  }

  # 8. far selezionare le colonne chiedendo i numeri
  col_names = c(col_names, '.')
  user.input = dlg_list(col_names, multiple = T, preselect= '.', title = 'Select sample, target, ct mean and ct sd')$res

  user.input = user.input[-length(user.input)]
  col2take = vector('numeric')
  for(i in seq_along(user.input)) {
    temp = strsplit(user.input[i], ' . ')
    col2take = c(col2take, as.numeric(temp[[1]][1]))
  }

  # 9. ridurre la tabella
  my_table = my_table[,col2take] # selezionare giuste colonne
  my_table = my_table[seq(1, dim(my_table)[1],2),]

  # 10. trasformare in dataset
  my_data = data.frame(my_table)

  # 11. creo vettore dei geni e vettore condizioni
  genes = unique(my_data[,2])
  conditions = unique(my_data[,1])

  # creo vettore segno ± e creo vettore nomi colonne per tabelle
  my_sign = rep('±', length = length(conditions))
  new_names = c(names(my_data)[1], names(my_data)[3], 'Sign',names(my_data)[4])


  # 12. creo tante tabelle, una corrispondente per ogni gene ed esporto come xlsx
  for(i in seq_along(genes)){
    temp_str1 = paste(genes[i], ' <- my_data[my_data[,2] == genes[i],-2]')
    temp_str2 = paste(genes[i], '<- data.frame(', genes[i], '[,1], ', genes[i], '[,2], my_sign[1:dim(',
                      genes[i], ')[1]], ', genes[i], '[,3])')
    temp_str3 = paste('names(', genes[i], ') <- new_names')
    temp_str4 = paste(genes[i], '[,2] <- as.numeric(', genes[i], '[,2])')
    temp_str5 = paste(genes[i], '[,4] <- as.numeric(', genes[i], '[,4])')
    temp_str6 = paste('write.xlsx(',genes[i], ', file = "', name_end,'_mod.xlsx" ,sheetName="', genes[i], '", append=TRUE, row.name=FALSE)')


    eval(parse(text=temp_str1))
    eval(parse(text=temp_str2))
    eval(parse(text=temp_str3))
    eval(parse(text=temp_str4))
    eval(parse(text=temp_str5))
    eval(parse(text=temp_str6))
  }

}

desc_table <- function(data, group.var = NULL, cont.var = NULL, cat.var = NULL, paired = F){

  # librerie
  suppressPackageStartupMessages(library(Rmisc))
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(svGUI))
  suppressPackageStartupMessages(library(svDialogs))

  # Controlli qualità
  # Data è un dataframe?
  if (!is.data.frame(data)) {
    stop("Data MUST be a dataframe or a tibble")
  }

  # Creo tabella vuota
  tabella <- tibble(a = character(), b = character(), c = character(), d = numeric(), e = character())

  # Variabile gruppi
  # vedo se c'è, altrimenti chiedo
  if (is_empty(group.var)) {
    group.var <- dlg_list(colnames(data), multiple = F, title = "Select grouping variable")$res
  }

  # vedo se ne coniene solo 2
  if (length( unique( unlist( data[, group.var]))) != 2) {
    stop("Grouping variable MUST have 2 different values")
  }

  # ricavo numerosità gruppi
  numerosita <- data %>%
    group_by_at(vars(one_of(group.var))) %>%
    summarize(n = n())

  # Creo nomi colonne tabella e li cambio
  group1_name <- paste(numerosita[1, 1], ' (n = ', as.character(numerosita[1, 2]), ')', sep = '')
  group2_name <- paste(numerosita[2, 1], ' (n = ', as.character(numerosita[2, 2]), ')', sep = '')
  nomi_colonne <- c("Variabile", group1_name, group2_name, "p-value", "method")
  colnames(tabella) <- nomi_colonne

  # Chiedo quali sono le continue
  if (is_empty(cont.var)) {
    cont.var <- dlg_list(c(colnames(data), "preselect"), multiple = T, title = "Select continue variable/s, if none are present, press 0",
                         preselect = "preselect")$res
    cont.var <- cont.var[-length(cont.var)]
  }

  # Chiedo quali sono le categoriche
  if (is_empty(cat.var)) {
    cat.var <- dlg_list(c(colnames(data), "preselect"), multiple = T, title = "Select categorical variable/s, if none are present, press 0",
                        preselect = "preselect")$res
    cat.var <- cat.var[-length(cat.var)]
  }


  # statistica sulle continue
  if (!is_empty(cont.var)) {
    for (i in seq_along(cont.var)) {

      # nome variabile
      nome_var <- cont.var[i]

      # tabella summary e valori ctrl e intervention
      summary_tab <- summarySE(data = data, measurevar = cont.var[i], groupvars = group.var)
      group1_val <- paste(as.character( round( summary_tab[1, cont.var[i]] , 3)), '±',
                          as.character( round( summary_tab[1, "sd"] , 3)))

      group2_val <- paste(as.character( round( summary_tab[2, cont.var[i]] , 3)), '±',
                          as.character( round( summary_tab[2, "sd"] , 3)))

      # controllo distribuzione normale
      shap <- tapply(unlist(data[, cont.var[i]]), unlist(data[, group.var]), shapiro.test)

      # vedo se fare bartlett
      if (shap[[1]]$p.value > 0.05 & shap[[2]]$p.value > 0.05) {
        bart <- bartlett.test(unlist(data[, cont.var[i]]), unlist(data[, group.var]))

        # vedo se fare parametrico t.test o welch
        if (bart$p.val > 0.05) {
          p_val  <- t.test(unlist(data[, cont.var[i]]) ~ unlist(data[, group.var]), var.equal = T, paired = paired)$p.val
          method <- t.test(unlist(data[, cont.var[i]]) ~ unlist(data[, group.var]), var.equal = T, paired = paired)$method

        } else {
          p_val  <- t.test(unlist(data[, cont.var[i]]) ~ unlist(data[, group.var]), var.equal = F, paired = paired)$p.val
          method <- t.test(unlist(data[, cont.var[i]]) ~ unlist(data[, group.var]), var.equal = F, paired = paired)$method
        }
      } else {

        # se non parametrico
        p_val  <- suppressWarnings(wilcox.test(unlist(data[,cont.var[i]]) ~ unlist(data[, group.var]), paired = paired)$p.val)
        method <- suppressWarnings(wilcox.test(unlist(data[,cont.var[i]]) ~ unlist(data[, group.var]), paired = paired)$method)

      }

      # creo riga da aggiungere
      riga_da_agg <- data.frame(nome_var, group1_val, group2_val, p_val, method)
      colnames(riga_da_agg) <- nomi_colonne

      # aggiungo riga
      tabella <- rbind(tabella, riga_da_agg)
    }
  }

  # statistica sulle categoriche
  if (!is_empty(cat.var)) {
    for (i in seq_along(cat.var)) {
      # creo tabelle contingenza
      cont_tab <- table(unlist(data[, cat.var[i]]), unlist(data[, group.var]))
      perc1 <- round(cont_tab[, 1]/sum(cont_tab[, 1])*100, 3)
      perc2 <- round(cont_tab[, 2]/sum(cont_tab[, 2])*100, 3)
      perc_tab <- cbind(perc1, perc2)

      # creo vettori con freq (%)
      group1_cat <- paste(cont_tab[, 1], ' (', perc_tab[, 1], ')', sep="")
      names(group1_cat) <- as.character(1:length(rownames(cont_tab)))

      group2_cat <- paste(cont_tab[, 2], ' (', perc_tab[, 2], ')', sep="")
      names(group2_cat) <- as.character(1:length(rownames(cont_tab)))


      # creo dataframe dai due vettori sopra per poi aggiungerli nella tabella generale
      tab_cat <- data.frame(var = rownames(cont_tab), Group1 = group1_cat, Group2 = group2_cat, p = NA, meth = "")
      colnames(tab_cat) <- nomi_colonne

      # vedo che test usare
      # condizioni
      num_tot <- sum(cont_tab)
      gr_cate <- nrow(cont_tab)
      val_att <- suppressWarnings( chisq.test( unlist( data[, cat.var[i]]), unlist( data[, group.var]))$expected)

      # calcolo
      if (gr_cate > 2 | num_tot > 100) {
        p_cat  <- suppressWarnings(chisq.test(cont_tab)$p.val)
        method_cat <- suppressWarnings(chisq.test(cont_tab)$method)
      } else {
        if (all(val_att >5)) {
          p_cat  <- suppressWarnings(chisq.test(cont_tab, correct = T)$p.val)
          method_cat <- suppressWarnings(chisq.test(cont_tab, correct = T)$method)
        } else {
          p_cat  <- suppressWarnings(fisher.test(cont_tab)$p.val)
          method_cat <- suppressWarnings(fisher.test(cont_tab)$method)
        }
      }

      # creo vettore solo nome variabile, p.val e method
      varname <- data.frame(cat.var[i], "", "", p_cat, method_cat)
      colnames(varname) <- nomi_colonne

      # unisco le tabelle
      tabella <- rbind(tabella, varname)
      tabella <- rbind(tabella, tab_cat)

    }
  }
  rownames(tabella) <- NULL

  # tolgo spazi nome method
  tabella$method <- trimws(tabella$method, "both")

  return(tabella)
}

desc_kable <- function(data, group.var = NULL, cont.var = NULL, cat.var = NULL, paired = F){

  # librerie
  suppressPackageStartupMessages(library(Rmisc))
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(svGUI))
  suppressPackageStartupMessages(library(svDialogs))
  suppressPackageStartupMessages(library(kableExtra))
  suppressPackageStartupMessages(library(knitr))

  # Controlli qualità
  # Data è un dataframe?
  if (!is.data.frame(data)) {
    stop("Data MUST be a dataframe or a tibble")
  }

  # Creo tabella vuota
  tabella <- tibble(a = character(), b = character(), c = character(), d = numeric(), e = character())

  # Variable gruppi
  # vedo se c'è, altrimenti chiedo
  if (is_empty(group.var)) {
    group.var <- dlg_list(colnames(data), multiple = F, title = "Select grouping Variable")$res
  }

  # vedo se ne coniene solo 2
  if (length( unique( unlist( data[, group.var]))) != 2) {
    stop("Grouping Variable MUST have 2 different values")
  }

  # ricavo numerosità gruppi
  numerosita <- data %>%
    group_by_at(vars(one_of(group.var))) %>%
    summarize(n = n())

  # Creo nomi colonne tabella e li cambio
  group1_name <- paste(numerosita[1, 1], ' (n = ', as.character(numerosita[1, 2]), ')', sep = '')
  group2_name <- paste(numerosita[2, 1], ' (n = ', as.character(numerosita[2, 2]), ')', sep = '')
  nomi_colonne <- c("Variable", group1_name, group2_name, "p-value", "method")
  colnames(tabella) <- nomi_colonne

  # Chiedo quali sono le continue
  if (is_empty(cont.var)) {
    cont.var <- dlg_list(c(colnames(data), "preselect"), multiple = T, title = "Select continue Variable/s, if none are present, press 0",
                         preselect = "preselect")$res
    cont.var <- cont.var[-length(cont.var)]
  }

  # Chiedo quali sono le categoriche
  if (is_empty(cat.var)) {
    cat.var <- dlg_list(c(colnames(data), "preselect"), multiple = T, title = "Select categorical Variable/s, if none are present, press 0",
                        preselect = "preselect")$res
    cat.var <- cat.var[-length(cat.var)]
  }


  # statistica sulle continue
  if (!is_empty(cont.var)) {
    for (i in seq_along(cont.var)) {

      # nome Variable
      nome_var <- paste("<b>", cont.var[i],"</b>", sep = "")

      # tabella summary e valori ctrl e intervention
      summary_tab <- summarySE(data = data, measurevar = cont.var[i], groupvars = group.var)
      group1_val <- paste(as.character( round( summary_tab[1, cont.var[i]] , 3)), '±',
                          as.character( round( summary_tab[1, "sd"] , 3)))

      group2_val <- paste(as.character( round( summary_tab[2, cont.var[i]] , 3)), '±',
                          as.character( round( summary_tab[2, "sd"] , 3)))

      # controllo distribuzione normale
      shap <- tapply(unlist(data[, cont.var[i]]), unlist(data[, group.var]), shapiro.test)

      # vedo se fare bartlett
      if (shap[[1]]$p.value > 0.05 & shap[[2]]$p.value > 0.05) {
        bart <- bartlett.test(unlist(data[, cont.var[i]]), unlist(data[, group.var]))

        # vedo se fare parametrico t.test o welch
        if (bart$p.val > 0.05) {
          p_val  <- t.test(unlist(data[, cont.var[i]]) ~ unlist(data[, group.var]), var.equal = T, paired = paired)$p.val
          method <- t.test(unlist(data[, cont.var[i]]) ~ unlist(data[, group.var]), var.equal = T, paired = paired)$method

        } else {
          p_val  <- t.test(unlist(data[, cont.var[i]]) ~ unlist(data[, group.var]), var.equal = F, paired = paired)$p.val
          method <- t.test(unlist(data[, cont.var[i]]) ~ unlist(data[, group.var]), var.equal = F, paired = paired)$method
        }
      } else {

        # se non parametrico
        p_val  <- suppressWarnings(wilcox.test(unlist(data[,cont.var[i]]) ~ unlist(data[, group.var]), paired = paired)$p.val)
        method <- suppressWarnings(wilcox.test(unlist(data[,cont.var[i]]) ~ unlist(data[, group.var]), paired = paired)$method)

      }

      # creo riga da aggiungere
      riga_da_agg <- c(nome_var, group1_val, group2_val, (p_val), method)
      names(riga_da_agg) <- nomi_colonne

      # aggiungo riga
      tabella <- rbind(tabella, t(data.frame(riga_da_agg)))
    }
  }

  # statistica sulle categoriche
  if (!is_empty(cat.var)) {
    for (i in seq_along(cat.var)) {
      # creo tabelle contingenza
      cont_tab <- table(unlist(data[, cat.var[i]]), unlist(data[, group.var]))
      perc1 <- round(cont_tab[, 1]/sum(cont_tab[, 1])*100, 3)
      perc2 <- round(cont_tab[, 2]/sum(cont_tab[, 2])*100, 3)
      perc_tab <- cbind(perc1, perc2)

      # creo vettori con freq (%)
      group1_cat <- paste(cont_tab[, 1], ' (', perc_tab[, 1], ')', sep="")
      names(group1_cat) <- as.character(1:length(rownames(cont_tab)))

      group2_cat <- paste(cont_tab[, 2], ' (', perc_tab[, 2], ')', sep="")
      names(group2_cat) <- as.character(1:length(rownames(cont_tab)))


      # creo dataframe dai vettori sopra per poi aggiungerli nella tabella generale
      tab_cat <- data.frame(var = rownames(cont_tab), Group1 = group1_cat, Group2 = group2_cat, p = rep("", length(rownames(cont_tab))), meth = NA)
      colnames(tab_cat) <- nomi_colonne

      # vedo che test usare
      # condizioni
      num_tot <- sum(cont_tab)
      gr_cate <- nrow(cont_tab)
      val_att <- suppressWarnings( chisq.test( unlist( data[, cat.var[i]]), unlist( data[, group.var]))$expected)

      # calcolo
      if (gr_cate > 2 | num_tot > 100) {
        p_cat  <- suppressWarnings(chisq.test(cont_tab)$p.val)
        method_cat <- suppressWarnings(chisq.test(cont_tab)$method)
      } else {
        if (all(val_att >5)) {
          p_cat  <- suppressWarnings(chisq.test(cont_tab, correct = T)$p.val)
          method_cat <- suppressWarnings(chisq.test(cont_tab, correct = T)$method)
        } else {
          p_cat  <- suppressWarnings(fisher.test(cont_tab)$p.val)
          method_cat <- suppressWarnings(fisher.test(cont_tab)$method)
        }
      }

      # creo vettore solo nome Variable, p.val e method
      varname <- data.frame(paste("<b>", cat.var[i],"</b>", sep = ""), "", "", p_cat, method_cat)
      colnames(varname) <- nomi_colonne

      # unisco le tabelle
      tabella <- rbind(tabella, varname)
      tabella <- rbind(tabella, tab_cat)

      # tolgo eventuali spazi nei nomi method
      tabella$method <- trimws(tabella$method, which = "both")

    }
  }
  rownames(tabella) <- NULL

  # tolgo spazi nome method
  tabella$method <- trimws(tabella$method, "both")

  # aggiungo simboli

    # trovo metodi unici
    un_meth <- unique(na.omit(tabella$method))

    # uso ciclo per aggiungere footnote marker alfabetico

    for (i in seq_along(un_meth)) {
      tabella$Variable[which(tabella$method == un_meth[i])] <- paste(tabella$Variable[which(tabella$method == un_meth[i])], footnote_marker_alphabet(i))
    }

    footnote_met <- paste(unique(na.omit(tabella$method)), ";", sep = "")

  # tolgo method
  tabella <- tabella[,-5]

  # modifico colonna p-value
  tabella <- tabella %>%
    mutate(`p-value` = case_when(
      .$`p-value` == ""  ~ "",
      as.numeric(.$`p-value`) <0.001 ~ "<0.001",
      as.numeric(.$`p-value`) <0.01 ~ "<0.01",
      TRUE ~ as.character(round(as.numeric(.$`p-value`), 3))
    ))

  # creo kable html

  text1 <- paste('tabella_html <- kbl(tabella, escape = F)%>%
    kable_classic(full_width = F)%>%
    footnote(general = "Statistical tests",
             alphabet = footnote_met)%>%
     add_header_above(c(" " = 1,', group.var ,' = 2, " " = 1))')

  eval(parse(text = text1))

  return(tabella_html)
}
