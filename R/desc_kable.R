#' @name desc_kable
#' @title Descriptive statistics kable
#' @description This function returns a descriptive statistics kable for both categorical and continues variable for a 2-group dataset.
#' It calculates also the p-value choosing the right statistical analysis.
#' @usage desc_kable(data, group.var = NULL, cont.var = NULL, cat.var = NULL, paired = F)
#' @param data Dataframe to use
#' @param group.var String of the grouping variable column name. The grouping variable MUST have only 2 different values
#' @param cont.var Chr string (or vector) of the continuous variable/s column name to take into account
#' @param cat.var Chr string (or vector) of the categorical variable/s column name to take into account
#' @param paired Logical, are continuous variables paired?
#' @returns A kable with the following columns: variable name - group_1 values (mean ± sd for continuous, n(\%) for categorical) - group_2 values (mean ± sd for continuous, n(\%) for categorical) - p-value
#' Moreover, statistical tests are summarized as footnote
#' @seealso \code{\link{desc_table}} for descriptive table
#' @author Matteo Miotto
#'
#' @importFrom Rmisc summarySE
#' @importFrom tibble tibble
#' @importFrom rlang is_empty
#' @importFrom svDialogs dlg_list
#' @import dplyr
#' @importFrom kableExtra kable_classic footnote add_header_above footnote_marker_alphabet kbl
#' @export



desc_kable <- function(data, group.var = NULL, cont.var = NULL, cat.var = NULL, paired = F){

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
    summarize(n = n()) %>%
    as.data.frame(stringAsFactors = F)

  # Creo nomi colonne tabella e li cambio
  group1_name <- paste(numerosita[1, 1], ' (n = ', as.character(numerosita[1, 2]), ')', sep = '')
  group2_name <- paste(numerosita[2, 1], ' (n = ', as.character(numerosita[2, 2]), ')', sep = '')
  nomi_colonne <- c("Variable", group1_name, group2_name, "p-value", "method")
  colnames(tabella) <- nomi_colonne

  # Chiedo quali sono le continue
  if (is.null(cont.var)) {
    cont.var <- dlg_list(c(colnames(data), "preselect"), multiple = T, title = "Select continue variable/s, if none are present, press 0",
                         preselect = "preselect")$res
    cont.var <- cont.var[-length(cont.var)]
    if (length(cont.var) == 0) {cont.var <- NA}
  }

  suppressWarnings(if (is.na(cont.var)) {cont.var <- NULL})

  # Chiedo quali sono le categoriche
  if (is.null(cat.var)) {
    cat.var <- dlg_list(c(colnames(data), "preselect"), multiple = T, title = "Select categorical variable/s, if none are present, press 0",
                        preselect = "preselect")$res
    cat.var <- cat.var[-length(cat.var)]
    if (length(cat.var) == 0) {cat.var <- NA}
  }

  suppressWarnings(if (is.na(cat.var)) {cat.var <- NULL})


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
      perc1 <- round(cont_tab[, 1]/sum(cont_tab[, 1])*100, 2)
      perc2 <- round(cont_tab[, 2]/sum(cont_tab[, 2])*100, 2)
      perc_tab <- cbind(perc1, perc2)

      # creo vettori con freq (%)
      group1_cat <- paste(cont_tab[, 1], ' (', perc_tab[, 1], '%)', sep="")
      names(group1_cat) <- as.character(1:length(rownames(cont_tab)))

      group2_cat <- paste(cont_tab[, 2], ' (', perc_tab[, 2], '%)', sep="")
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
