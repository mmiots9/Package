#' @name pcr.lab.th
#' @title Create a PCR result file suitable for thesis lab
#' @description This function returns an Excel file in which each sheet represent one gene, and in each sheet the table have the following columns: condition, ctmean, ± sign, ct SD
#' @usage pcr.lab.th()
#' @returns Excel file
#' @author Matteo Miotto
#' @importFrom readxl read_excel
#' @importFrom svDialogs dlg_list
#' @importFrom xlsx write.xlsx

#' @export




pcr.lab.th = function(){
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
