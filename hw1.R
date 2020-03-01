#!/usr/bin/env Rscript

### Usage:      Rscript --vanilla hw1.R <input file> <score file>
### Example:    Rscript --vanilla hw1.R input.txt blosum62.txt
### Note:       Smith-Waterman Algorithm

### Rscript --vanilla ~/Documents/Yale_Med/PhD/Gersteins_class/HW1/hw1.R ~/Documents/Yale_Med/PhD/Gersteins_class/HW1/input.txt ~/Documents/Yale_Med/PhD/Gersteins_class/HW1/blosum62.txt
### This is one way to read in arguments in R We need to read input file and score file.
args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("At least two arguments must be supplied (inputFile, scoreFile).n", call.=FALSE)
} else if (length(args)>=2) {
  # default gap penalties for opening and extension gap, respectively
  args[3] = -2
  args[4] = -1
}

## Specifying author and email
p <- c(person("Maryam", "Zekavat", role = "aut", email = "maryam.zekavat@yale.edu"))

## Implementing Smith-Waterman Algorithm
runSW <- function(inputFile, scoreFile, openGap = -2, extGap = -1) {
  
  ### Loading sequence File:
  inFile = read.table(inputFile,header=F, as.is=T, stringsAsFactors=F,sep = "", comment.char = '',quote = "\"") 
  print("Aligning Sequences:")
  print(inFile)
  print(paste("openGap:", openGap, "| extGap", extGap))
  
  ### Initializing a scoring matrix data frame (newDF) where the first row referances the Sequence 1 letters and the first column refers to Sequence 2. 
  Seq_1 = inFile[1,1]
  Seq_2 = inFile[2,1]
  
  Seq_1_strings = unlist(strsplit(Seq_1, ""))
  Seq_2_strings = unlist(strsplit(Seq_2, ""))
 
  Num_Rows = length(Seq_2_strings)+2 # Note: No colnames or rownames are used, which explains the need for "+2" here
  Num_Cols = length(Seq_1_strings)+2
  N=Num_Rows*Num_Cols
  
  newDF = data.frame(matrix(vector(mode='numeric', length=N), nrow=Num_Rows, ncol=Num_Cols),check.names=FALSE)
  newDF[1,] = c("","", Seq_1_strings)
  newDF[c(2:Num_Rows),1]= c("",Seq_2_strings)
  
  ### Making similarity matrix using the scoreFile
  Scores_DF = read.table(scoreFile,header=T, as.is=T, stringsAsFactors=F,sep = "", comment.char = '',quote = "\"") 
  print("Making similarity matrix using scoreFile...")
  for (i in 3:Num_Cols){
    Seq_1_letter = newDF[1,i]
    for (j in 3:Num_Rows){
      Seq_2_letter = newDF[j,1]
      Score = Scores_DF[which(colnames(Scores_DF) == Seq_1_letter), which(rownames(Scores_DF) == Seq_2_letter)]
      newDF[j,i] = as.numeric(Score)
    }
  }
  
  ### iterating through Smith-Waterman Algorithm and making a new dataframe called SW_DF:
  SW_DF = newDF
  print("Running Smith-Waterman...")

  openGap=as.numeric(openGap)
  extGap = as.numeric(extGap)

  for (i in 3:Num_Rows){
    Colgap_vector = rep(0, (i-2))
    for (cols in 1:(i-2)){Colgap_vector[cols] = openGap + extGap*(cols-1)}
    
    for (j in 3:Num_Cols){
      Rowgap_vector = rep(0, j-2)
      for (cols in 1:(j-2)){Rowgap_vector[cols] = openGap + extGap*(cols-1)}
      
      #For column neighbors, finding the current cell value by taking the max of the following elements: 
      InitialScore = as.numeric(SW_DF[i,j]) + as.numeric(SW_DF[i-1,j-1])
      Max_colneighbor = max(rev(as.numeric(SW_DF[c(2:(i-1)),j])) + Colgap_vector) 
      Max_rowneighbor = max(rev(as.numeric(SW_DF[i,c(2:(j-1))])) + Rowgap_vector) 
      cur_cellVal = max(InitialScore, Max_colneighbor, Max_rowneighbor, 0)
      SW_DF[i,j] = cur_cellVal
      
    }
  }
  
  ### now finding the Alignment Score, the max of the whole matrix:
  MAX_Matrix = max(data.matrix(SW_DF[c(3:Num_Rows),c(3:Num_Cols)]))
  print(paste("Alignment Score:", MAX_Matrix))
  SW_DF_withNA = SW_DF
  SW_DF_withNA[1,] = NA
  SW_DF_withNA[,1] = NA
  SW_Matrix = data.matrix(SW_DF_withNA)
  
  MAX_indx = which(SW_Matrix == MAX_Matrix, arr.ind = T)

  ### Finding the matrix element correspondng to the Max Alignment Score:
  row_i = MAX_indx[1]
  col_i = MAX_indx[2]

  ### Initializing the alignment data frame (AlignmentDF):
  #Will first setup the AlignmentDF, then at the very end will reverse each column and transpose the DF such that it is in the desired human-readable format (3 rows).
  #The AlignmentDF contains 3 fields: Seq_1_alig, Aligned, Seq_2_alig.
  #The Seq_1_alig and Seq_2_alig fields refer to the first and second sequence alignments, respectively 
  #The Aligned vector field refers to the "|" denoting sequence match.

  AlignmentDF = data.frame(Seq_1_alig =c(), Aligned = c(), Seq_2_alig = c())
  
  ### First, Dealing with the non-aligning ends after the max alignment:
  #Will need to paste the flanking non-aligned sequences into their respective seq_alignment vectors and add enough spaces in the other vector such that the alignments line up.
  #Find which free non-aligning end is longer:
    Seq1_nonaligningLength = Num_Cols - (col_i)
    Seq2_nonaligningLength = Num_Rows - (row_i)
    
    #If sequence 1 right-sided non-aligning end is longer, then:
    if(Seq1_nonaligningLength > Seq2_nonaligningLength){
      if (Seq2_nonaligningLength > 0) {
        cur_Seq_1_alig = c(Seq_1_strings[(Num_Cols-2):(col_i - 1)], c(")"))  #the -2 are because of the two additional rows included in the matrices
        
        #accounting additional spaces needed in seq2 for the human-readable alignment to work: 
        Spaces = rep(" ", Seq1_nonaligningLength - Seq2_nonaligningLength)
        
        cur_Seq_2_alig = c(Seq_2_strings[(NumRows-2):(row_i-1)], Spaces, c(")") )
        cur_Seq_2_alig = cur_Seq_2_alig[which(!is.na(cur_Seq_2_alig))]
        curAligned = rep(" ", Seq1_nonaligningLength+1)
        
        AlignmentDF = rbind(AlignmentDF, data.frame(Seq_1_alig = cur_Seq_1_alig, Aligned = curAligned, Seq_2_alig = cur_Seq_2_alig))
      }
      else{ #if Seq2_nonaligningLength ==0
        cur_Seq_1_alig = c(Seq_1_strings[(Num_Cols-2):(col_i - 1)], c(")"))  #the -2 are because of the two additional rows included in the matrices
        
        #accounting additional spaces needed in seq2 for the human-readable alignment to work: 
        Spaces = rep(" ", Seq1_nonaligningLength - Seq2_nonaligningLength)
        
        cur_Seq_2_alig = c(Spaces, c(")") )
        cur_Seq_2_alig = cur_Seq_2_alig[which(!is.na(cur_Seq_2_alig))]
        curAligned = rep(" ", Seq1_nonaligningLength+1)
        
        AlignmentDF = rbind(AlignmentDF, data.frame(Seq_1_alig = cur_Seq_1_alig, Aligned = curAligned, Seq_2_alig = cur_Seq_2_alig))
      }
    
    } 
    if(Seq1_nonaligningLength < Seq2_nonaligningLength){
      if (Seq1_nonaligningLength > 0) {
        cur_Seq_2_alig = c(Seq_2_strings[(Num_Rows-2):(row_i - 1)], c(")")) #the -2 are because of the two additional cols included in the matrices
        
        #accounting additional spaces needed in seq2 for the human-readable alignment to work: 
        Spaces = rep(" ", Seq2_nonaligningLength - Seq1_nonaligningLength)
        
        cur_Seq_1_alig = c(Seq_1_strings[(NumCols-2):(col_i-1)], Spaces, c(")") )
        cur_Seq_1_alig = cur_Seq_1_alig[which(!is.na(cur_Seq_1_alig))]
        curAligned = rep(" ", Seq2_nonaligningLength+1) #the plus one accounts for the ")"
        
        AlignmentDF = rbind(AlignmentDF, data.frame(Seq_1_alig = cur_Seq_1_alig, Aligned = curAligned, Seq_2_alig = cur_Seq_2_alig))
      }
      else{ #if Seq1_nonaligningLength ==0
        cur_Seq_2_alig = c(Seq_2_strings[(Num_Rows-2):(row_i - 1)], c(")")) #the -2 are because of the two additional cols included in the matrices
        
        #accounting additional spaces needed in seq2 for the human-readable alignment to work: 
        Spaces = rep(" ", Seq2_nonaligningLength - Seq1_nonaligningLength)
        
        cur_Seq_1_alig = c(Spaces, c(")") )
        cur_Seq_1_alig = cur_Seq_1_alig[which(!is.na(cur_Seq_1_alig))]
        curAligned = rep(" ", Seq2_nonaligningLength+1) #the plus one accounts for the ")"
        
        AlignmentDF = rbind(AlignmentDF, data.frame(Seq_1_alig = cur_Seq_1_alig, Aligned = curAligned, Seq_2_alig = cur_Seq_2_alig))
      }
    } 
    if(Seq1_nonaligningLength == Seq2_nonaligningLength){
      if(Seq1_nonaligningLength >0){ #if the two ends are equal in length but longer than 0:
          cur_Seq_2_alig = c(Seq_2_strings[(Num_Rows-2):(row_i - 1)], c(")") )#the -2 are because of the two additional cols included in the matrices
          cur_Seq_1_alig = c(Seq_1_strings[(NumCols-2):(col_i-1)], c(")"))
          cur_Seq_1_alig = cur_Seq_1_alig[which(!is.na(cur_Seq_1_alig))]
          curAligned = rep(" ", Seq1_nonaligningLength+1) #the plus one accounts for the ")"
          
          AlignmentDF = rbind(AlignmentDF, data.frame(Seq_1_alig = cur_Seq_1_alig, Aligned = curAligned, Seq_2_alig = cur_Seq_2_alig))
      }
      else{ #if the two ends ==0
          cur_Seq_2_alig = c(")")  #the -1 are because of the two additional cols included in the matrices
          cur_Seq_1_alig = c(")") 
          curAligned = rep(" ", Seq1_nonaligningLength+1) #the plus one accounts for the ")"
          
          AlignmentDF = rbind(AlignmentDF, data.frame(Seq_1_alig = cur_Seq_1_alig, Aligned = curAligned, Seq_2_alig = cur_Seq_2_alig))
      }
    } 

### Adding the first best alignment:  
row_i = MAX_indx[1]
col_i = MAX_indx[2]
cur_Seq_1_alig = SW_DF[1, col_i]
cur_Seq_2_alig = SW_DF[row_i, 1]
curAligned = ifelse(cur_Seq_1_alig == cur_Seq_2_alig, "|", " ")
AlignmentDF = rbind(AlignmentDF, data.frame(Seq_1_alig = cur_Seq_1_alig, Aligned = curAligned, Seq_2_alig = cur_Seq_2_alig))


# Finding the Best Sequence Alignment:
# Begin moving up, left, or diagonally within the matrix till we hit a zero
# Move in the direction of the max_cell. If in case of multliple max_cells (ie: ties), choose left first, then diagonal, then up, to maximize the aligned sequence length. 
# Continuously adding to the alignment dataframe each time.  
# For moves up and left there will be an insertion/deletion, so adding "-" to the seq1 with up mves and "-" to seq2 with left moves
 orig_AlignmentDF = AlignmentDF
 row_i = MAX_indx[1]
 col_i = MAX_indx[2]

 while (SW_Matrix[(row_i-1), (col_i-1)] != 0){
    # Values of the current cell, the cells up, left, diagonally up, and the max
    up       = SW_Matrix[(row_i - 1), col_i]
    left     = SW_Matrix[row_i, (col_i - 1)]
    diagn    = SW_Matrix[(row_i - 1), (col_i - 1)]
    max_cell = max(up, left, diagn)
    
    if (left == max_cell) {
      col_i = col_i - 1
      cur_Seq_1_alig = SW_DF[1, col_i]
      cur_Seq_2_alig = "-"
      curAligned = ifelse(cur_Seq_1_alig == cur_Seq_2_alig, "|", " ")
      AlignmentDF = rbind(AlignmentDF, data.frame(Seq_1_alig = cur_Seq_1_alig, Aligned = curAligned, Seq_2_alig = cur_Seq_2_alig))
    } else if (diagn == max_cell) {
      row_i = row_i - 1
      col_i = col_i - 1
      cur_Seq_1_alig = SW_DF[1, col_i]
      cur_Seq_2_alig = SW_DF[row_i, 1]
      curAligned = ifelse(cur_Seq_1_alig == cur_Seq_2_alig, "|", " ")
      AlignmentDF = rbind(AlignmentDF, data.frame(Seq_1_alig = cur_Seq_1_alig, Aligned = curAligned, Seq_2_alig = cur_Seq_2_alig))
    } else if (up == max_cell) {
      row_i = row_i - 1
      cur_Seq_1_alig = "-"
      cur_Seq_2_alig = SW_DF[row_i, 1]
      curAligned = ifelse(cur_Seq_1_alig == cur_Seq_2_alig, "|", " ")
      AlignmentDF = rbind(AlignmentDF, data.frame(Seq_1_alig = cur_Seq_1_alig, Aligned = curAligned, Seq_2_alig = cur_Seq_2_alig))
    } 
    
  }
  
 midAlignment=AlignmentDF
 mid_col_i = col_i
 mid_row_i = row_i
 
 ### Now adding the remaining free end at the beginning
 Seq1_nonaligningLength = col_i-3
 Seq2_nonaligningLength = row_i-3
 
 #If sequence 1 right-sided non-aligning end is longer, then:
 if(Seq1_nonaligningLength > Seq2_nonaligningLength){
   if (Seq2_nonaligningLength > 0) {
     cur_Seq_1_alig = c(c("("), Seq_1_strings[Seq1_nonaligningLength:1]) 
     
     #accounting additional spaces needed in seq2 for the human-readable alignment to work: 
     Spaces = rep(" ", Seq1_nonaligningLength - Seq2_nonaligningLength)
     
     cur_Seq_2_alig = c(c("("), Seq_2_strings[Seq2_nonaligningLength:1], Spaces )
     cur_Seq_2_alig = cur_Seq_2_alig[which(!is.na(cur_Seq_2_alig))]
     curAligned = rep(" ", Seq1_nonaligningLength+1)
     
     AlignmentDF = rbind(AlignmentDF, data.frame(Seq_1_alig = cur_Seq_1_alig, Aligned = curAligned, Seq_2_alig = cur_Seq_2_alig))
   }
   else{ #if Seq2_nonaligningLength ==0
     cur_Seq_1_alig = c(c("("), Seq_1_strings[Seq1_nonaligningLength:1]) 
     
     #accounting additional spaces needed in seq2 for the human-readable alignment to work: 
     Spaces = rep(" ", Seq1_nonaligningLength - Seq2_nonaligningLength)
     
     cur_Seq_2_alig = c(c("("),Spaces )
     cur_Seq_2_alig = cur_Seq_2_alig[which(!is.na(cur_Seq_2_alig))]
     curAligned = rep(" ", Seq1_nonaligningLength+1)
     
     AlignmentDF = rbind(AlignmentDF, data.frame(Seq_1_alig = cur_Seq_1_alig, Aligned = curAligned, Seq_2_alig = cur_Seq_2_alig))
   }
 }
 if(Seq1_nonaligningLength < Seq2_nonaligningLength){
     if (Seq1_nonaligningLength > 0) {
       cur_Seq_2_alig = c(c("("), Seq_2_strings[Seq2_nonaligningLength:1]) 
       
       #accounting additional spaces needed in seq2 for the human-readable alignment to work: 
       Spaces = rep(" ", Seq2_nonaligningLength - Seq1_nonaligningLength)
       
       cur_Seq_1_alig = c(c("("), Seq_1_strings[Seq1_nonaligningLength:1], Spaces )
       cur_Seq_1_alig = cur_Seq_1_alig[which(!is.na(cur_Seq_1_alig))]
       curAligned = rep(" ", Seq2_nonaligningLength+1) #the plus one accounts for the ")"
       
       AlignmentDF = rbind(AlignmentDF, data.frame(Seq_1_alig = cur_Seq_1_alig, Aligned = curAligned, Seq_2_alig = cur_Seq_2_alig))
     }
     else{ #if Seq1_nonaligningLength ==0
       cur_Seq_2_alig = c(c("("), Seq_2_strings[Seq2_nonaligningLength:1]) 
       
       #accounting additional spaces needed in seq2 for the human-readable alignment to work: 
       Spaces = rep(" ", Seq2_nonaligningLength - Seq1_nonaligningLength)
       
       cur_Seq_1_alig = c(c("("), Spaces )
       cur_Seq_1_alig = cur_Seq_1_alig[which(!is.na(cur_Seq_1_alig))]
       curAligned = rep(" ", Seq2_nonaligningLength+1) #the plus one accounts for the ")"
       
       AlignmentDF = rbind(AlignmentDF, data.frame(Seq_1_alig = cur_Seq_1_alig, Aligned = curAligned, Seq_2_alig = cur_Seq_2_alig))
     }
   } 
 if(Seq1_nonaligningLength == Seq2_nonaligningLength){
     if(Seq1_nonaligningLength >0){ #if the two ends are equal in length but longer than 0:
       cur_Seq_2_alig = c(c("(") , Seq_2_strings[Seq2_nonaligningLength:1])
       cur_Seq_1_alig = c(c("(") , Seq_1_strings[Seq1_nonaligningLength:1])
       cur_Seq_1_alig = cur_Seq_1_alig[which(!is.na(cur_Seq_1_alig))]
       curAligned = rep(" ", Seq1_nonaligningLength+1) #the plus one accounts for the ")"
       
       AlignmentDF = rbind(AlignmentDF, data.frame(Seq_1_alig = cur_Seq_1_alig, Aligned = curAligned, Seq_2_alig = cur_Seq_2_alig))
     }
     else{ #if the two ends ==0
       cur_Seq_2_alig = c("(")  #the -1 are because of the two additional cols included in the matrices
       cur_Seq_1_alig = c("(") 
       curAligned = rep(" ", Seq1_nonaligningLength+1) #the plus one accounts for the ")"
       
       AlignmentDF = rbind(AlignmentDF, data.frame(Seq_1_alig = cur_Seq_1_alig, Aligned = curAligned, Seq_2_alig = cur_Seq_2_alig))
     }
   } 
   
 ### Reversing each field within AlignmentDF such that it is in the original order of the input sequences, and transposing it such that we have 3 rows instead of 3 columns.
 Final_AlignmentDF = data.frame(Seq_1_alig = rev(AlignmentDF$Seq_1_alig), Aligned = rev(AlignmentDF$Aligned), Seq_2_alig = rev(AlignmentDF$Seq_2_alig))
 tFinal_AlignmentDF = t(Final_AlignmentDF)
 
 ### Writing otuput files
 outfile_components = unlist(strsplit(inputFile, "/"))
 
 prefixFileName = paste(outfile_components[1:(length(outfile_components)-1)], collapse="/")
 outfile_components_Scores = paste(prefixFileName, "/AlignmentScores_", MAX_Matrix,".output_scoreFile.tsv", sep="")
 outfile_components_Alignment= paste(prefixFileName, "/FinalAlignment_", MAX_Matrix,".output_Alignment.txt", sep="")
 
 write.table(SW_DF,outfile_components_Scores,
             col.names = F, row.names = F, quote = F, sep = "\t")
 
 write.table(tFinal_AlignmentDF,outfile_components_Alignment,
             col.names = F, row.names = F, quote = F, sep = "")
 
 print("Final Alignment:")
 print(paste(tFinal_AlignmentDF[1,], collapse=""))
 print(paste(tFinal_AlignmentDF[2,], collapse=""))
 print(paste(tFinal_AlignmentDF[3,], collapse=""))
 
 print(paste("Score Matrix Output:", outfile_components_Scores))
 print(paste("Final Alignment Output:", outfile_components_Alignment))
 print("Done.")
}

## Run the main function and generate results
runSW(inputFile=args[1], scoreFile=args[2], openGap=args[3], extGap=args[4])