options(warn = -1)
library(Peptides)
library(peptider)
library(tools)
library(log4r, warn.conflicts = FALSE)

exit <- function() 
{ 
    invokeRestart("abort") 
}

tryCatch(
    expr = {
        my_logfile = "logs.txt"
        my_console_appender = console_appender(layout = default_log_layout())
        my_file_appender = file_appender(my_logfile, append = TRUE, layout = default_log_layout())
        my_logger <- log4r::logger(threshold = "INFO", appenders= list(my_console_appender, my_file_appender))
        
        input <- commandArgs(trailingOnly = TRUE)
        
        if(length(input) != 1)
        {
            log4r::error(my_logger, "Incorrect number of parameters passed! Allowed number of parameters: 1")
            exit()
        }
        
        if(!file.exists(input[1]))
        {
            log4r::error(my_logger, "File does not exist!")
            exit()
        }
        
        if (file_ext(input[1]) != "csv")
        {
            log4r::error(my_logger, "The file passed has wrong file extension! Allowed file extension: .csv")
            exit()
        }
        
        df <- read.csv(input, header = TRUE)
        
        if (ncol(df) != 2)
        {
            log4r::error(my_logger, "Input file must contain two columns only!")
            exit()
        }
        
        i <- 0
        for(seq in df$Peptide.Sequence)
        {
            i <- i + 1

            f1 <- lengthpep(seq)
            f2 <- aIndex(seq)
            f3 <- boman(seq)
            f4 <- hmoment(seq)
            f5 <- instaIndex(seq)
            f6 <- ppeptide(seq, libscheme = "NNK", N = 10^8)
            f7 <- getNofNeighbors(seq)
            f8 <- hydrophobicity(seq)
            f9 <- mw(seq)
            f10 <- pI(seq)
            f11 <- charge(seq)
            f12 <- codons(seq, libscheme = "NNK", flag = FALSE)
            f13 <- makowski(f1, "NNK")
            
            out <- data.frame(f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13)
            
            f14 <- unlist(aaComp(seq))
            for(j in f14)
            {
                out <- cbind(out, j)
            }
            
            f15 <- unlist(kideraFactors(seq))
            for(j in f15)
            {
                out <- cbind(out, j)
            }
            
            out <- cbind(df$Peptide.Sequence[i], out, df$Target[i])
            colnames(out) <- c("Sequence", "Length", "Aliphatic Index", "Boman Index", "Hydrophobic Moment", "Instability Index", "Probability of Detection", "No of Neighbors", "Hydrophobicity Index", "Molecular Weight", "Isoelectric Point", "Theoretical Net Charge", "No of Codon Representations", "Makowski Diversity Index", "Tiny", "Small", "Aliphatic", "Aromatic", "Non Polar", "Polar", "Charged", "Basic", "Acidic", "Mole% Tiny", "Mole% Small", "Mole% Aliphatic", "Mole% Aromatic", "Mole% Non Polar", "Mole% Polar", "Mole% Charged", "Mole% Basic", "Mole% Acidic", "KF1", "KF2", "KF3", "KF4", "KF5", "KF6", "KF7", "KF8", "KF9", "KF10", "Target")
            if(i == 1)
            {
                write.table(out, "output.csv", sep = ",", row.names = F, col.names = T, append = T)
            }else
            {
                write.table(out, "output.csv", sep = ",", row.names = F, col.names = F, append = T)
            }
        }
    }
) 