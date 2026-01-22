# Open Agilent Tapestation RIN data from RNA extracted from plant meristem tissue sections;
# Sections were embedded in paraffin, deparaffinized with Xylene, and dehydrated using EtOH


plotRIN <- function(temp = c("56", "62"), ...) {


    # Show error message if no temperature is given
    if ((missing(temp)) || (!is.element(temp, c("56", "62"))))

        stop ("Please choose temp value of either '56' or '62'",
            call. = TRUE
            )


    # Set file path for input files and read table
    if (temp == "56") {

        tpth = file.path(in_dir, "Tapestation", "Fixation and embedding test_paraffin_samples_56.csv")
        ts_data <- read.table(tpth, sep = ";", dec = ".", header = TRUE, stringsAsFactors = FALSE)

    } else if (temp == "62") {

        tpth = file.path(in_dir, "Tapestation", "Fixation and embedding test_paraffin_samples_62.csv")
        ts_data <- read.table(tpth, sep = ";", dec = ".", header = TRUE, stringsAsFactors = FALSE)
    }


    # return_list <- list("ts_data" = ts_data, "temp" = temp)
    # return(return_list)
    # }
    # return_objects <- plotRIN(temp = "62")
    # list2env(return_objects, envir = .GlobalEnv)

    # Show message
    message("Reading data...")


    #------------ Process data: get average and sd for RNA concentration and RIN values ------------
    
    # Show message
    message("Starting analysis...")


    # Get average from replicates
    
    if (temp == "62") {
        numb <- 2

    } else if (temp == "56") {
        numb <- 3
    }


    getRepl <- function(x) { 

        split(x, 
            rep(1:(nrow(x)/numb), 
                each = numb)
            )
    }

    repl_lst <- getRepl(ts_data)


    plt_df <- do.call(rbind, lapply(repl_lst, function(x) {

        rna_mean <- sum(x[,3])/numb
        rna_sd <- sd(x[,3])

        rin_mean <- sum(x[,4])/numb
        rin_sd <- sd(x[,4])

        x_df <- data.frame(sample_type = gsub('.{2}$', '', x[2,2]), RNA_avg = rna_mean, 
                           RNA_sd = rna_sd, RIN_avg = rin_mean, RIN_sd = rin_sd)
        return(x_df)

        }))



    # Plot data
    plotData <- function(df) {

        fname <- sprintf('%s.pdf', paste(deparse(substitute(df)), temp , sep="_"))

        if (temp == 62) {

            plot_mar <- c(0.5, 0.5, 1, 1)
            plt_title <- "GO enrichment global q1"
            df$sample_type <- gsub("Col-0_contro", "Col-0_snap", df$sample_type)
            df$sample_type <- gsub("EtOH_70%", "70% EtOH", df$sample_type)
            df$sample_type <- gsub("EtOH_100%", "100% EtOH", df$sample_type)
            df$sample_type <- gsub("Aceton_100%", "100% Aceton", df$sample_type)
            df$sample_type <- gsub("EtOH_Acetic Acid_3:1", "EtOH_AA_3:1", df$sample_type)

        } else if (temp == 56) {

            plot_mar <- c(0.5, 0.5, 1, 1)
            plt_title <- "GO enrichment global q1"

        }


        p <- ggplot(df, aes(x = RIN_avg, y = RNA_avg)) +
        geom_point() + 
        scale_x_continuous(expand = c(0, 0), limits = c(0, 10)) + 
        labs(x = "RIN", y = "RNA concentration (ng)") + 
        theme(panel.background = element_blank(), 
            axis.ticks.length = unit(0.25, "cm"), 
            axis.ticks = element_line(colour = "black", size = 1.0), 
            axis.line = element_line(colour = 'black', size = 1.0), 
            plot.margin = unit(plot_mar, "cm"), 
            plot.title = element_text(size = 22.75, margin = margin(t = 0, r = 0, b = 9, l = 0), hjust = 0.5),
            axis.title.y = element_text(size = 22.75, margin = margin(t = 0, r = 7.0, b = 0, l = 10), 
                colour = "black", face = "bold"), 
            axis.title.x = element_text(size = 22.75, margin = margin(t = 0.5, r = 0, b = 8.15, l = 0), 
                colour = "black", face = "bold"), 
            axis.text.x = element_text(size = 18.8, margin = margin(t = 3.5, b = 7), colour = "grey20"), 
            axis.text.y = element_text(size = 19.0, angle = 0, margin = margin(l = 10, r = -2), colour = "black")
        )

        ggsave(file = file.path(out_dir, "plots", fname), plot = p, 
               width = 11.5, height = 6.5, dpi = 300, units = c("in"), limitsize = FALSE)


    }


    plotData(plt_df)





}



