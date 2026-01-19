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


    #------------ Combine DevSeq core ortholog expression tables with GOslim data -------------
    
    # Show message
    message("Starting analysis...")









}



