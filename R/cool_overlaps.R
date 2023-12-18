cool_overlaps <- function(
    data1,
    data2,
    columns1 = c("Chromosome", "Start_Position", "End_Position"),
    columns2 = c("Chromosome", "Start_Position", "End_Position"),
    type = "any"
){

    # Ensure all columns provided for overlap are present in the data frame
    if(! length(columns1) == length(intersect(columns1, colnames(data1)))){
        stop(
            "Not all of the requested columns for overlap in data1 are present."
        )
    }

    if(! length(columns2) == length(intersect(columns2, colnames(data2)))){
        stop(
            "Not all of the requested columns for overlap in data2 are present."
        )
    }

    # What is the name of the column in columns1 that specifies start and end?
    start1 <- columns1[grepl("start", columns1, ignore.case = TRUE)]
    end1 <- columns1[grepl("end", columns1, ignore.case = TRUE)]

    # What is the name of the column in columns1 that specifies start and end?
    start2 <- columns1[grepl("start", columns2, ignore.case = TRUE)]
    end2 <- columns1[grepl("end", columns2, ignore.case = TRUE)]

    # Prepare for overlap
    overlap <- dplyr::inner_join(
        data1,
        data2,
        by = structure(names = columns1, .Data = columns2),
        relationship = "many-to-many"
    )

    # Return matches based on mode
    if(type == "any"){
        message(
            "Running in default mode of any..."
        )
        overlap <- overlap %>%
            dplyr::filter(
               !!sym(start2) <= !!sym(end1)
            )
    } else if (type == "start"){
        message(
            "Running in the mode start..."
        )
        overlap <- overlap %>%
            dplyr::filter(
               !!sym(start1) == !!sym(start2)
            )
    } else if (type == "end"){
        message(
            "Running in the mode end..."
        )
        overlap <- overlap %>%
            dplyr::filter(
               !!sym(end1) == !!sym(end2)
            )
    } else if (type == "within"){
        message(
            "Running in the mode within..."
        )
        overlap <- overlap %>%
            dplyr::filter(
               (!!sym(start1) >= !!sym(start2)) & (!!sym(end1) <= !!sym(end2))
            )
    } else if (type == "equal"){
        message(
            "Running in the mode equal..."
        )
        overlap <- overlap %>%
            dplyr::filter(
               (!!sym(start1) == !!sym(start2)) & (!!sym(end1) == !!sym(end2))
            )
    } else {
        message(
            "You have requested mode that is not supported."
        )
        stop(
            "Please supply one of any, start, end, within, or equal with type."
        )
    }

    return(overlap)
}
