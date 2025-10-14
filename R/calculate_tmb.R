#' @title Calculate tumour mutation burden.
#'
#' @description This function implements tumor mutation burden calculation.
#' TODO: add more details.
#'
#' @param maf_data Incoming data frame representing simple somatic
#' mutations in maf format. The minimal required columns are
#' Tumor_Sample_Barcode, Hugo_Symbol, Chromosome, Start_Position, End_Position,
#' Variant_Classification.
#' @param regions_bed Optionally, specify the data frame with sequenced panel.
#' Expected to be in bed format with first 3 columns indicating chromosome,
#' start, and end positions. This function is agnostic to the naming of these
#' columns. If there are any other columns present, they are ignored.
#' @param projection The genome build of projection for the data in the incoming
#' maf. Default is "grch37".
#' @param subset_to_nonSyn Logical argument to control whether the maf file
#' should be subset to non-synonymous mutations only. Defaults to TRUE (perform
#' subsetting).
#' @param vc_nonSyn Vector of annotations from Variant_Classification column
#' that should be considered non-synonymous.
#'
#' @return data frame
#'
#' @examples
#' # obtain maf data
#' library(GAMBLR.open)
#' maf1 <- get_ssm_by_samples()
#'
#' coding_counts = calculate_tmb(maf1)
#' 
#' dplyr::filter(coding_counts,Tumor_Sample_Barcode=="OCI-Ly10")
#' 
#' all_counts = calculate_tmb(
#'     maf1,
#'     subset_to_nonSyn = FALSE
#' )
#' dplyr::filter(all_counts,Tumor_Sample_Barcode=="OCI-Ly10")
#'
#' @import dplyr GAMBLR.data
#' @export
#'
calculate_tmb <- function(
    maf_data,
    regions_bed,
    projection = "grch37",
    subset_to_nonSyn = TRUE,
    vc_nonSyn = NULL
){

    if(missing(maf_data)){
        stop(
            "Please provide maf data to operate on."
        )
    }

    # Handle possibility of duplicates
    maf_data <- maf_data %>%
        dplyr::distinct(
            Tumor_Sample_Barcode,
            Hugo_Symbol,
            Chromosome,
            Start_Position,
            End_Position,
            Variant_Classification,
            .keep_all = TRUE
        )

    # get all samples in maf
    sample_ids <- maf_data %>%
        dplyr::distinct(Tumor_Sample_Barcode)

    if(subset_to_nonSyn){
        # Define non-synonymous variants
        if(is.null(vc_nonSyn)){
            vc_nonSynonymous <- vc_nonSynonymous
        }else{
            vc_nonSynonymous <- vc_nonSyn
        }

        # Subset maf to only non-synonymous variants
        maf_data <- maf_data %>%
            filter(
                Variant_Classification %in% vc_nonSynonymous
            )
    }

    # assume WGS for the TSB denominator
    # Will be reset below if custom bed panel is provided
    if(projection == "grch37"){
        chromosomes <- GAMBLR.data::chromosome_arms_grch37
    }else{
        chromosomes <- GAMBLR.data::chromosome_arms_hg38
    }

    denominator_size <- chromosomes %>%
        dplyr::group_by(chromosome) %>%
        dplyr::summarize(length = max(end)) %>%
        ungroup %>%
        dplyr::summarize(sum(length)) %>%
        pull

    # Optionally subset to panel
    if(!missing(regions_bed)){
        # Expect bed format but be flexible about column names
        columns <- colnames(regions_bed)[1:3]

        overlap <- cool_overlaps(
            data1 = maf_data,
            data2 = regions_bed,
            columns2 = columns
        )

        denominator_size <- regions_bed %>%
            dplyr::mutate(
                length = base::get(columns[3]) - base::get(columns[2])
            ) %>%
            dplyr::summarize(sum(length)) %>%
            pull
    }

    # Now calculate N of mutations and TMB
    tmb <- maf_data %>%
        dplyr::count(Tumor_Sample_Barcode) %>%
        dplyr::mutate(total_perMB = 1000000 * n/denominator_size) %>%
        select(
            Tumor_Sample_Barcode,
            total = n,
            total_perMB
        )

    tmb <- left_join(
        sample_ids,
        tmb
    ) %>%
        replace(is.na(.), 0)  %>%
        dplyr::arrange(total)

    return(tmb)
}
