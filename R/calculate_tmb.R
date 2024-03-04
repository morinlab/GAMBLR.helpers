#' @title Return TMB counts.
#'
#' @description This function implements tumor mutation burden calculation.
#' TODO: add more details.
#'
#' @param maf_data Incoming data frame representing simple somatic
#' mutations in maf format.
#' @param vc_nonSyn Vector of annotations from Variant_Classification column
#' that should be considered non-synonymous.
#'
#'
#'
#' @return data frame
#'
#' @examples
#' # obtain maf data
#' maf1 <- get_coding_ssm(
#'     these_sample_ids = "DOHH-2"
#' )
#'
#' TODO
#'
#' @import dplyr
#' @export
#'
calculate_tmb <- function(
    maf_data,
    regions_bed,
    this_seq_type = "genome",
    projection = "grch37",
    subset_to_nonSyn = TRUE,
    verbose = FALSE,
    these_sample_ids = NULL,
    these_samples_metadata = NULL,
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

    # get samples with the dedicated helper function
    metadata <- id_ease(
        these_samples_metadata = these_samples_metadata,
        these_sample_ids = these_sample_ids,
        verbose = verbose,
        this_seq_type = this_seq_type
    )

    sample_ids <- metadata$sample_id

    maf_data <- maf_data %>%
        dplyr::filter(Tumor_Sample_Barcode %in% sample_ids)

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

        overlap <- GAMBLR.helpers::cool_overlaps(
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
        dplyr::mutate(total_perMB = n/denominator_size) %>%
        select(
            Tumor_Sample_Barcode,
            total = n,
            total_perMB)

    tmb <- left_join(
        metadata %>%
            dplyr::select(Tumor_Sample_Barcode),
        tmb
    ) %>%
        replace(is.na(.), 0)  %>%
        dplyr::arrange(total)

    return(tmb)
}
