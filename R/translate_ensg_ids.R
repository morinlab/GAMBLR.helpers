#' @title Translate ENSG ids.
#'
#' @description Translate gene ids in ENSG formats (with or without trailing
#'      .<number>) to a human-readable format of Hugo symbols. Note that there
#'      are more ENSG ids then there are hugo symbols: 20,075 of ENSG ids don't
#'      have associated Hugo symbol. Those are novel transcripts, lncRNAS,
#'      pseudogenes etc, and they will be dropped from the resulting output of
#'      this function.
#'
#'
#' @param df an incoming data frame with ENSG ids either in column or rownames.
#' @param ids_in_rownames boolean specifying whether ENSG ids are in rownames or
#'      in a designated column. Default is FALSE (the ENSG ids are a column).
#' @param ensg_column_in string specifying column name where to find the ENSG
#'      ids. Only used when `ids_in_rownames` is set to FALSE. Default is
#'      "gene_id".
#' @param symbol_column_out string specifying column name where to output the
#'      ENSG ids. Only used when `ids_in_rownames` is set to FALSE. Default is
#'      "gene_id".
#'
#' @return df
#'
#' @import dplyr tibble GAMBLR.data
#' @export
#'
#' @examples
#' library(tibble)
#' test_this <- GAMBLR.data::gencode_to_symbol %>%
#'     select(1) %>%
#'     head(10) %>%
#'     mutate(test = "test")
#'
#' translate_ensg_ids(
#'     df = test_this %>%
#'         column_to_rownames("ensembl_gene_id"),
#'     ids_in_rownames = TRUE
#' )
#' 
#' translate_ensg_ids(
#'     df = test_this
#' )
#' 
#' translate_ensg_ids(
#'     df = test_this %>%
#'         select(
#'             test_col = ensembl_gene_id,
#'             everything()
#'         ),
#'     ensg_column_in = "test_col"
#' )
#' 
#' translate_ensg_ids(
#'     df = test_this %>%
#'         select(
#'             test_col = ensembl_gene_id,
#'             everything()
#'         ),
#'     ensg_column_in = "test_col",
#'     symbol_column_out = "this, Ha!"
#' )
#' 
#' 
translate_ensg_ids = function(
    df,
    ensg_column_in = "ensembl_gene_id",
    symbol_column_out = "hgnc_symbol",
    ids_in_rownames = FALSE
){
    if(ids_in_rownames){
        # When ENSG is in rownames, convert them to a column
        df <- df %>%
            as.data.frame %>%
            rownames_to_column("ensembl_gene_id")
    }else{
        # Otherwise, rename the column for consistency
        df <- df %>%
            rename(
                "ensembl_gene_id" := !!ensg_column_in
            )

    }

    # Remove PAR_Y from the converter
    converter <- GAMBLR.data::gencode_to_symbol %>%
        filter(!grepl("PAR_Y", gene_id)) %>%
        select(ensembl_gene_id, hgnc_symbol)
    
    # Handle possible .<number> in the incoming data
    df <- df %>%
        mutate(
            ensembl_gene_id = sub("\\..*", "", ensembl_gene_id)
        )
    
    # Make the translation
    df <- left_join(
        df,
        converter
    ) %>%
        select(-ensembl_gene_id)
    
    # Format return df
    if(ids_in_rownames){
        df <- df %>%
            column_to_rownames("hgnc_symbol")
    }else{
        df <- df %>%
            select(hgnc_symbol, everything()) %>%
            rename(
                !!symbol_column_out := "hgnc_symbol"
            )
    }
    return(df)
}
