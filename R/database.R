cnames = c("CHROM_A", "START_A", "END_A", "CHROM_B", "START_B", "END_B", "NAME", "SOMATIC_SCORE", "STRAND_A", "STRAND_B", "TYPE", "FILTER", "VAF_tumour", "VAF_normal", "DP_tumour", "DP_normal", "tumour_sample_id", "normal_sample_id", "pair_status")



coding_class = c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Silent", "Splice_Region", "Splice_Site", "Targeted_Region", "Translation_Start_Site")



maf_header = c("Hugo_Symbol"=1,"Entrez_Gene_Id"=2,"Center"=3,"NCBI_Build"=4,"Chromosome"=5,"Start_Position"=6,"End_Position"=7,"Strand"=8,"Variant_Classification"=9,"Variant_Type"=10,"Reference_Allele"=11,"Tumor_Seq_Allele1"=12,"Tumor_Seq_Allele2"=13,"dbSNP_RS"=14,"dbSNP_Val_Status"=15,"Tumor_Sample_Barcode"=16,"Matched_Norm_Sample_Barcode"=17,"Match_Norm_Seq_Allele1"=18,"Match_Norm_Seq_Allele2"=19,"Tumor_Validation_Allele1"=20,"Tumor_Validation_Allele2"=21,"Match_Norm_Validation_Allele1"=22,"Match_Norm_Validation_Allele2"=23,"Verification_Status"=24,"Validation_Status"=25,"Mutation_Status"=26,"Sequencing_Phase"=27,"Sequence_Source"=28,"Validation_Method"=29,"Score"=30,"BAM_File"=31,"Sequencer"=32,"Tumor_Sample_UUID"=33,"Matched_Norm_Sample_UUID"=34,"HGVSc"=35,"HGVSp"=36,"HGVSp_Short"=37,"Transcript_ID"=38,"Exon_Number"=39,"t_depth"=40,"t_ref_count"=41,"t_alt_count"=42,"n_depth"=43,"n_ref_count"=44,"n_alt_count"=45,"all_effects"=46,"Allele"=47,"Gene"=48,"Feature"=49,"Feature_type"=50,"Consequence"=51,"cDNA_position"=52,"CDS_position"=53,"Protein_position"=54,"Amino_acids"=55,"Codons"=56,"Existing_variation"=57,"ALLELE_NUM"=58,"DISTANCE"=59,"STRAND_VEP"=60,"SYMBOL"=61,"SYMBOL_SOURCE"=62,"HGNC_ID"=63,"BIOTYPE"=64,"CANONICAL"=65,"CCDS"=66,"ENSP"=67,"SWISSPROT"=68,"TREMBL"=69,"UNIPARC"=70,"RefSeq"=71,"SIFT"=72,"PolyPhen"=73,"EXON"=74,"INTRON"=75,"DOMAINS"=76,"AF"=77,"AFR_AF"=78,"AMR_AF"=79,"ASN_AF"=80,"EAS_AF"=81,"EUR_AF"=82,"SAS_AF"=83,"AA_AF"=84,"EA_AF"=85,"CLIN_SIG"=86,"SOMATIC"=87,"PUBMED"=88,"MOTIF_NAME"=89,"MOTIF_POS"=90,"HIGH_INF_POS"=91,"MOTIF_SCORE_CHANGE"=92,"IMPACT"=93,"PICK"=94,"VARIANT_CLASS"=95,"TSL"=96,"HGVS_OFFSET"=97,"PHENO"=98,"MINIMISED"=99,"GENE_PHENO"=100,"FILTER"=101,"flanking_bps"=102,"vcf_id"=103,"vcf_qual"=104,"gnomAD_AF"=105,"gnomAD_AFR_AF"=106,"gnomAD_AMR_AF"=107,"gnomAD_ASJ_AF"=108,"gnomAD_EAS_AF"=109,"gnomAD_FIN_AF"=110,"gnomAD_NFE_AF"=111,"gnomAD_OTH_AF"=112,"gnomAD_SAS_AF"=113,"vcf_pos"=114,"gnomADg_AF"=115,"blacklist_count"=116)



#' @title Add ICGC metadata.
#'
#' @description Layer on ICGC metadata from a supplemental table to fill in missing COO.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::get_gambl_metadata], not meant for out-of-package usage.
#'
#' @param incoming_metadata A metadata table (probably output from `get_gambl_metadata`).
#'
#' @return Metadata with layered information (ICGC).
#'
#' @import dplyr readr stringr
#'
#' @noRd
#'
#' @examples
#' my_meta = get_gambl_metadata()
#' icgc_metadata = add_icgc_metadata(incoming_metadata = my_meta)
#'
add_icgc_metadata = function(incoming_metadata){
  repo_base = check_config_value(config::get("repo_base"))
  icgc_publ_file = paste0(repo_base,"data/metadata/raw_metadata/MALY_DE_tableS1.csv")
  icgc_publ = suppressMessages(suppressWarnings(read_csv(icgc_publ_file)))
  icgc_publ = icgc_publ[,c(1:20)]
  #fix commas as decimals
  icgc_publ = mutate(icgc_publ, purity = str_replace(purity, ",", "."))
  icgc_publ = mutate(icgc_publ, sex = str_to_upper(sex))

  icgc_raw_path = paste0(repo_base,"data/metadata/raw_metadata/ICGC_MALY_seq_md.tsv")

  #check for missingness
  if(!file.exists(icgc_raw_path)){
    print(paste("missing: ", icgc_raw_path))
    message("Have you cloned the GAMBL repo and added the path to this directory under the local section of your config?")
  }

  icgc_raw = suppressMessages(read_tsv(icgc_raw_path))

  icgc_raw = icgc_raw %>%
    dplyr::select(-compression, -bam_available, -read_length, -time_point, -unix_group, -ffpe_or_frozen, -link_name)  %>%
    dplyr::filter(tissue_status %in% c("tumor", "tumour"))

  icgc_all = left_join(icgc_raw, icgc_publ,by = "ICGC_ID") %>%
    dplyr::select(-tissue_status, -seq_type, -protocol, -seq_source_type, -data_path, -genome_build, -RNA_available) %>%
    dplyr::select(sample_id, ICGC_ID, pathology.x, pathology.y, COO, molecular_BL, MYC_sv, BCL2_sv, BCL6_sv) %>%
    dplyr::rename("ICGC_MYC_sv" = "MYC_sv") %>%
    dplyr::rename("ICGC_BCL2_sv" = "BCL2_sv") %>%
    dplyr::rename("ICGC_BCL6_sv" = "BCL6_sv") %>%
    dplyr::rename("detailed_pathology" = "pathology.x") %>%
    dplyr::rename("ICGC_PATH" = "pathology.y")

  #join with all metadata to fill in blanks
  #all_meta=get_gambl_metadata()
  rejoined = left_join(incoming_metadata, icgc_all,by = "sample_id") %>%
    mutate(COO_final = case_when(!is.na(COO_consensus) ~ COO_consensus, COO != "n.a." & COO != "TypeIII" ~ COO, TRUE ~ "NA")) %>%
    dplyr::select(-COO)
  return(rejoined)
}



add_prps_result = function(incoming_metadata){
  prps_res = suppressMessages(read_tsv("/projects/rmorin/projects/gambl-repos/gambl-rmorin/results/icgc_dart/derived_and_curated_metadata/outputs/BL_dhitsig_PRPS.tsv"))
  colnames(prps_res)[1] = "sample_id"
  prps_res = dplyr::select(prps_res, sample_id, PRPS_score, PRPS_class)

  #need to associate each sample with a patient ID then annotate the metadata based on patient ID
  patient_meta_g = get_gambl_metadata(seq_type_filter = "genome") %>%
    dplyr::select(sample_id, patient_id)

  patient_meta_r = get_gambl_metadata(seq_type_filter = "mrna") %>%
    dplyr::select(sample_id, patient_id)

  patient_meta = bind_rows(patient_meta_g, patient_meta_r)
}



#' @title Append To Table.
#'
#' @description Housekeeping function to add results to a table.
#'
#' @details INTERNAL FUNCTION, not meant for out-of-package usage.
#'
#' @param table_name The name of the database table to update/populate.
#' @param data_df A dataframe of values to load into the table.
#'
#' @return A table.
#'
#' @import RMariaDB DBI
#'
#' @noRd
#'
#' @examples
#' table_up = append_to_table("my_table", "my_df")
#'
append_to_table = function(table_name,
                           data_df){

  db = check_config_value(config::get("database_name"))
  con = DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)
  dbWriteTable(con, table_name, table_data, append = TRUE)
}



#' @title Get GAMBL Outcomes.
#'
#' @description Get the patient-centric clinical metadata.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::get_gambl_metadata], not meant for out-of-package usage.
#'
#' @param patient_ids Vector of patient IDs.
#' @param time_unit Return follow-up times in one of three time units: year, month or day. Default is "year".
#' @param censor_cbioportal Optionally request the censoring to be encoded in the specific style required by cBioPortal. Default is FALSE.
#' @param complete_missing Optionally fill in any gaps to ensure we have values for every patient (censor at 0 if missing). Default is FALSE.
#' @param from_flatfile Optionally set to FALSE to use the database to get the survival data. Default is TRUE.
#'
#' @return Data frame with one row for each patient_id.
#'
#' @import tidyr dplyr readr RMariaDB DBI
#'
#' @noRd
#'
#' @examples
#' outcome_df = get_gambl_outcomes()
#'
get_gambl_outcomes = function(patient_ids,
                              time_unit = "year",
                              censor_cbioportal = FALSE,
                              complete_missing = FALSE,
                              from_flatfile = TRUE){

  if(from_flatfile){
    outcome_flatfile = paste0(check_config_value(config::get("repo_base")), check_config_value(config::get("table_flatfiles")$outcomes))

    #check for missingness
    if(!file.exists(outcome_flatfile)){
      print(paste("missing: ", outcome_flatfile))
      message("Have you cloned the GAMBL repo and added the path to this directory under the local section of your config?")
    }

    all_outcome = suppressMessages(read_tsv(outcome_flatfile))

  }else{
    db = check_config_value(config::get("database_name"))
    con = DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)
    all_outcome = dplyr::tbl(con, "outcome_metadata") %>%
      as.data.frame()

    DBI::dbDisconnect(con)
  }
  if(!missing(patient_ids)){
    all_outcome = all_outcome %>%
      dplyr::filter(patient_id %in% patient_ids)

    if(complete_missing){
      #add NA values and censored outcomes for all missing patient_ids
      all_outcome = all_outcome %>%
        complete(patient_id = patient_ids, fill = list(OS_YEARS = 0, PFS_years = 0, TTP_YEARS = 0, DSS_YEARS = 0, CODE_OS = 0, CODE_PFS = 0, CODE_DSS = 0, CODE_TTP = 0))
    }
  }
  if(time_unit == "month"){
    all_outcome = all_outcome %>%
      mutate(OS_MONTHS = OS_YEARS * 12)

    all_outcome = all_outcome %>%
      mutate(PFS_MONTHS = PFS_YEARS * 12)

    all_outcome = all_outcome %>%
      mutate(TTP_MONTHS = TTP_YEARS * 12)

    all_outcome = all_outcome %>%
      mutate(DSS_MONTHS = DSS_YEARS * 12)

    all_outcome = all_outcome %>%
      dplyr::select(-c("OS_YEARS", "PFS_YEARS", "TTP_YEARS", "DSS_YEARS"))

  }else if(time_unit == "day"){
    all_outcome = all_outcome %>%
      mutate(OS_DAYS = OS_YEARS * 365)

    all_outcome = all_outcome %>%
      mutate(PFS_DAYS = PFS_YEARS * 365)

    all_outcome = all_outcome %>%
      mutate(TTP_DAYS = TTP_YEARS * 365)

    all_outcome = all_outcome %>%
      mutate(DSS_DAYS = DSS_YEARS * 365)

    all_outcome = all_outcome %>%
      dplyr::select(-c("OS_YEARS", "PFS_YEARS", "TTP_YEARS", "DSS_YEARS"))
  }

  #if necessary, convert the censoring into the cBioPortal format for OS and PFS
  if(censor_cbioportal){
    all_outcome$OS_STATUS = as.character(all_outcome$CODE_OS)
    all_outcome = all_outcome %>%
      mutate(OS_STATUS = case_when(OS_STATUS == "0" ~ "0:LIVING", OS_STATUS == "1"~"1:DECEASED"))

    all_outcome$DFS_STATUS = as.character(all_outcome$CODE_PFS)
    all_outcome = all_outcome %>%
      mutate(DFS_STATUS = case_when(DFS_STATUS == "0" ~ "0:DiseaseFree", DFS_STATUS == "1"~"1:Recurred/Progressed"))

    all_outcome = all_outcome %>%
      mutate(all_outcome, DFS_MONTHS = PFS_MONTHS)
  }
  all_outcome = all_outcome %>%
    mutate(is_adult = ifelse(age < 20, "Pediatric", "Adult"))

  return(all_outcome)
}
