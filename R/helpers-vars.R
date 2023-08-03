#adding coding_vc to global enviroment
#' @export
coding_vc = c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Region", "Splice_Site", "Targeted_Region", "Translation_Start_Site")

#Global variable specifying what metadata columns are absolutely required
#' @export
required_cols = c("sample_id","patient_id","pathology","seq_type","genome_build","pairing_status","Tumor_Sample_Barcode")

#' @export
cnames = c("CHROM_A", "START_A", "END_A", "CHROM_B", "START_B", "END_B", "NAME", "SOMATIC_SCORE", "STRAND_A", "STRAND_B", "TYPE", "FILTER", "VAF_tumour", "VAF_normal", "DP_tumour", "DP_normal", "tumour_sample_id", "normal_sample_id", "pair_status")

#' @export
coding_class = c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Silent", "Splice_Region", "Splice_Site", "Targeted_Region", "Translation_Start_Site")

#' @export
maf_header = c("Hugo_Symbol"=1,"Entrez_Gene_Id"=2,"Center"=3,"NCBI_Build"=4,"Chromosome"=5,"Start_Position"=6,"End_Position"=7,"Strand"=8,"Variant_Classification"=9,"Variant_Type"=10,"Reference_Allele"=11,"Tumor_Seq_Allele1"=12,"Tumor_Seq_Allele2"=13,"dbSNP_RS"=14,"dbSNP_Val_Status"=15,"Tumor_Sample_Barcode"=16,"Matched_Norm_Sample_Barcode"=17,"Match_Norm_Seq_Allele1"=18,"Match_Norm_Seq_Allele2"=19,"Tumor_Validation_Allele1"=20,"Tumor_Validation_Allele2"=21,"Match_Norm_Validation_Allele1"=22,"Match_Norm_Validation_Allele2"=23,"Verification_Status"=24,"Validation_Status"=25,"Mutation_Status"=26,"Sequencing_Phase"=27,"Sequence_Source"=28,"Validation_Method"=29,"Score"=30,"BAM_File"=31,"Sequencer"=32,"Tumor_Sample_UUID"=33,"Matched_Norm_Sample_UUID"=34,"HGVSc"=35,"HGVSp"=36,"HGVSp_Short"=37,"Transcript_ID"=38,"Exon_Number"=39,"t_depth"=40,"t_ref_count"=41,"t_alt_count"=42,"n_depth"=43,"n_ref_count"=44,"n_alt_count"=45,"all_effects"=46,"Allele"=47,"Gene"=48,"Feature"=49,"Feature_type"=50,"Consequence"=51,"cDNA_position"=52,"CDS_position"=53,"Protein_position"=54,"Amino_acids"=55,"Codons"=56,"Existing_variation"=57,"ALLELE_NUM"=58,"DISTANCE"=59,"STRAND_VEP"=60,"SYMBOL"=61,"SYMBOL_SOURCE"=62,"HGNC_ID"=63,"BIOTYPE"=64,"CANONICAL"=65,"CCDS"=66,"ENSP"=67,"SWISSPROT"=68,"TREMBL"=69,"UNIPARC"=70,"RefSeq"=71,"SIFT"=72,"PolyPhen"=73,"EXON"=74,"INTRON"=75,"DOMAINS"=76,"AF"=77,"AFR_AF"=78,"AMR_AF"=79,"ASN_AF"=80,"EAS_AF"=81,"EUR_AF"=82,"SAS_AF"=83,"AA_AF"=84,"EA_AF"=85,"CLIN_SIG"=86,"SOMATIC"=87,"PUBMED"=88,"MOTIF_NAME"=89,"MOTIF_POS"=90,"HIGH_INF_POS"=91,"MOTIF_SCORE_CHANGE"=92,"IMPACT"=93,"PICK"=94,"VARIANT_CLASS"=95,"TSL"=96,"HGVS_OFFSET"=97,"PHENO"=98,"MINIMISED"=99,"GENE_PHENO"=100,"FILTER"=101,"flanking_bps"=102,"vcf_id"=103,"vcf_qual"=104,"gnomAD_AF"=105,"gnomAD_AFR_AF"=106,"gnomAD_AMR_AF"=107,"gnomAD_ASJ_AF"=108,"gnomAD_EAS_AF"=109,"gnomAD_FIN_AF"=110,"gnomAD_NFE_AF"=111,"gnomAD_OTH_AF"=112,"gnomAD_SAS_AF"=113,"vcf_pos"=114,"gnomADg_AF"=115,"blacklist_count"=116)

#' @export
colour_aliases = list("COO_consensus" = "coo", "COO" = "coo", "DHITsig_consensus" = "coo",
                      "pathology" = "pathology", "analysis_cohort" = "pathology", "group" = "pathology",
                      "FL_group" = "FL", "lymphgen" = "lymphgen", "lymphgen_with_cnv" = "lymphgen",
                      "bcl2_ba" = "pos_neg", "BCL2_status" = "pos_neg", "myc_ba" = "pos_neg",
                      "bcl6_ba" = "pos_neg", "EBV_status_inf"="EBV_status",
                      "manta_BCL2_sv" = "pos_neg", "manual_BCL2_sv" = "pos_neg", "manta_MYC_sv" = "pos_neg")

#' @export
rainfall_conv = c("T>C", "T>C", "C>T", "C>T", "T>A", "T>A", "T>G", "T>G", "C>A", "C>A", "C>G", "C>G", "InDel")
names(rainfall_conv) = c('A>G', 'T>C', 'C>T', 'G>A', 'A>T', 'T>A', 'A>C', 'T>G', 'C>A', 'G>T', 'C>G', 'G>C', 'InDel')
