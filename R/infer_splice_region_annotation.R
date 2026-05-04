#' @title Infer Splice Region Annotation
#'
#' @description Extracts cDNA position and splice offset (e.g., c.123+1G>A) to calculate the 
#' affected amino acid and assign a protein label (e.g., p.Splice_Donor+1@AA41).
#'
#' @details The function uses ceiling(abs(cdna_pos)/3) to determine the Amino Acid position.
#' Labels are assigned as Splice_Donor (+), Splice_Acceptor (-), or NearSplice (no offset).
#'
#' @param maf_data A dataframe containing HGVSc and HGVSp_Short columns.
#'
#' @return A dataframe with updated HGVSp_Short labels based on splice region.
#'
#' @import dplyr
#' @import stringr
#' @export
#'
#' @examples
#' #' #this is an example
#' ###For your reference, this line is exactly 100 characters. Do not exceed 100 characters per line

function_name = function(a_parameter,
                         another_parameter){
}

library(dplyr)
library(stringr)

input <- GAMBLR.open::get_coding_ssm()
maf <- input

infer_splice_region_annotation <- function(
    maf_data
){
    # Fix the splice region annotations
    maf <- maf_data %>%
    mutate(
      # Extract base cDNA coordinate
      cdna_pos = as.numeric(str_extract(HGVSc, "(?<=c\\.)-?\\d+")),
      
      # Extract splice offset (+N or -N)
      splice_offset = str_extract(HGVSc, "(?<=\\d)[+-]\\d+"),
      
      # Infer AA position (absolute for UTR cases)
      Amino_Acid_Position = ceiling(abs(cdna_pos) / 3),
      
      # Create protein label from splice direction
      HGVSp_Short = case_when(
        
        !is.na(HGVSp_Short) ~ HGVSp_Short,
        
        !is.na(splice_offset) & str_starts(splice_offset, "\\+") ~
          paste0("p.Splice_Donor", splice_offset, "@AA", Amino_Acid_Position),
        
        !is.na(splice_offset) & str_starts(splice_offset, "\\-") ~
          paste0("p.Splice_Acceptor", splice_offset, "@AA", Amino_Acid_Position),
        
        TRUE ~
          paste0("p.NearSplice@AA", Amino_Acid_Position)
      )
    ) %>%
    select(-cdna_pos, -splice_offset, -Amino_Acid_Position)
    
  return(maf) 
}

infer_splice_region_annotation(maf_data = input)
input %>%
  select(HGVSp_Short) %>% 
  head(100) %>%
  as.data.frame()
  # will display all the values included ones N/A

nonpopulated <- input %>% 
  as.data.frame() %>%
  filter(is.na(HGVSp_Short))
  # filter(): looks at the condition, sees if specific row matches the condition, if true will use, if not, discard
  # is.na: function that checks whether value is empty or not 

outcome <- infer_splice_region_annotation(maf_data = nonpopulated)
outcome %>%
  select(HGVSp_Short) %>%
  head(100) %>%
  as.data.frame()

populated <- input %>%
  as.data.frame() %>%
  filter(!is.na(HGVSp_Short))

output <- infer_splice_region_annotation(maf_data = populated)
output %>%
  select(HGVSp_Short) %>%
  head(100) %>%
  as.data.frame()

identical(populated, output)
all.equal(
  nonpopulated,  
  outcome
)
# identical() and all.equal() is same, but latter will provide more verbose output

