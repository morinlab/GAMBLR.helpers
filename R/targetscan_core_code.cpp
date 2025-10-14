#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;
// [[Rcpp::export]]
DataFrame targetscan_rcpp(DataFrame maf, DataFrame mirna){
  StringVector maf_chr = maf["Chromosome"];
  NumericVector maf_start = maf["Start_Position"];
  NumericVector maf_end = maf["End_Position"];
  StringVector mi_chr = mirna["Chromosome"];
  NumericVector mi_start = mirna["Start_Position"];
  NumericVector mi_end = mirna["End_Position"];
  StringVector mi_rna = mirna["miRNA"];
  NumericVector mi_site = mirna["sites"];
  StringVector slc_chr;
  NumericVector slc_start;
  NumericVector slc_end;
  StringVector slc_mirna;
  NumericVector slc_sites;
  for(int i=0; i<maf.nrows(); ++i){
    for(int j=0; j<mirna.nrows(); ++j){
      if (maf_chr[i] == mi_chr[j] && ((mi_start[j] <= maf_start[i] && mi_end[j] >= maf_start[i]) || (mi_start[j] >= maf_start[i] && mi_start[j] <= maf_end[i]))){
        slc_chr.push_back(maf_chr[i]);
        slc_start.push_back(maf_start[i]);
        slc_end.push_back(maf_end[i]);
        slc_mirna.push_back(mi_rna[j]);
        slc_sites.push_back(mi_site[j]);
      }
    }
  }
  DataFrame df = DataFrame::create( Named("Chromosome") = slc_chr,
                                    Named("Start_Position") = slc_start,
                                    Named("End_Position") = slc_end,
                                    Named("miRNA") = slc_mirna, 
                                    Named("sites") = slc_sites);
  return(df);
}