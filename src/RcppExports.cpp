// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// calculate_allele_spectrum_cpp
NumericMatrix calculate_allele_spectrum_cpp(NumericVector v1, double step_size, bool progress_bar);
RcppExport SEXP _isoSIM_calculate_allele_spectrum_cpp(SEXP v1SEXP, SEXP step_sizeSEXP, SEXP progress_barSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type v1(v1SEXP);
    Rcpp::traits::input_parameter< double >::type step_size(step_sizeSEXP);
    Rcpp::traits::input_parameter< bool >::type progress_bar(progress_barSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_allele_spectrum_cpp(v1, step_size, progress_bar));
    return rcpp_result_gen;
END_RCPP
}
// create_population_cpp
List create_population_cpp(int pop_size, int number_of_founders, int total_runtime, double morgan, bool progress_bar, bool track_junctions);
RcppExport SEXP _isoSIM_create_population_cpp(SEXP pop_sizeSEXP, SEXP number_of_foundersSEXP, SEXP total_runtimeSEXP, SEXP morganSEXP, SEXP progress_barSEXP, SEXP track_junctionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type pop_size(pop_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_founders(number_of_foundersSEXP);
    Rcpp::traits::input_parameter< int >::type total_runtime(total_runtimeSEXP);
    Rcpp::traits::input_parameter< double >::type morgan(morganSEXP);
    Rcpp::traits::input_parameter< bool >::type progress_bar(progress_barSEXP);
    Rcpp::traits::input_parameter< bool >::type track_junctions(track_junctionsSEXP);
    rcpp_result_gen = Rcpp::wrap(create_population_cpp(pop_size, number_of_founders, total_runtime, morgan, progress_bar, track_junctions));
    return rcpp_result_gen;
END_RCPP
}
// create_isofemale_line_cpp
List create_isofemale_line_cpp(NumericVector v, int pop_size, int total_runtime, double morgan, bool progress_bar);
RcppExport SEXP _isoSIM_create_isofemale_line_cpp(SEXP vSEXP, SEXP pop_sizeSEXP, SEXP total_runtimeSEXP, SEXP morganSEXP, SEXP progress_barSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type pop_size(pop_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type total_runtime(total_runtimeSEXP);
    Rcpp::traits::input_parameter< double >::type morgan(morganSEXP);
    Rcpp::traits::input_parameter< bool >::type progress_bar(progress_barSEXP);
    rcpp_result_gen = Rcpp::wrap(create_isofemale_line_cpp(v, pop_size, total_runtime, morgan, progress_bar));
    return rcpp_result_gen;
END_RCPP
}
// create_two_populations_migration_cpp
List create_two_populations_migration_cpp(int pop_size, int number_of_founders, int total_runtime, double morgan, double migration, bool progress_bar);
RcppExport SEXP _isoSIM_create_two_populations_migration_cpp(SEXP pop_sizeSEXP, SEXP number_of_foundersSEXP, SEXP total_runtimeSEXP, SEXP morganSEXP, SEXP migrationSEXP, SEXP progress_barSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type pop_size(pop_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_founders(number_of_foundersSEXP);
    Rcpp::traits::input_parameter< int >::type total_runtime(total_runtimeSEXP);
    Rcpp::traits::input_parameter< double >::type morgan(morganSEXP);
    Rcpp::traits::input_parameter< double >::type migration(migrationSEXP);
    Rcpp::traits::input_parameter< bool >::type progress_bar(progress_barSEXP);
    rcpp_result_gen = Rcpp::wrap(create_two_populations_migration_cpp(pop_size, number_of_founders, total_runtime, morgan, migration, progress_bar));
    return rcpp_result_gen;
END_RCPP
}
// test_fish_functions
void test_fish_functions();
RcppExport SEXP _isoSIM_test_fish_functions() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    test_fish_functions();
    return R_NilValue;
END_RCPP
}
// create_population_selection_cpp
List create_population_selection_cpp(NumericMatrix select, int pop_size, int number_of_founders, int total_runtime, double morgan, bool progress_bar, bool track_frequency, bool multiplicative_selection);
RcppExport SEXP _isoSIM_create_population_selection_cpp(SEXP selectSEXP, SEXP pop_sizeSEXP, SEXP number_of_foundersSEXP, SEXP total_runtimeSEXP, SEXP morganSEXP, SEXP progress_barSEXP, SEXP track_frequencySEXP, SEXP multiplicative_selectionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type select(selectSEXP);
    Rcpp::traits::input_parameter< int >::type pop_size(pop_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_founders(number_of_foundersSEXP);
    Rcpp::traits::input_parameter< int >::type total_runtime(total_runtimeSEXP);
    Rcpp::traits::input_parameter< double >::type morgan(morganSEXP);
    Rcpp::traits::input_parameter< bool >::type progress_bar(progress_barSEXP);
    Rcpp::traits::input_parameter< bool >::type track_frequency(track_frequencySEXP);
    Rcpp::traits::input_parameter< bool >::type multiplicative_selection(multiplicative_selectionSEXP);
    rcpp_result_gen = Rcpp::wrap(create_population_selection_cpp(select, pop_size, number_of_founders, total_runtime, morgan, progress_bar, track_frequency, multiplicative_selection));
    return rcpp_result_gen;
END_RCPP
}
// select_population_cpp
List select_population_cpp(Rcpp::NumericVector v1, Rcpp::NumericMatrix selectM, int population_size, int run_time, double morgan, bool progress_bar, bool track_frequency, bool multiplicative_selection);
RcppExport SEXP _isoSIM_select_population_cpp(SEXP v1SEXP, SEXP selectMSEXP, SEXP population_sizeSEXP, SEXP run_timeSEXP, SEXP morganSEXP, SEXP progress_barSEXP, SEXP track_frequencySEXP, SEXP multiplicative_selectionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type v1(v1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type selectM(selectMSEXP);
    Rcpp::traits::input_parameter< int >::type population_size(population_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type run_time(run_timeSEXP);
    Rcpp::traits::input_parameter< double >::type morgan(morganSEXP);
    Rcpp::traits::input_parameter< bool >::type progress_bar(progress_barSEXP);
    Rcpp::traits::input_parameter< bool >::type track_frequency(track_frequencySEXP);
    Rcpp::traits::input_parameter< bool >::type multiplicative_selection(multiplicative_selectionSEXP);
    rcpp_result_gen = Rcpp::wrap(select_population_cpp(v1, selectM, population_size, run_time, morgan, progress_bar, track_frequency, multiplicative_selection));
    return rcpp_result_gen;
END_RCPP
}
// simulate_cpp
List simulate_cpp(Rcpp::NumericVector input_population, NumericMatrix select, int pop_size, int number_of_founders, int total_runtime, double morgan, bool progress_bar, bool track_frequency, NumericVector track_markers, bool track_junctions, bool multiplicative_selection);
RcppExport SEXP _isoSIM_simulate_cpp(SEXP input_populationSEXP, SEXP selectSEXP, SEXP pop_sizeSEXP, SEXP number_of_foundersSEXP, SEXP total_runtimeSEXP, SEXP morganSEXP, SEXP progress_barSEXP, SEXP track_frequencySEXP, SEXP track_markersSEXP, SEXP track_junctionsSEXP, SEXP multiplicative_selectionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type input_population(input_populationSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type select(selectSEXP);
    Rcpp::traits::input_parameter< int >::type pop_size(pop_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_founders(number_of_foundersSEXP);
    Rcpp::traits::input_parameter< int >::type total_runtime(total_runtimeSEXP);
    Rcpp::traits::input_parameter< double >::type morgan(morganSEXP);
    Rcpp::traits::input_parameter< bool >::type progress_bar(progress_barSEXP);
    Rcpp::traits::input_parameter< bool >::type track_frequency(track_frequencySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type track_markers(track_markersSEXP);
    Rcpp::traits::input_parameter< bool >::type track_junctions(track_junctionsSEXP);
    Rcpp::traits::input_parameter< bool >::type multiplicative_selection(multiplicative_selectionSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_cpp(input_population, select, pop_size, number_of_founders, total_runtime, morgan, progress_bar, track_frequency, track_markers, track_junctions, multiplicative_selection));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_isoSIM_calculate_allele_spectrum_cpp", (DL_FUNC) &_isoSIM_calculate_allele_spectrum_cpp, 3},
    {"_isoSIM_create_population_cpp", (DL_FUNC) &_isoSIM_create_population_cpp, 6},
    {"_isoSIM_create_isofemale_line_cpp", (DL_FUNC) &_isoSIM_create_isofemale_line_cpp, 5},
    {"_isoSIM_create_two_populations_migration_cpp", (DL_FUNC) &_isoSIM_create_two_populations_migration_cpp, 6},
    {"_isoSIM_test_fish_functions", (DL_FUNC) &_isoSIM_test_fish_functions, 0},
    {"_isoSIM_create_population_selection_cpp", (DL_FUNC) &_isoSIM_create_population_selection_cpp, 8},
    {"_isoSIM_select_population_cpp", (DL_FUNC) &_isoSIM_select_population_cpp, 8},
    {"_isoSIM_simulate_cpp", (DL_FUNC) &_isoSIM_simulate_cpp, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_isoSIM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
