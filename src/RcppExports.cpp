// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// calculate_allele_spectrum_cpp
NumericMatrix calculate_allele_spectrum_cpp(NumericVector v1, int numFounders, double step_size);
RcppExport SEXP _isoSIM_calculate_allele_spectrum_cpp(SEXP v1SEXP, SEXP numFoundersSEXP, SEXP step_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type v1(v1SEXP);
    Rcpp::traits::input_parameter< int >::type numFounders(numFoundersSEXP);
    Rcpp::traits::input_parameter< double >::type step_size(step_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_allele_spectrum_cpp(v1, numFounders, step_size));
    return rcpp_result_gen;
END_RCPP
}
// simulate_from_population
List simulate_from_population(std::string file_name, int total_runtime, double morgan, int number_of_markers, int seed);
RcppExport SEXP _isoSIM_simulate_from_population(SEXP file_nameSEXP, SEXP total_runtimeSEXP, SEXP morganSEXP, SEXP number_of_markersSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type file_name(file_nameSEXP);
    Rcpp::traits::input_parameter< int >::type total_runtime(total_runtimeSEXP);
    Rcpp::traits::input_parameter< double >::type morgan(morganSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_markers(number_of_markersSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_from_population(file_name, total_runtime, morgan, number_of_markers, seed));
    return rcpp_result_gen;
END_RCPP
}
// create_population_cpp
List create_population_cpp(int pop_size, int number_of_founders, int total_runtime, double morgan, int seed, bool writeToFile);
RcppExport SEXP _isoSIM_create_population_cpp(SEXP pop_sizeSEXP, SEXP number_of_foundersSEXP, SEXP total_runtimeSEXP, SEXP morganSEXP, SEXP seedSEXP, SEXP writeToFileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type pop_size(pop_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_founders(number_of_foundersSEXP);
    Rcpp::traits::input_parameter< int >::type total_runtime(total_runtimeSEXP);
    Rcpp::traits::input_parameter< double >::type morgan(morganSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< bool >::type writeToFile(writeToFileSEXP);
    rcpp_result_gen = Rcpp::wrap(create_population_cpp(pop_size, number_of_founders, total_runtime, morgan, seed, writeToFile));
    return rcpp_result_gen;
END_RCPP
}
// create_isofemale_line_cpp
List create_isofemale_line_cpp(NumericVector v, int pop_size, int total_runtime, double morgan, int seed);
RcppExport SEXP _isoSIM_create_isofemale_line_cpp(SEXP vSEXP, SEXP pop_sizeSEXP, SEXP total_runtimeSEXP, SEXP morganSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type pop_size(pop_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type total_runtime(total_runtimeSEXP);
    Rcpp::traits::input_parameter< double >::type morgan(morganSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(create_isofemale_line_cpp(v, pop_size, total_runtime, morgan, seed));
    return rcpp_result_gen;
END_RCPP
}
// create_two_populations_cpp
List create_two_populations_cpp(int pop_size, int number_of_founders, int total_runtime, double morgan, int seed, double overlap, bool writeToFile);
RcppExport SEXP _isoSIM_create_two_populations_cpp(SEXP pop_sizeSEXP, SEXP number_of_foundersSEXP, SEXP total_runtimeSEXP, SEXP morganSEXP, SEXP seedSEXP, SEXP overlapSEXP, SEXP writeToFileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type pop_size(pop_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_founders(number_of_foundersSEXP);
    Rcpp::traits::input_parameter< int >::type total_runtime(total_runtimeSEXP);
    Rcpp::traits::input_parameter< double >::type morgan(morganSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< double >::type overlap(overlapSEXP);
    Rcpp::traits::input_parameter< bool >::type writeToFile(writeToFileSEXP);
    rcpp_result_gen = Rcpp::wrap(create_two_populations_cpp(pop_size, number_of_founders, total_runtime, morgan, seed, overlap, writeToFile));
    return rcpp_result_gen;
END_RCPP
}
// create_two_populations_migration_cpp
List create_two_populations_migration_cpp(int pop_size, int number_of_founders, int total_runtime, double morgan, int seed, double migration, bool writeToFile);
RcppExport SEXP _isoSIM_create_two_populations_migration_cpp(SEXP pop_sizeSEXP, SEXP number_of_foundersSEXP, SEXP total_runtimeSEXP, SEXP morganSEXP, SEXP seedSEXP, SEXP migrationSEXP, SEXP writeToFileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type pop_size(pop_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_founders(number_of_foundersSEXP);
    Rcpp::traits::input_parameter< int >::type total_runtime(total_runtimeSEXP);
    Rcpp::traits::input_parameter< double >::type morgan(morganSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< double >::type migration(migrationSEXP);
    Rcpp::traits::input_parameter< bool >::type writeToFile(writeToFileSEXP);
    rcpp_result_gen = Rcpp::wrap(create_two_populations_migration_cpp(pop_size, number_of_founders, total_runtime, morgan, seed, migration, writeToFile));
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
// select_population_cpp
List select_population_cpp(Rcpp::NumericVector v1, Rcpp::NumericVector selectM, double s, int population_size, int run_time, double morgan, int seed, bool writeToFile);
RcppExport SEXP _isoSIM_select_population_cpp(SEXP v1SEXP, SEXP selectMSEXP, SEXP sSEXP, SEXP population_sizeSEXP, SEXP run_timeSEXP, SEXP morganSEXP, SEXP seedSEXP, SEXP writeToFileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type v1(v1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type selectM(selectMSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type population_size(population_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type run_time(run_timeSEXP);
    Rcpp::traits::input_parameter< double >::type morgan(morganSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< bool >::type writeToFile(writeToFileSEXP);
    rcpp_result_gen = Rcpp::wrap(select_population_cpp(v1, selectM, s, population_size, run_time, morgan, seed, writeToFile));
    return rcpp_result_gen;
END_RCPP
}
// create_population_selection_cpp
List create_population_selection_cpp(int pop_size, int number_of_founders, int total_runtime, double morgan, Rcpp::NumericVector select_matrix, double selection, int seed, bool write_to_file);
RcppExport SEXP _isoSIM_create_population_selection_cpp(SEXP pop_sizeSEXP, SEXP number_of_foundersSEXP, SEXP total_runtimeSEXP, SEXP morganSEXP, SEXP select_matrixSEXP, SEXP selectionSEXP, SEXP seedSEXP, SEXP write_to_fileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type pop_size(pop_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_founders(number_of_foundersSEXP);
    Rcpp::traits::input_parameter< int >::type total_runtime(total_runtimeSEXP);
    Rcpp::traits::input_parameter< double >::type morgan(morganSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type select_matrix(select_matrixSEXP);
    Rcpp::traits::input_parameter< double >::type selection(selectionSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< bool >::type write_to_file(write_to_fileSEXP);
    rcpp_result_gen = Rcpp::wrap(create_population_selection_cpp(pop_size, number_of_founders, total_runtime, morgan, select_matrix, selection, seed, write_to_file));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_isoSIM_calculate_allele_spectrum_cpp", (DL_FUNC) &_isoSIM_calculate_allele_spectrum_cpp, 3},
    {"_isoSIM_simulate_from_population", (DL_FUNC) &_isoSIM_simulate_from_population, 5},
    {"_isoSIM_create_population_cpp", (DL_FUNC) &_isoSIM_create_population_cpp, 6},
    {"_isoSIM_create_isofemale_line_cpp", (DL_FUNC) &_isoSIM_create_isofemale_line_cpp, 5},
    {"_isoSIM_create_two_populations_cpp", (DL_FUNC) &_isoSIM_create_two_populations_cpp, 7},
    {"_isoSIM_create_two_populations_migration_cpp", (DL_FUNC) &_isoSIM_create_two_populations_migration_cpp, 7},
    {"_isoSIM_test_fish_functions", (DL_FUNC) &_isoSIM_test_fish_functions, 0},
    {"_isoSIM_select_population_cpp", (DL_FUNC) &_isoSIM_select_population_cpp, 8},
    {"_isoSIM_create_population_selection_cpp", (DL_FUNC) &_isoSIM_create_population_selection_cpp, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_isoSIM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
