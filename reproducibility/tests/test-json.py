from snif import *
from snif2tex import *

json_inference_parameters = InferenceParameters(
    data_source = './.results/validation_main/XQ/NN/sc06NN_C_s445_c06NN_dA_w100_XQ_simulated-model-032.json',
    source_type = SourceType.JSON,
    IICR_type = IICRType.Exact,
    ms_reference_size = 0,
    ms_simulations = int(1e6),
    psmc_mutation_rate = 1.25e-8,
    psmc_number_of_sequences = 100,
    psmc_length_of_sequences = int(2e6),	
    infer_scale = False,
    data_cutoff_bounds = (0.01, 100),
    data_time_intervals = 64,
    distance_function = ErrorFunction.ApproximatePDF,
    distance_parameter = 1,
    distance_max_allowed = 1e-10,
    distance_computation_interval = (0.01, 100),
    rounds_per_test_bounds = (500, 5000),
    repetitions_per_test = 1,
    number_of_components = 6,	
    bounds_islands = (2, 50),
    bounds_migrations_rates = (0.05, 60),
    bounds_deme_sizes = (1, 1),
    bounds_event_times = (0.01, 100),
    bounds_effective_size = (0, 0)
)

json_settings = Settings(
    static_library_location = './libs/libsnif.so',
    custom_filename_tag = 'states-inc',
    output_directory = './.results/validation_main/XQ/NN/selected',
    default_output_dirname = '_SNIF_results'
)

basename = infer(
    inf = json_inference_parameters, 
    settings = json_settings,
)

config = Configuration(
    SNIF_basename = basename,
    plot_width_cm = 13,
    plot_height_cm = 6,
    IICR_plots_style  = OutputStyle.Full,
    PDF_plots_style = OutputStyle.Excluded,
    CDF_plots_style = OutputStyle.Excluded,
    islands_plot_style = OutputStyle.Excluded,
    Nref_plot_style = OutputStyle.Excluded,
    test_numbers = "all",
    one_file_per_test = False,
    versus_plot_style = OutputStyle.Excluded,
    CG_style = OutputStyle.Excluded,
    CG_size_history = False,
    CG_source = '',
    CG_source_legend = '',
    Nref_histograms_bins = 100,
    islands_histograms_bins = 100,
    time_histograms_bins = 100,
    migration_histograms_bins = 100,
    size_histograms_bins = 100,
    scaling_units = TimeScale.Coalescent,
    generation_time = 25,
)

TeXify(config)