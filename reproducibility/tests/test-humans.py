from snif import *

for c in [5]:
    for w in [1, 0.5, 0.2]:
        for restrict in [False, True]:

            low_dci = 5e4/25 if restrict else 1e4/25
            h_inference_parameters = InferenceParameters(
                data_source = './.data/human-psmc/',
                source_type = SourceType.PSMC,
                IICR_type = IICRType.Exact,
                ms_reference_size = 5000,
                ms_simulations = int(1e5),
                psmc_mutation_rate = 1.25e-8,
                psmc_number_of_sequences = 100,
                psmc_length_of_sequences = int(2e6),	
                infer_scale = True,
                data_cutoff_bounds = (1e4/25, 2e7/25),
                data_time_intervals = 64,
                distance_function = ErrorFunction.ApproximatePDF,
                distance_parameter = w,
                distance_max_allowed = 7e3,
                distance_computation_interval = (low_dci, 6e6/25),
                rounds_per_test_bounds = (50, 50),
                repetitions_per_test = 1,
                number_of_components = c,	
                bounds_islands = (2, 100),
                bounds_migrations_rates = (0.05, 20),
                bounds_deme_sizes = (1, 1),
                bounds_event_times = (low_dci, 6e6/25),
                bounds_effective_size = (100, 10_000)
            )

            h_settings = Settings(
                static_library_location = './libs/libsnif.so',
                custom_filename_tag = '_yr' if restrict else '_nr',
                output_directory = './.results/humans_4',
                default_output_dirname = '_SNIF_results'
            )

            infer(inf = h_inference_parameters, settings = h_settings)