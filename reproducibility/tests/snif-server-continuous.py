def get_sim_sizes(c, scale, sizes):
    if sizes:
        if scale:
            return [(1, 1)] + [(0.1, 10)] * (c - 1)
        else:
            return [(0.1, 10)] * c

    return [(1, 1)] * c

def get_inf_sizes(c, scale, sizes):
    if sizes:
        if scale:
            return [(1, 1)] + [(0.05, 20)] * (c - 1)
        else:
            return [(0.05, 20)] * c

    return [(1, 1)] * c

def get_sim_interval(scale):
    return (2e2, 1e5) if scale else (0.1, 50)

def get_inf_interval(scale):
    return (2e1, 2e5) if scale else (1e-2, 1e2)

if __name__ == '__main__':
    import argparse
    from snif import *
    
    parser = argparse.ArgumentParser()
    parser.add_argument('components', help = 'Number of components for simulation and inference')
    parser.add_argument('tests', help = 'The number of tests to perform')
    parser.add_argument('rounds', help = 'Maximum inference rounds to perform per each test')
    args = parser.parse_args()
    
    c = int(args.components)
    s = int(args.tests)
    r = int(args.rounds)

    for sizes in [True, False]:
        for scale in [True, False]:
            sim_parameters = SimulationParameters(
                sampling_strategy = SamplingStrategy.Continuous,
                total_samples = s,
                simulate_scale = scale,
                repetitions_per_scenario = 1,
                number_of_components = c,
                islands = (2, 40),
                migrations_rates = [(0.1, 50)] * c,
                deme_sizes = get_sim_sizes(c, scale, sizes),
                event_times = [get_sim_interval(scale)] * (c - 1),
                effective_size = (1000, 1000) if scale else (0, 0),
                # data_cutoff_bounds = "N/A",
                # data_time_intervals = "N/A",
                # ms_samples = "N/A",
                # ms_reference_size = "N/A",
            )

            inf_params = InferenceParameters(
                data_source = './.data/human-psmc/',
                source_type = SourceType.Simulation,
                IICR_type = IICRType.Exact,
                infer_scale = scale,
                data_cutoff_bounds = get_inf_interval(scale),
                data_time_intervals = 64,
                distance_function = ErrorFunction.ApproximatePDF,
                distance_parameter = 1.0,
                distance_max_allowed = 1e-7 if scale else 1e-10,
                distance_computation_interval = get_inf_interval(scale),
                rounds_per_test_bounds = (2, r),
                number_of_components = c,	
                bounds_islands = (2, 50),
                bounds_migrations_rates = [(0.05, 60)] * c,
                bounds_deme_sizes = get_inf_sizes(c, scale, sizes),
                bounds_event_times = [get_inf_interval(scale)] * (c - 1),
                bounds_effective_size = (100, 10000) if scale else (0, 0),
                ms_simulations = 0,
                ms_reference_size = 5000,
                psmc_mutation_rate = 1.25e-8,
                psmc_number_of_sequences = 100,
                psmc_length_of_sequences = int(2e7),
            )

            opt_params = OptimizationParameters(
                strategy = OptimizationStrategies.best2exp.name,
                maxiter = 1000,
                popsize = 15,
                tol = 0.01,
                mutation = (0.5, 1),
                recombination = 0.7,
                workers = 1
            )

            from numpy.random import random

            settings = Settings(
                static_library_location = './libs/libsnif.so',
                custom_filename_tag = '{:06d}'.format(int(random() * 1e6)),
                output_directory = './.results/main-validations',
                default_output_dirname = '_SNIF_results'
            )

            basename = simulate_and_infer(
                inf = inf_params, 
                sim = sim_parameters, 
                settings = settings,
                opt = opt_params
            )