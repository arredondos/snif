from snif import *
import subprocess

strategies = [
	OptimizationStrategies.best1bin, 
	OptimizationStrategies.best1exp, 
	OptimizationStrategies.rand1exp, 
	OptimizationStrategies.randtobest1exp, 
	OptimizationStrategies.currenttobest1exp, 
	OptimizationStrategies.best2exp, 
	OptimizationStrategies.rand2exp, 
	OptimizationStrategies.randtobest1bin, 
	OptimizationStrategies.currenttobest1bin, 
	OptimizationStrategies.best2bin, 
	OptimizationStrategies.rand2bin, 
	OptimizationStrategies.rand1bin, 
]

#for i in range(7, 12):
opt_params = OptimizationParameters(
    strategy = OptimizationStrategies.best2exp.name,
    maxiter = 10000,
    popsize = 15,
    tol = 0.01,
    mutation = (0.5, 1),
    recombination = 0.7,
    workers = 1
)

inf_params = InferenceParameters(
    data_source = './.results/simulations2/difficult-models/191202-165044_sc04NN_D_s500_c04YN_dA_w100_Q_simulated-model-144.json',
    source_type = SourceType.JSON,
    IICR_type = IICRType.Exact,
    ms_reference_size = 5000,
    ms_simulations = int(2e4),
    psmc_mutation_rate = 1.25e-8,
    psmc_number_of_sequences = 100,
    psmc_length_of_sequences = int(2e7),
    infer_scale = False,
    data_cutoff_bounds = (1e-2, 1e2),
    data_time_intervals = 64,
    distance_function = ErrorFunction.ApproximatePDF,
    distance_parameter = 1.0,
    distance_max_allowed = 1e-10,
    distance_computation_interval = (1e-2, 1e2),
    rounds_per_test_bounds = (2, 100),
    repetitions_per_test = 1,
    number_of_components = 4,	
    bounds_islands = (2, 50),
    bounds_migrations_rates = (0.05, 60),
    bounds_deme_sizes = (1, 1),
    bounds_event_times = (1e-2, 1e2),
    bounds_effective_size = (0, 0)
)

settings = Settings(
    static_library_location = './libs/libsnif.so',
    custom_filename_tag = 'max-iter-2',
    output_directory = './.results/strategies-tests/round 2-B',
    default_output_dirname = '_SNIF_results'
)

basename = infer(inf = inf_params, opt = opt_params, settings = settings)

# basename = snif.simulate_and_infer(inf = inf_params, sim=sim_parameters, settings = settings)
# subprocess.run(['python3', 'snif2tex.py', basename])