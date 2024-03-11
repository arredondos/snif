from ctypes import cdll, c_int, c_double, c_char_p
from scipy.optimize import differential_evolution
from collections import namedtuple
from math import log10, isnan
from enum import Enum, auto
from sys import float_info
from os import path, mkdir
from warnings import warn
import atexit

class SourceType(Enum):
	PSMC, MSMC, JSON, MSCommand, Simulation, IICR = range(6)

class IICRType(Enum):
	Exact, T_sim, Seq_sim = range(3)

class ErrorFunction(Enum):
	Visual, ExactPDF, ApproximatePDF = range(3)

InferenceParameters = namedtuple('InferenceParameters', [
	'data_source',
	'source_type',
	'IICR_type',
	'infer_scale',
	'data_cutoff_bounds',
	'data_time_intervals',
	'distance_function',
	'distance_parameter',
	'distance_max_allowed',
	'distance_computation_interval',
	'rounds_per_test_bounds',
	'repetitions_per_test',
	'number_of_components',
	'bounds_islands',
	'bounds_migrations_rates',
	'bounds_deme_sizes',
	'bounds_event_times',
	'bounds_effective_size',
	'psmc_mutation_rate',
    'psmc_number_of_sequences',
    'psmc_length_of_sequences',
	'ms_reference_size',
	'ms_simulations'
])

# if infer_scale == False, then
# data_cutoff_bounds, distance_computation_interval, and
# bounds_event_times are interpreted in the coalescent 
# time scale. Otherwise, they are interpreted in generations.
default_inference_parameters = InferenceParameters(
	data_source = './.data/human-psmc/Dai_Upper.psmc',
	source_type = SourceType.PSMC,
	IICR_type = IICRType.Exact,
	ms_reference_size = 5000,
	ms_simulations = int(1e5),
	psmc_mutation_rate = 1.25e-8,
    psmc_number_of_sequences = 100,
    psmc_length_of_sequences = int(2e6),	
	infer_scale = True,
	data_cutoff_bounds = (3e2, 2e5),
	data_time_intervals = 64,
	distance_function = ErrorFunction.ApproximatePDF,
	distance_parameter = 1,
	distance_max_allowed = 7e3,
	distance_computation_interval = (4e2, 1e5),
	rounds_per_test_bounds = (1, 3),
	repetitions_per_test = 1,
	number_of_components = 3,	
	bounds_islands = (2, 100),
	bounds_migrations_rates = (0.05, 20),
	bounds_deme_sizes = (1, 1),
	bounds_event_times = (4e2, 1e5),
	bounds_effective_size = (100, 100000)
)

class SamplingStrategy(Enum):
	Continuous, Discrete = range(2)

SimulationParameters = namedtuple('SimulationParameters', [
	'sampling_strategy',
	'total_samples',
	'simulate_scale',
	'number_of_components',
	'islands',
	'migrations_rates',
	'deme_sizes',
	'event_times',
	'effective_size',
])

# if simulate_scale == False, then
# data_cutoff_bounds and event_times 
# are interpreted in the coalescent 
# time scale. Otherwise, they are 
# interpreted in generations.
default_simulation_parameters = SimulationParameters(
	sampling_strategy = SamplingStrategy.Continuous,
	total_samples = 5,
	simulate_scale = False,
	number_of_components = 4,
	islands = (2, 30),
	migrations_rates = (0.05, 10),
	deme_sizes = (1, 1),
	event_times = (1e-2, 1e2),
	effective_size = (100, 100000),
)

Settings = namedtuple('Settings', [
	'static_library_location',
	'custom_filename_tag',
	'output_directory',
	'default_output_dirname'
])

default_settings = Settings(
	static_library_location = './libs/libsnif.so',
	custom_filename_tag = '',
	output_directory = './.results/',
	default_output_dirname = '_SNIF_results'
)

class OptimizationStrategies(Enum):
	best1bin = auto() 
	best1exp = auto() 
	rand1exp = auto() 
	randtobest1exp = auto() 
	currenttobest1exp = auto() 
	best2exp = auto() 
	rand2exp = auto() 
	randtobest1bin = auto() 
	currenttobest1bin = auto() 
	best2bin = auto() 
	rand2bin = auto() 
	rand1bin = auto() 

OptimizationParameters = namedtuple('OptimizationParameters', [
	'strategy',
	'maxiter',
	'popsize',
	'tol',
	'mutation',
	'recombination',
	'workers'
])

default_optimization_parameters = OptimizationParameters(
	 strategy = OptimizationStrategies.best2exp.name,
	 maxiter = 1000,
	 popsize = 15,
	 tol = 0.01,
	 mutation = (0.5, 1),
	 recombination = 0.7,
	 workers = 1
)

## Public Interface ##

def infer(
	inf = default_inference_parameters, 
    opt = default_optimization_parameters,
	settings = default_settings):

	print('\nWelcome to SNIF -- an Inferential Framework for Demographic Structure')
	
	source_filenames = _get_source_filenames(inf.data_source, inf.source_type)
	k = len(source_filenames)
	if not (k > 0):
		print("could not find any files of type {} in {}".format(inf.source_type.name, inf.data_source))
		return

	file_string = 'files' if k > 1 else 'file'
	print('found {} source {} of type {}'.format(k, file_string, inf.source_type.name))
	
	(infer_sizes, bounds) = _build_bounds(
		inf.number_of_components, 
		inf.bounds_islands,
		inf.bounds_event_times, 
		inf.bounds_migrations_rates, 
		inf.bounds_deme_sizes,
		inf.infer_scale,
		inf.bounds_effective_size
	)	
	
	basename = _generate_inference_basename(
		inf.infer_scale, 
		infer_sizes,
		inf.number_of_components, 
		inf.distance_function,
		inf.distance_parameter, 
		inf.source_type, 
		inf.IICR_type,
		settings.custom_filename_tag
	)

	output_directory = _build_output_directory(
		inf.data_source,
		settings.output_directory, 
		settings.default_output_dirname
	)

	_load_library(settings.static_library_location)
	
	print('initializing static library..')
	snif_lib.initialize(
		c_char_p(basename.encode('utf-8')),
		c_char_p(output_directory.encode('utf-8')),
		c_int(0),
		c_int(inf.number_of_components),
		c_int(inf.distance_function.value),
		c_double(inf.distance_parameter),
		c_int(1 if inf.infer_scale else 0),
		c_int(1 if infer_sizes else 0),
		c_int(inf.data_time_intervals),
		c_double(inf.data_cutoff_bounds[0]),
		c_double(inf.data_cutoff_bounds[1]),
		c_double(inf.distance_computation_interval[0]),
		c_double(inf.distance_computation_interval[1])
	)

	objective_function = _build_objective_function(
		inf.number_of_components, 
		infer_sizes, 
		inf.infer_scale
	)

	for i in range(k):
		source_filename = source_filenames[i]
		binary_filename = source_filename.encode('utf-8')
		print('\nloading file {} of {}: {}'.format(i + 1, k, source_filename))

		if inf.source_type == SourceType.PSMC :
			snif_lib.load_smc_problem_from_output(
				c_char_p(binary_filename),
				c_double(inf.psmc_mutation_rate)
			)
		
		elif inf.source_type == SourceType.MSMC:
			raise NotImplementedError("The MSMC parser is not implemented yet")

		elif inf.source_type == SourceType.JSON:
			if inf.IICR_type == IICRType.Exact:
				snif_lib.load_exact_problem_from_description(
					c_char_p(binary_filename)
				)
			elif inf.IICR_type == IICRType.T_sim:
				snif_lib.load_t2_problem_from_description(
					c_char_p(binary_filename),
					c_int(inf.ms_simulations)
				)
			elif inf.IICR_type == IICRType.Seq_sim:
				snif_lib.load_smc_problem_from_description(
					c_char_p(binary_filename),
					c_double(inf.psmc_mutation_rate),
					c_int(inf.psmc_number_of_sequences),
					c_int(inf.psmc_length_of_sequences)
				)
			else:
				raise ValueError("Invalid IICR type string")

		elif inf.source_type == SourceType.MSCommand:			
			if inf.IICR_type == IICRType.Exact:
				raise ValueError("Cannot build an exact IICR from an ms command!")
			elif inf.IICR_type == IICRType.T_sim:
				snif_lib.load_t2_problem_from_mscommand(
					c_char_p(binary_filename),
					c_double(inf.ms_reference_size)
				)
			elif inf.IICR_type == IICRType.Seq_sim:
				snif_lib.load_smc_problem_from_mscommand(
					c_char_p(binary_filename),
					c_double(inf.psmc_mutation_rate)
				)
			else:
				raise ValueError("Invalid IICR type string")
			
		elif inf.source_type == SourceType.IICR:
			if inf.IICR_type != IICRType.Exact:
				raise ValueError(f"Cannot build a {inf.IICR_type} IICR from a curve!")
			
			snif_lib.load_exact_problem_from_curve(
				c_char_p(binary_filename),
			)

		else:
			raise ValueError("Invalid data type string")
		
		reps = inf.repetitions_per_test
		for scenario_round in range(reps):
			if reps > 1:
				print("repetition", scenario_round + 1, "of", reps)

			_run_one_inference(
				bounds,
				inf.rounds_per_test_bounds, 
				inf.distance_max_allowed, 
				objective_function, 
				opt._asdict()
			)

	inference_request = { 
		'inference_parameters' : inf._asdict(),
		'optimization_parameters' : opt._asdict(),
		'settings' : settings._asdict(),
		'source_filenames' : source_filenames,
		'final_output_directory' : output_directory,
		'basename' : basename,
	}
	inference_request['inference_parameters']['source_type'] = inf.source_type.name
	inference_request['inference_parameters']['IICR_type'] = inf.IICR_type.name
	inference_request['inference_parameters']['distance_function'] = inf.distance_function.name
	
	_write_settings_file(basename, output_directory, inference_request)
	return path.join(output_directory, basename)

def simulate_and_infer(
	inf = default_inference_parameters, 
	sim = default_simulation_parameters,
    opt = default_optimization_parameters,
	settings = default_settings):
	
	from itertools import islice	
	
	print('\nWelcome to SNIF -- an Inferential Framework for Demographic Structure')

	(infer_sizes, inf_bounds) = _build_bounds(
		inf.number_of_components, 
		inf.bounds_islands,
		inf.bounds_event_times, 
		inf.bounds_migrations_rates, 
		inf.bounds_deme_sizes,
		inf.infer_scale,
		inf.bounds_effective_size
	)	

	if sim.sampling_strategy == SamplingStrategy.Discrete:
		simulate_sizes = sim.deme_sizes != [1.0]
		
		tests_generator = _pre_selected_values_generator(
			sim.number_of_components,
			sim.simulate_scale,
			sim.islands,
			sim.event_times,
			sim.migrations_rates,
			sim.deme_sizes,
			sim.effective_size
		)
	else:
		(simulate_sizes, sim_bounds) = _build_bounds(
			sim.number_of_components, 
			sim.islands,
			sim.event_times, 
			sim.migrations_rates, 
			sim.deme_sizes,
			sim.simulate_scale,
			sim.effective_size
		)

		tests_generator = _random_values_generator(
			sim.number_of_components,
			sim.simulate_scale,
			simulate_sizes,
			sim_bounds
		)

	scenarios = list(islice(tests_generator, sim.total_samples)) 	
	print('successfully generated {} scenarios for simulation'.format(sim.total_samples))
	
	basename = _generate_simulation_and_inference_basename(
		sim.simulate_scale,
		simulate_sizes,
		sim.number_of_components,
		sim.sampling_strategy,
		sim.total_samples,
		inf.infer_scale, 
		infer_sizes,
		inf.number_of_components, 
		inf.distance_function,
		inf.distance_parameter, 
		inf.source_type, 
		inf.IICR_type,
		settings.custom_filename_tag
	)

	output_directory = _build_output_directory(
		path.curdir,
		settings.output_directory, 
		settings.default_output_dirname
	)

	_load_library(settings.static_library_location)
	print('initializing static library..')
	snif_lib.initialize(
		c_char_p(basename.encode('utf-8')),
		c_char_p(output_directory.encode('utf-8')),
		c_int(sim.number_of_components),
		c_int(inf.number_of_components),
		c_int(inf.distance_function.value),
		c_double(inf.distance_parameter),
		c_int(1 if inf.infer_scale else 0),
		c_int(1 if infer_sizes else 0),
		c_int(inf.data_time_intervals),
		c_double(inf.data_cutoff_bounds[0]),
		c_double(inf.data_cutoff_bounds[1]),
		c_double(inf.distance_computation_interval[0]),
		c_double(inf.distance_computation_interval[1])
	)

	objective_function = _build_objective_function(
		inf.number_of_components, 
		infer_sizes, 
		inf.infer_scale
	)

	test_number = 0
	c_arr = c_double * sim.number_of_components
	for scenario in scenarios:
		test_number += 1
		print("\nsimulating scenario", test_number, "of", sim.total_samples)

		if inf.IICR_type == IICRType.Exact:
			snif_lib.load_exact_simulated_problem(
				c_int(scenario[0]),
				c_arr(*scenario[1]),
				c_arr(*scenario[2]),
				c_arr(*scenario[3]),
				c_double(scenario[4])
			)

		elif inf.IICR_type == IICRType.T_sim:
			snif_lib.load_t2_simulated_problem(
				c_int(scenario[0]),
				c_arr(*scenario[1]),
				c_arr(*scenario[2]),
				c_arr(*scenario[3]),
				c_double(scenario[4]),
				c_int(inf.ms_simulations)
			)

		elif inf.IICR_type == IICRType.Seq_sim:
			snif_lib.load_smc_simulated_problem(
				c_int(scenario[0]),
				c_arr(*scenario[1]),
				c_arr(*scenario[2]),
				c_arr(*scenario[3]),
				c_double(scenario[4]),
				c_double(inf.psmc_mutation_rate),
				c_int(inf.psmc_number_of_sequences),
				c_int(inf.psmc_length_of_sequences)
			)

		else:
			raise ValueError("Invalid IICR type string")
		
		reps = inf.repetitions_per_test
		for scenario_round in range(reps):
			if reps > 1:
				print("repetition", scenario_round + 1, "of", reps)

			_run_one_inference(
				inf_bounds,
				inf.rounds_per_test_bounds, 
				inf.distance_max_allowed, 
				objective_function, 
				opt._asdict()
			)
			
	inference_request = { 
		'simulation_parameters' : sim._asdict(),
		'inference_parameters' : inf._asdict(),
		'optimization_parameters' : opt._asdict(),
		'settings' : settings._asdict(),
		'final_output_directory' : output_directory,
		'basename' : basename,
	}
	inference_request['inference_parameters']['source_type'] = inf.source_type.name
	inference_request['inference_parameters']['IICR_type'] = inf.IICR_type.name
	inference_request['inference_parameters']['distance_function'] = inf.distance_function.name
	#inference_request['simulation_parameters']['IICR_type'] = sim.IICR_type.name
	inference_request['simulation_parameters']['sampling_strategy'] = sim.sampling_strategy.name
	
	_write_settings_file(basename, output_directory, inference_request)
	return path.join(output_directory, basename)

## Private Interface ##

def _load_library(static_library_location):
	print('loading static library..')
	global snif_lib
	snif_lib = cdll.LoadLibrary(static_library_location)
	snif_lib.evaluate.restype = c_double
	snif_lib.close_problem.restype = c_double
	atexit.register(_finalization)

def _get_source_filenames(data_source, source_type):
	if path.isfile(data_source):
		return [data_source]

	if not path.isdir(data_source):
		raise ValueError('Invalid data source string')
	
	mask = ''
	if source_type == SourceType.PSMC:
		mask = '*.psmc'
	elif source_type == SourceType.MSMC:
		mask = '*.msmc'
	elif source_type == SourceType.MSCommand:
		mask = '*.ms'
	elif source_type == SourceType.JSON:
		mask = '*.json'
	else: 
		raise ValueError('Invalid source data type:', source_type)
	
	from glob import glob
	pattern = path.join(data_source, mask)
	file_list = glob(pattern)	
	file_list.sort()

	return file_list

def _build_bounds(
	number_of_components, 
	islands, 
    event_times, 
	migrations_rates, 
	deme_sizes, 
	infer_scale, 
	effective_size):

	if not (number_of_components > 0):
		raise ValueError('incorrect number of components:', number_of_components)

	all_bounds = isinstance(event_times, list) and \
			     isinstance(migrations_rates, list) and \
				 isinstance(deme_sizes, list) and \
				 len(event_times) == number_of_components - 1 and \
				 len(migrations_rates) == number_of_components and \
				 len(deme_sizes) == number_of_components
	
	if all_bounds:
		return _parse_all_bounds(
			islands, 
			event_times, 
			migrations_rates, 
			deme_sizes, 
			infer_scale, 
			effective_size
		)

	some_bounds = isinstance(event_times, tuple) and \
  				  isinstance(migrations_rates, tuple) and \
  				  isinstance(deme_sizes, tuple)

	if some_bounds:
		return _parse_some_bounds(
			number_of_components, 
			islands, 
			event_times, 
			migrations_rates, 
			deme_sizes, 
			infer_scale, 
			effective_size
		)

	raise ValueError('inconsistent specifications for the bounds and the number of components')

def _parse_all_bounds(
	islands, 
    event_times, 
	migrations_rates, 
	deme_sizes, 
	infer_scale, 
	effective_size):

	from math import log10

	(n_low, n_high) = islands
	if (n_low < 2) or (n_high < n_low):
		raise ValueError('bad bounds for the number of islands: ', islands)

	bounds = [ (n_low - 0.49, n_high + 0.49) ]

	for time in event_times:
		(t_low, t_high) = time
		if (t_low < 0) or (t_high < t_low):
			raise ValueError('bad bounds for the event times: ', time)
		bounds.append((log10(t_low), log10(t_high)))

	for mig in migrations_rates:
		(m_low, m_high) = mig
		if (m_low < 0) or (m_high < m_low):
			raise ValueError('bad bounds for the migration rates: ', mig)
		bounds.append(mig)

	infer_sizes = False
	for size in deme_sizes:
		(s_low, s_high) = size
		if size != (1, 1):
			infer_sizes = True

		if (s_low < 0) or (s_high < s_low):
			raise ValueError('bad bounds for the deme sizes: ', size)

	if infer_sizes:
		bounds.extend(deme_sizes)

	if infer_scale:
		(n0_low, n0_high) = effective_size
		if (n0_low < 0) or (n0_high < n0_low):
			raise ValueError('bad bounds for the reference size: ', effective_size)
		
		bounds.append(effective_size)

	return (infer_sizes, bounds)

def _parse_some_bounds(
	number_of_components,
	islands, 
    event_times, 
	migrations_rates, 
	deme_sizes, 
	infer_scale, 
	effective_size):
	
	from math import log10

	(n_low, n_high) = islands
	if (n_low < 2) or (n_high < n_low):
		raise ValueError('bad bounds for the number of islands: ', islands)
	bounds = [ (n_low - 0.49, n_high + 0.49) ]

	(t_low, t_high) = event_times
	if (t_low < 0) or (t_high < t_low):
		raise ValueError('bad bounds for the event times: ', event_times)
	bounds.extend([(log10(t_low), log10(t_high))] * (number_of_components - 1))

	(m_low, m_high) = migrations_rates
	if (m_low < 0) or (m_high < m_low):
		raise ValueError('bad bounds for the migration rates: ', migrations_rates)
	bounds.extend([migrations_rates] * number_of_components)

	(s_low, s_high) = deme_sizes
	if (s_low < 0) or (s_high < s_low):
		raise ValueError('bad bounds for the deme sizes: ', deme_sizes)

	infer_sizes = not (deme_sizes == (1, 1))
	if infer_sizes:
		bounds.extend([deme_sizes] * number_of_components)

	if infer_scale:
		(n0_low, n0_high) = effective_size
		if (n0_low < 0) or (n0_high < n0_low):
			raise ValueError('bad bounds for the reference size: ', effective_size)		
		bounds.append(effective_size)

	return (infer_sizes, bounds)

def _generate_inference_basename(
	infer_scale, 
	infer_sizes, 
	number_of_components,
    distance_function, 
	distance_parameter, 
	source_type, 
	IICR_type,
	custom_filename_tag):

	fn = _generate_datetime_tag()
	
	fn += _generate_inference_tag(
		infer_scale, 
		infer_sizes, 
		number_of_components,
		distance_function, 
		distance_parameter, 
		source_type,
		IICR_type,
		custom_filename_tag
	)
	
	return fn

def _generate_simulation_and_inference_basename(
	simulate_scale,
	simulate_sizes,
	simulation_components,
	sampling_strategy,
	total_samples,
	infer_scale, 
	infer_sizes, 
	number_of_components,
    distance_function, 
	distance_parameter, 
	source_type, 
	IICR_type,
	custom_filename_tag):

	fn = _generate_datetime_tag()
	
	fn += _generate_simulation_tag(
		simulate_scale,
		simulate_sizes,
		simulation_components,
		sampling_strategy,
		total_samples
	)
	
	fn += _generate_inference_tag(
		infer_scale, 
		infer_sizes, 
		number_of_components,
		distance_function, 
		distance_parameter, 
		source_type,
		IICR_type,
		custom_filename_tag
	)

	return fn

def _generate_datetime_tag():
	from time import strftime
	return strftime("%y%m%d-%H%M%S")

def _generate_simulation_tag(
	simulate_scale,
	simulate_sizes,
	simulation_components,
	sampling_strategy,
	total_samples):

	scale_string = 'Y' if simulate_scale else 'N'
	sizes_string = 'Y' if simulate_sizes else 'N'
	fn = '_sc{:02d}{}{}'.format(simulation_components, scale_string, sizes_string)
	
	sampling_string = 'D' if sampling_strategy == SamplingStrategy.Discrete else 'C'
	fn += '_{}'.format(sampling_string)
	fn += '_s{:03d}'.format(int(total_samples))

	return fn

def _generate_inference_tag(
	infer_scale,
	infer_sizes,
	number_of_components,
	distance_function,
	distance_parameter,
	source_type,
	IICR_type,
	custom_filename_tag):

	scale_string = 'Y' if infer_scale else 'N'
	sizes_string = 'Y' if infer_sizes else 'N'
	fn = '_c{:02d}{}{}'.format(number_of_components, scale_string, sizes_string)
	
	distance_strings = ['V', 'E', 'A']
	fn += '_d{}'.format(distance_strings[distance_function.value])
	
	fn += '_w{:03d}'.format(int(round(100 * distance_parameter)))

	source_strings = ['S', 'S', 'J', 'M', 'X', 'C']
	iicr_type_strings = ['Q', 'T', 'S']
	fn += '_{}{}'.format(
		source_strings[source_type.value],
		iicr_type_strings[IICR_type.value],
	)

	if not custom_filename_tag == '':
		fn += '_{}'.format(str(custom_filename_tag))
	
	return fn

def _random_values_generator(
	number_of_components,
	scaled,
	size_changes,
	bounds):

	c = number_of_components

	islands = bounds[0]
	log_times = bounds[1 : c]
	migrations = bounds[c : 2*c]

	sizes = [(1, 1)] * c
	if size_changes:
		sizes = bounds[2*c : 3*c]
	
	ref_sizes = (0, 0)
	if scaled:
		ref_sizes = bounds[-1]

	from numpy.random import randint, uniform
	while True:
		n = randint(round(islands[0]), round(islands[1]) + 1)
		t_log = [uniform(b[0], b[1]) for b in log_times]
		t = [0.0] + [10 ** x for x in sorted(t_log)]
		M = [uniform(b[0], b[1]) for b in migrations]
		s = [uniform(b[0], b[1]) for b in sizes]
		n0 = uniform(ref_sizes[0], ref_sizes[1]) if scaled else 0
		#if scaled: t = [2 * n0 * x for x in t]
		
		yield (n, t, M, s, n0)

def _pre_selected_values_generator(
	number_of_components,
	scaled,
	islands,
	event_times,
	migrations_rates,
	deme_sizes,
	effective_size):

	from random import sample, choice
	c = number_of_components

	while True:
		t_indexes = sample(range(len(event_times)), c - 1)
		t_indexes.sort()
		t_map = _get_neighbors_map(t_indexes)
		if(any(t_map)):
			continue
		
		m_indexes = [choice(range(len(migrations_rates))) for _ in range(c)]
		m_map = _get_neighbors_map(m_indexes)

		s_indexes = [choice(range(len(deme_sizes))) for _ in range(c)]
		s_map = _get_neighbors_map(s_indexes)
		combined_map = [x and y for (x, y) in zip(m_map, s_map)]
		if(any(combined_map)):
			continue

		n = choice(islands)
		t = [0] + [event_times[i] for i in t_indexes]
		m = [migrations_rates[i] for i in m_indexes]
		s = [deme_sizes[i] for i in s_indexes]
		n0 = choice(effective_size)
		if scaled: t = [2 * n0 * x for x in t]

		yield (n, t, m, s, n0)

def _get_neighbors_map(l):
	neighbors_map = [False] * len(l)

	for i in range(1, len(l)):
		if abs(l[i] - l[i - 1]) <= 1:
			neighbors_map[i] = True
	
	return neighbors_map

def _run_one_inference(
	bounds, 
    rounds_per_test_bounds, 
	distance_max_allowed, 
	objective_function, 
	opt):

	import sys
	
	(min_repeats, max_repeats) = rounds_per_test_bounds
	current_distance = float_info.max
	round_number = 0

	snif_lib.open_problem()
	while True:
		repeat = (round_number < min_repeats) or (current_distance > distance_max_allowed)
		if repeat:
			if round_number < max_repeats:
				round_number += 1
				print('  round {} of {}'.format(round_number, max_repeats))
			else:
				print('  closing problem...')
				break
		else:
			print('  closing problem with early convergence...')
			break
		snif_lib.begin_round()
		current_distance = differential_evolution(objective_function, bounds, **opt).fun
	
	fit_value = snif_lib.close_problem()
	print('Inference completed. Best fit value:', fit_value)

def _build_objective_function(number_of_components, infer_sizes, infer_scale):
	l = 2 * number_of_components
	if infer_sizes: l += number_of_components
	if infer_scale: l += 1
	
	c_arr = c_double * l
	def objective(parameters):
		return snif_lib.evaluate(c_arr(*parameters))
	
	return objective

def _build_output_directory(data_source, output_directory, default_output_dirname):
	if output_directory == '':
		output_directory = _default_output_directory(data_source, default_output_dirname)
	elif not path.isdir(output_directory):
		try:
			mkdir(output_directory)
		except OSError:
			print('could not use output directory {}. Using default directory'.format(output_directory))
			output_directory = _default_output_directory(data_source, default_output_dirname)
	
	return output_directory

def _default_output_directory(data_source, default_output_dirname):
	source_dir = path.dirname(data_source)
	default_dir = path.join(source_dir, default_output_dirname)
	if not path.isdir(default_dir):
		try:
			mkdir(default_dir)
		except OSError:
			print('could not create default directory {}. Using current directory instead'.format(default_dir))
			default_dir = path.curdir()
	
	return default_dir

def _write_settings_file(basename, output_directory, payload):
	settings_filename = path.join(output_directory, basename + '_settings.json')
	
	from json import dump
	with open(settings_filename, 'w') as json_file:
		dump(payload, json_file, indent = 4)

def _parse_basename(basename, datetime = True):
	chunks = basename.split('_')
	description = []
	i = 0

	if datetime:
		# date and time
		cc = chunks[i] 
		import datetime
		dt = datetime.datetime(
			year = int('20' + cc[:2]), 
			month = int(cc[2:4]),
			day = int(cc[4:6]),
			hour = int(cc[7:9]),
			minute = int(cc[9:11]),
			second = int(cc[11:13]))
		description.append(("Date", str(dt.date())))
		description.append(("Time", str(dt.time())))
		i += 1

	cc = chunks[i]
	if cc[0] == 's':
		# simulation components
		description.append(("sim. components", int(cc[2:-2])))
		description.append(("sim. scale", cc[-2] == 'Y'))
		description.append(("sim. sizes", cc[-1] == 'Y'))

		# sampling method
		i += 1
		cc = chunks[i]
		options = ["D", "C"]
		labels = ["Discrete", "Continuous"]
		description.append(("Sampling method", labels[options.index(cc)]))

		# sample count
		i += 1
		cc = chunks[i]
		description.append(("Sample size", int(cc[1:])))

		# next chunk
		i += 1
		cc = chunks[i]
	
	description.append(("inf. components", int(cc[1:-2])))
	description.append(("inf. scale", cc[-2] == 'Y'))
	description.append(("inf. sizes", cc[-1] == 'Y'))

	# distance
	i += 1
	cc = chunks[i]
	options = ["A", "V", "E"]
	labels = ["Approx. pdf", "Visual", "Exact pdf"]
	description.append(("Distance type", labels[options.index(cc[1])]))

	# omega
	i += 1
	cc = chunks[i]
	description.append(("Omega parameter", float(cc[1:]) / 100))

	# source and IICR type
	i += 1
	cc = chunks[i]

	options = ["S", "J", "M", "X"]
	labels = ["SMC data", "JSON description", "ms command", "simulation"]
	description.append(("Source type", labels[options.index(cc[0])]))

	options = ["Q", "T", "S"]
	labels = ["Exact IICR", "T-sim IICR", "Seq-sim IICR"]
	description.append(("IICR type", labels[options.index(cc[1])]))

	# tag
	i += 1
	tag = '_'.join(chunks[i:]) if len(chunks) > i else ""
	description.append(("Custom tag", tag))

	return description

def _finalization():
	snif_lib.finalize()
	print('\nstatic SNIF library terminated')

if __name__ == '__main__':
	infer()
	#simulate_and_infer()