from snif2tex import *

bins_per_cm = 10
v_cm = 4.5
h_cm = 8.5

config = Configuration(
    SNIF_basename = './.results/validation_main/XT/sc05YN_D_s099_c05YN_dA_w100_XT',
    plot_width_cm = h_cm,
    plot_height_cm = v_cm,
	IICR_plots_style  = OutputStyle.Excluded,
	PDF_plots_style = OutputStyle.Excluded,
	CDF_plots_style = OutputStyle.Excluded,
    islands_plot_style = OutputStyle.Full,
    Nref_plot_style = OutputStyle.Full,
    test_numbers = "all",
    one_file_per_test = False,
    versus_plot_style = OutputStyle.Excluded,
	CG_style = OutputStyle.Full,
    CG_size_history = False,
	CG_source = '',
	CG_source_legend = '',
	Nref_histograms_bins = int(v_cm * bins_per_cm),
	islands_histograms_bins = int(v_cm * bins_per_cm),
	time_histograms_bins = int(h_cm * bins_per_cm),
	migration_histograms_bins = int(v_cm * bins_per_cm),
	size_histograms_bins = int(v_cm * bins_per_cm),
    scaling_units = TimeScale.Generations,
    generation_time = 20,
)

TeXify(config)