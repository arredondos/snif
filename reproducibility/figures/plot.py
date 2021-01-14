from snif2tex import *

basenames = [
    # './.results/validation_main/XQ/NN/sc01NN_C_s400_c01NN_dA_w100_XQ',
    # './.results/validation_main/XQ/NN/sc02NN_C_s400_c02NN_dA_w100_XQ',
    # './.results/validation_main/XQ/NN/sc03NN_C_s400_c03NN_dA_w100_XQ',
    # './.results/validation_main/XQ/NN/sc04NN_C_s487_c04NN_dA_w100_XQ',
    # './.results/validation_main/XQ/NN/sc05NN_C_s445_c05NN_dA_w100_XQ',
    # './.results/validation_main/XQ/NN/sc06NN_C_s445_c06NN_dA_w100_XQ',

    # './.results/validation_main/XQ/NY/sc01NY_C_s400_c01NY_dA_w100_XQ',
    # './.results/validation_main/XQ/NY/sc02NY_C_s400_c02NY_dA_w100_XQ',
    # './.results/validation_main/XQ/NY/sc03NY_C_s391_c03NY_dA_w100_XQ',
    # './.results/validation_main/XQ/NY/sc04NY_C_s430_c04NY_dA_w100_XQ',
    # './.results/validation_main/XQ/NY/sc05NY_C_s438_c05NY_dA_w100_XQ',
    # './.results/validation_main/XQ/NY/sc06NY_C_s416_c06NY_dA_w100_XQ',

    './.results/validation_main/XQ/YN/sc01YN_C_s400_c01YN_dA_w100_XQ',
    './.results/validation_main/XQ/YN/sc02YN_C_s400_c02YN_dA_w100_XQ',
    './.results/validation_main/XQ/YN/sc03YN_C_s391_c03YN_dA_w100_XQ',
    './.results/validation_main/XQ/YN/sc04YN_C_s445_c04YN_dA_w100_XQ',
    './.results/validation_main/XQ/YN/sc05YN_C_s439_c05YN_dA_w100_XQ',
    './.results/validation_main/XQ/YN/sc06YN_C_s435_c06YN_dA_w100_XQ',

    # './.results/validation_main/XQ/YY/sc02YY_C_s400_c02YY_dA_w100_XQ',
    # './.results/validation_main/XQ/YY/sc03YY_C_s391_c03YY_dA_w100_XQ',
    # './.results/validation_main/XQ/YY/sc04YY_C_s405_c04YY_dA_w100_XQ',
    # './.results/validation_main/XQ/YY/sc05YY_C_s411_c05YY_dA_w100_XQ',
    # './.results/validation_main/XQ/YY/sc06YY_C_s405_c06YY_dA_w100_XQ',
]

for bn in basenames:

    config = Configuration(
        SNIF_basename = bn,
        plot_width_cm = 3,
        plot_height_cm = 3,
        IICR_plots_style  = OutputStyle.Excluded,
        PDF_plots_style = OutputStyle.Excluded,
        CDF_plots_style = OutputStyle.Excluded,
        islands_plot_style = OutputStyle.Excluded,
        Nref_plot_style = OutputStyle.Excluded,
        test_numbers = "all",
        one_file_per_test = True,
        versus_plot_style = OutputStyle.Full,
        CG_style = OutputStyle.Excluded,
        CG_size_history = False,
        CG_source = '',
        CG_source_legend = '',
        Nref_histograms_bins = 100,
        islands_histograms_bins = 100,
        time_histograms_bins = 100,
        migration_histograms_bins = 100,
        size_histograms_bins = 100,
        scaling_units = TimeScale.Generations,
        generation_time = 20,
    )

    TeXify(config)