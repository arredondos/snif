module source.app;

import core.runtime;
import std.stdio;
import std.datetime.stopwatch;
import std.datetime.systime;
import std.string;
import std.path;
import std.range;
import std.algorithm;
import std.conv;
import std.math;

import nsnic;
import numeric;

private:

enum string csvSeparator = ",";
enum double failDistance = double.max;
enum IICRType {
    Exact,
    T_sim,
    Seq_sim
}

// simulation and inference parameters
bool sourceModelIsKnown,
     inferScale,
     inferDemeSizes,
     simulating;

int inferenceComponents, 
    simulationComponents,
    plotIntervals,
    plotHDIntervals;

double data_cutoff_low, 
       data_cutoff_high, 
       data_weight_low, 
       data_weight_high, 
       distanceParameter;

string basename, 
       currentProblemFilename;
       
ErrorStrategy distanceType;

// local state
int problemCount, 
    roundCount, 
    evaluationCount, 
    failedEvaluationCount;

double bestDistance;
double[] fitTime, 
         fitIICR, 
         plotTime, 
         plotIICR;

Scenario bestScenario;
StopWatch stopwatch;
NSNIC sourceModel, preallocatedModel;
File csvTable, 
     csvPlots;

void resetProblemGlobals(inout string filename = null) {
	roundCount = 0;
	evaluationCount = 0;
	failedEvaluationCount = 0;
	bestDistance = failDistance;
    if (filename) {
        currentProblemFilename = filename;
    }
}

void check_for_model() {
    // look for json description of currentProblemFilename. 
    // if any, parse and build the sourceModel and call 
    // sourceModel.allocateDistanceComputation(fitTime);
    // else:
    sourceModelIsKnown = false;
    sourceModel = NSNIC.init;
}

void description_problem_intro(inout string filename) {
    resetProblemGlobals(filename);
    sourceModelIsKnown = true;
    sourceModel = NSNIC(currentProblemFilename);
}

void simulated_problem_intro(int n, double* times_p, double* migrations_p, double* sizes_p, double n0) {
    resetProblemGlobals();
    auto times = times_p[0..simulationComponents].dup;
    auto migrations = migrations_p[0..simulationComponents].dup;
    auto sizes = sizes_p[0..simulationComponents].dup;    
    sourceModel = NSNIC(n, times, migrations, sizes, n0);
}

int smc_problem_body_and_outro(double mu, int sequences, int sites) {
    import ms : getMSCommand, runPSMCAnalysis;
    auto sourceScenario = sourceModel.getScenario(false);
    auto msCommand = getMSCommand(sourceScenario.islands,
        sourceScenario.eventTimes, sourceScenario.migrationRates,
        sourceScenario.demeSizes, sourceScenario.referenceSize,
        mu, sequences, sites, true, true);

    auto psmc_filename = runPSMCAnalysis(msCommand, basename ~ "_" ~ (problemCount + 1).to!string);
    return smc_problem_outro(psmc_filename, mu, false);    
}

int smc_problem_outro(string psmc_filename, double mu, bool checkSource) {
    import ms : parsePSMCFile;

    double[] readTime, readIICR;
    parsePSMCFile(psmc_filename, mu, readTime, readIICR);    
    constrainData(readTime, readIICR, data_cutoff_low, data_cutoff_high, plotTime, plotIICR);
    constrainData(readTime, readIICR, data_weight_low, data_weight_high, fitTime, fitIICR);
    preallocatedModel.allocateDistanceComputation(fitTime);
    if(checkSource) {
        check_for_model();
    }
    else {
        sourceModel.allocateDistanceComputation(fitTime);
    }
    stopwatch.start();
	return 0;
}

int t2_problem_outro(int simulations) {
    import ms : getMSCommand, generateTsimIICR;
    auto sourceScenario = sourceModel.getScenario(false);
    auto msCommand = getMSCommand(sourceScenario.islands,
        sourceScenario.eventTimes, sourceScenario.migrationRates,
        sourceScenario.demeSizes, simulations);    

    generateTsimIICR(msCommand, data_cutoff_low, data_cutoff_high, plotIntervals, 
	    sourceScenario.referenceSize, plotTime, plotIICR, basename ~ "_" ~ (problemCount + 1).to!string);

    constrainData(plotTime, plotIICR, data_weight_low, data_weight_high, fitTime, fitIICR);
    preallocatedModel.allocateDistanceComputation(fitTime);
    sourceModel.allocateDistanceComputation(fitTime);
    stopwatch.start();
    return 0;
}

int exact_problem_outro() {
    plotTime = GetLogTimeVector(data_cutoff_low, data_cutoff_high, plotIntervals);        
    plotIICR = new double[] (plotTime.length);
    sourceModel.iicr(plotTime, plotIICR);
    constrainData(plotTime, plotIICR, data_weight_low, data_weight_high, fitTime, fitIICR);
    preallocatedModel.allocateDistanceComputation(fitTime);
    sourceModel.allocateDistanceComputation(fitTime);
    stopwatch.start();
    return 0;
}

void constrainData(double[] time, double[] iicr, double lower, double upper, 
                   out double[] constrainedTime, out double[] constrainedIICR) {

    immutable frontDrop = time.assumeSorted.lowerBound(lower).length;
    immutable backDrop  = time.assumeSorted.upperBound(upper).length;
    constrainedTime = time.drop(frontDrop).dropBack(backDrop).array;
    constrainedIICR = iicr.drop(frontDrop).dropBack(backDrop).array;    
}

void logCurve(R)(string name, R curve) {
	csvPlots.write(name ~ "-" ~ problemCount.to!string ~ csvSeparator);
    csvPlots.writeln(curve.map!(to!string).joiner(csvSeparator));
}

void logCurve(string name) {
	csvPlots.writeln(name ~ "-" ~ problemCount.to!string);
}

extern(C) export:

int initialize(
        char* p_basename, 
        char* p_targetDirectory, 
        int simulationComponents,
        int inferenceComponents, 
        int distanceType, 
        double distanceParameter, 
        int inferScale,
        int inferDemeSizes,
        int plotIntervals, 
        double data_cutoff_low, 
        double data_cutoff_high, 
        double data_weight_low, 
        double data_weight_high) {	

    Runtime.initialize();
    writeln("D Runtime initialized succesfully");

    basename = buildPath(p_targetDirectory.to!string, p_basename.to!string);
    .simulationComponents = simulationComponents;
    .inferenceComponents = inferenceComponents;
    .distanceType = distanceType.to!ErrorStrategy;
    .distanceParameter = distanceParameter;
    .inferScale = cast(bool) inferScale;
    .inferDemeSizes = cast(bool) inferDemeSizes;
    .data_cutoff_low = data_cutoff_low;
    .data_cutoff_high = data_cutoff_high;
    .data_weight_low = data_weight_low;
    .data_weight_high = data_weight_high;
    .plotIntervals = plotIntervals;
    .plotHDIntervals = plotHDIntervals;
    problemCount = 0;

    if(!(simulationComponents > 0)) {
        .simulationComponents = 0;
        simulating = false;
    }
    else {
        simulating = true;
        sourceModelIsKnown = true;
        currentProblemFilename = "N/A";
    }
    preallocatedModel = NSNIC(inferenceComponents);
    bestScenario = preallocatedModel.getScenario(true);

    import std.path : buildPath;

    csvTable = File(basename ~ ".csv", "w");
    csvPlots = File(basename ~ "_curves.csv", "w");

    auto header = ["id", "d_omega", "d_visual", "rounds", "evaluations", "failed evals", "duration"];
	header ~= "inf. n";
    header ~= iota(inferenceComponents).map!(i => "inf. t" ~ i.to!string).array;
	header ~= iota(inferenceComponents).map!(i => "inf. M" ~ i.to!string).array;
	header ~= iota(inferenceComponents).map!(i => "inf. s" ~ i.to!string).array;
	header ~= "inf. N_ref";
    if(simulating) {
        header ~= "sim. n";
        header ~= iota(simulationComponents).map!(i => "sim. t" ~ i.to!string).array;
        header ~= iota(simulationComponents).map!(i => "sim. M" ~ i.to!string).array;
        header ~= iota(simulationComponents).map!(i => "sim. s" ~ i.to!string).array;
        header ~= "sim. N_ref";
    }
	header ~= "source filename";
    
	csvTable.writeln("SEP=" ~ csvSeparator);
	csvTable.writeln(header.joiner(csvSeparator));
	csvTable.flush();

    return cast(int)(header.length);
}

int load_smc_problem_from_output(char* filename, double mu) {	
    resetProblemGlobals(filename.to!string);
    return smc_problem_outro(currentProblemFilename, mu, true);
}
int load_smc_problem_from_mscommand(char* filename, double mu) {
    resetProblemGlobals(filename.to!string);

	import ms : runPSMCAnalysis;
	auto msCommand = File(currentProblemFilename, "r").readln().splitter.array;
    auto psmc_filename = runPSMCAnalysis(msCommand, basename ~ "_" ~ (problemCount + 1).to!string);
    return smc_problem_outro(psmc_filename, mu, true);
}
int load_smc_problem_from_description(char* filename, double mu, int sequences, int sites) {
    description_problem_intro(filename.to!string);
    return smc_problem_body_and_outro(mu, sequences, sites);
}
int load_t2_problem_from_mscommand(char* filename, double referenceSize) {	
    resetProblemGlobals(filename.to!string);

	import ms : generateTsimIICR;
	auto msCommand = File(currentProblemFilename, "r").readln().splitter.array;
    generateTsimIICR(msCommand, data_cutoff_low, data_cutoff_high, 
        plotIntervals, referenceSize, plotTime, plotIICR, basename ~ "_" ~ (problemCount + 1).to!string);

    constrainData(plotTime, plotIICR, data_weight_low, data_weight_high, 
        fitTime, fitIICR);

    preallocatedModel.allocateDistanceComputation(fitTime);
    check_for_model();
    stopwatch.start();
	return 0;
}
int load_t2_problem_from_description(char* filename, int simulations) {
    description_problem_intro(filename.to!string);    
    return t2_problem_outro(simulations);
}
int load_exact_problem_from_description(char* filename) {
    description_problem_intro(filename.to!string);
    return exact_problem_outro();
}
int load_exact_simulated_problem(int n, double* times_p, double* migrations_p, double* sizes_p, double n0) {	
    simulated_problem_intro(n, times_p, migrations_p, sizes_p, n0);
    return exact_problem_outro();
}
int load_t2_simulated_problem(int n, double* times_p, double* migrations_p, double* sizes_p, double n0, 
    int simulations) {	

    simulated_problem_intro(n, times_p, migrations_p, sizes_p, n0);
    return t2_problem_outro(simulations);
}
int load_smc_simulated_problem(int n, double* times_p, double* migrations_p, double* sizes_p, double n0, 
    double mu, int sequences, int sites) {	
        
    simulated_problem_intro(n, times_p, migrations_p, sizes_p, n0);
    return smc_problem_body_and_outro(mu, sequences, sites);
}

int begin_round() {
    ++roundCount;
    if(roundCount % 10 == 0) {
        writefln(" >>> %s; %.12f; %s; [%(%.12f, %)]; [%(%.12f, %)];", roundCount, bestDistance, 
        bestScenario.islands, bestScenario.migrationRates, bestScenario.eventTimes);
    }
    return 0;
}

int open_problem() {
    ++problemCount;
    resetProblemGlobals();
    return 0;
}

@nogc
double evaluate(double* parameters) {
    ++evaluationCount;

    import std.math : isNaN;
	if(parameters[0].isNaN) {
		++failedEvaluationCount;
		return failDistance;
	}
    
    preallocatedModel.rebuild(parameters, inferDemeSizes, inferScale);
    auto distance = preallocatedModel.distance(fitIICR, distanceType, distanceParameter, sourceModel);

	if(!(distance < failDistance)) {
		++failedEvaluationCount;
		return failDistance;
	}
    
    if(evaluationCount == 1 || distance < bestDistance) {
        bestDistance = distance;
        preallocatedModel.copyScenario(bestScenario, true);
    }

    return distance;
}

double close_problem() {	
    stopwatch.stop();

    auto inferredModel = NSNIC(bestScenario);
    auto temp_iicr = new double[] (plotTime.length);
    auto temp_cdf = new double[] (plotTime.length);
    auto temp_pdf = new double[] (plotTime.length);
    
    logCurve("time", plotTime);

    getApproximateCDF(plotTime, plotIICR, temp_cdf);
    getApproximatePDF(plotTime, plotIICR, temp_pdf);
    logCurve("source-iicr", plotIICR);
    logCurve("source-cdf", temp_cdf);
    logCurve("source-pdf", temp_pdf);
    
    inferredModel.iicr(plotTime, temp_iicr);
    logCurve("inferred-iicr", temp_iicr);
    logCurve("inferred-cdf", plotTime.map!(t => inferredModel.cdf(t)));
    logCurve("inferred-pdf", plotTime.map!(t => inferredModel.pdf(t)));

    if(sourceModelIsKnown){
        sourceModel.iicr(plotTime, temp_iicr);
        logCurve("model-iicr", temp_iicr);
        logCurve("model-cdf", plotTime.map!(t => sourceModel.cdf(t)));
        logCurve("model-pdf", plotTime.map!(t => sourceModel.pdf(t)));
    }
    else {
        logCurve("model-iicr");
        logCurve("model-cdf");
        logCurve("model-pdf");
    }

    csvPlots.flush();

    inferredModel.allocateDistanceComputation(fitTime, true);
    immutable visualDistance = inferredModel.distance(
        fitIICR, 
        ErrorStrategy.Visual, 
        double.nan, 
        inferredModel
    );

	auto line = [
        problemCount.to!string, 
        (distanceType == ErrorStrategy.Visual) ? "N/A" : bestDistance.to!string, 
        visualDistance.to!string, 
        roundCount.to!string, 
        evaluationCount.to!string, 
        failedEvaluationCount.to!string,
        (stopwatch.peek.total!"msecs" / 1e3).to!string
    ];

    line ~= bestScenario.islands.to!string;
    line ~= bestScenario.eventTimes.map!(to!string).array;
    line ~= bestScenario.migrationRates.map!(to!string).array;
    line ~= bestScenario.demeSizes.map!(to!string).array;
    line ~= bestScenario.referenceSize.to!string;

    if(simulating) {
        auto sourceScenario = sourceModel.getScenario(true);
        line ~= sourceScenario.islands.to!string;
        line ~= sourceScenario.eventTimes.map!(to!string).array;
        line ~= sourceScenario.migrationRates.map!(to!string).array;
        line ~= sourceScenario.demeSizes.map!(to!string).array;
        line ~= sourceScenario.referenceSize.to!string;
    }    
	line ~= currentProblemFilename;

	csvTable.writeln(line.joiner(csvSeparator));
	csvTable.flush();
    
    stopwatch.reset();
    return bestDistance;
}

int finalize() {	
	Runtime.terminate();
	return 0;
}