module ms;

import std.stdio;
import std.process;
import std.range;
import std.range.interfaces;
import std.algorithm;
import std.conv;
import std.string;
import std.math;
import std.exception: enforce;

import numeric : GetLogTimeVector;
import nsnic : Scenario;

///
// this is the French PSMC ratio
enum double rhoThetaRatio = 0.14; // 1.0 / 5.0;

/// builds ms command for simulating T2 samples
auto getMSCommand(int n, double[] times, double[] migrs, double[] sizes, 
	int simulations, bool sameDeme = true, bool useSCRM = false) {
	
	auto msCommand = useSCRM ? ["./scrm"] : ["./ms"];
	msCommand ~= ["2", to!string(simulations), "-T", "-L", "-I", n.to!string];
	
	if(sameDeme)
		msCommand ~= ["2"] ~ repeat("0", n - 1).array ~ [migrs[0].to!string];
	else
		msCommand ~= ["1", "1"] ~ repeat("0", n - 2).array ~ [migrs[0].to!string];

	assert(times[0] == 0.0);
	
	if(sizes[0] != 1.0)
		msCommand ~= ["-eN", "0", sizes[0].to!string];

	foreach(i, t; times.drop(1)) {
		msCommand ~= ["-eM", (0.5 * t).to!string, migrs[i + 1].to!string];
		msCommand ~= ["-eN", (0.5 * t).to!string, sizes[i + 1].to!string];
	}

	return msCommand;
}

/// builds ms command for simulating sequences
auto getMSCommand(int n, double[] times, double[] migrs, double[] sizes, double referenceSize,
	double mu, int sequences, int sites, bool sameDeme = true, bool useSCRM = true) {
	
	immutable theta = sites * 4.0 * referenceSize * mu;
	immutable rho = theta * rhoThetaRatio;
	immutable display_digits = 8;

	auto msCommand = useSCRM ? ["./scrm"] : ["./ms"];
	msCommand ~= ["2", to!string(sequences), "-t", theta.to!string, "-r", rho.to!string, 
				  sites.to!string, "-p", display_digits.to!string, "-I", n.to!string];
	
	if(sameDeme)
		msCommand ~= ["2"] ~ repeat("0", n - 1).array ~ [migrs[0].to!string];
	else
		msCommand ~= ["1", "1"] ~ repeat("0", n - 2).array ~ [migrs[0].to!string];

	assert(times[0] == 0.0);
	
	if(sizes[0] != 1.0)
		msCommand ~= ["-eN", "0", sizes[0].to!string];

	foreach(i, t; times.drop(1)) {
		msCommand ~= ["-eM", (0.5 * t).to!string, migrs[i + 1].to!string];
		msCommand ~= ["-eN", (0.5 * t).to!string, sizes[i + 1].to!string];
	}

	return msCommand;
}

///
void generateTsimIICR(string[] msCommand, double x_low, double x_high, int x_intervals, 
	double referenceSize, ref double[] plotTime, ref double[] plotIICR, string basename = null) {

	auto t2_samples = getT2Samples(msCommand, basename);
	generateTsimIICR(t2_samples, x_low, x_high, x_intervals, referenceSize, plotTime, plotIICR);
	cropSimulatedIICR(plotTime, plotIICR);
}

///
auto runPSMCAnalysis(string[] msCommand, string basename) {
	immutable ms_output_fn = basename ~ ".msout";
	immutable psmc_input_fn = basename ~ ".psmcfa";
	immutable psmc_output_fn = basename ~ ".psmc";	

	// run ms command and save its output
	write("simulating sequences...");
	auto ms = pipeProcess(msCommand, Redirect.stdout);
	scope(exit) 
		wait(ms.pid);

	auto ms_output = File(ms_output_fn, "w");
	foreach(line; ms.stdout.byLine) {
		ms_output.writeln(line);
	}
	ms_output.close();
	
	// generate psmc input file
	auto pythonCommand = ["python", "./ms2psmcfa.py", ms_output_fn];
	write(" done.\npreparing PSMC input...");
	auto py = pipeProcess(pythonCommand, Redirect.stdout);
	scope(exit) 
		wait(py.pid);

	auto psmc_input = File(psmc_input_fn, "w");
	foreach(line; py.stdout.byLine) {
		psmc_input.writeln(line);
	}
	psmc_input.close();

	// run psmc
	auto psmcCommand = ["./psmc", "-N25", "-t15", "-r10", "-p", "4+25*2+4+6", "-o", psmc_output_fn, psmc_input_fn];
	write(" done.\nrunning PSMC analysis...");
	auto psmc = pipeProcess(psmcCommand, Redirect.stdout);
	wait(psmc.pid);

	write(" done.\n");
	return psmc_output_fn;
}

///
void parsePSMCFile(string psmc_filename, inout double mu, ref double[] samples, ref double[] lambdas) {
    auto psmc_file = File(psmc_filename, "r");
	auto lines = psmc_file.byLineCopy.map!(line => line.chomp).array;
	auto block = lines.splitter("//").tail(2).front;
	psmc_file.close();

	double theta, rho;
	samples = new double[] (0);
	lambdas = new double[] (0);
	foreach(line; block) {
		if(line.startsWith("TR")) {
			auto words = line.splitter('\t').drop(1);
			theta = words.front.to!double;
			
			words.popFront;			
			rho = words.front.to!double;
		}
		else if(line.startsWith("RS")) {
			auto words = line.splitter('\t').drop(2);
			samples ~= words.front.to!double;
			
			words.popFront;			
			lambdas ~= words.front.to!double;
		}
	}

	// scaling the PSMC
	immutable s = 100; // bin size
	immutable N0 = theta / (4 * mu * s);
	samples[] *= 2 * N0;
	lambdas[] *= N0;
}

///
void parseIICRFile(string filename, ref double[] time, ref double[] iicr) {
    auto file = File(filename, "r");
	auto lines = file.byLineCopy.map!(line => line.chomp).array;
	file.close();

	time = lines[0].splitter('\t').map!(to!double).array;
	iicr = lines[1].splitter('\t').map!(to!double).array;
	enforce(time.length == iicr.length, "Invalid IICR file: arrays length mismatch.");
}

private:

void generateTsimIICR(double[] samples, double x_low, double x_high, int intervals, 
	double referenceSize, ref double[] timeVector, ref double[] iicr) {

	immutable simulations = samples.length;
	immutable time_scale = (referenceSize == 0.0) ? 1.0 : 2 * referenceSize;
	immutable iicr_scale = (referenceSize == 0.0) ? 1.0 : referenceSize;
	
	// generate histograms
	timeVector = GetLogTimeVector(x_low / time_scale, x_high / time_scale, intervals);
	auto histogram = getHistogram(assumeSorted(timeVector), samples);

	auto shiftedTimeVector = new double[intervals + 1];
	shiftedTimeVector[0] = 0.0;
	shiftedTimeVector[1..$] = (timeVector[1..$] + timeVector[0..$-1]) * 0.5;
	auto shiftedHistogram = getHistogram(assumeSorted(shiftedTimeVector), samples);
	
	// generate F (distribution), f (density) and the IICR	
	auto normalizedHist = new double[intervals];
	normalizedHist[] = histogram[] / (simulations);
	auto distrib = [0.0] ~ normalizedHist.cumulativeFold!"a + b".array;
	
	auto normalizedShiftHist = new double[intervals];
	normalizedShiftHist[] = shiftedHistogram[] / (simulations);
	auto density = new double[intervals + 1];
	density[0..$-1] =  normalizedShiftHist[0..$] / (shiftedTimeVector[1..$] - shiftedTimeVector[0..$-1]);
	density[$ - 1] = 0.0;
	
	iicr = new double[intervals + 1];
	iicr[] = iicr_scale * (1.0 - distrib[]) / density[];
	timeVector[] *= time_scale;
}

void cropSimulatedIICR(ref double[] plotTime, ref double[] plotIICR) {
	immutable len = plotTime.length;
    assert(len > 0);
	assert(plotIICR.length == len);
	
	// int frontKeep = cast(int) plotIICR.length;
    // foreach(i, x; plotIICR) {
    //     if(!(x > 0 && x < double.max)) {
    //         frontKeep = cast(int) i;
    //         break;
    //     }
    // }

    // plotTime = plotTime[0 .. frontKeep - 1];
    // plotIICR = plotIICR[0 .. frontKeep - 1];

	auto frontDrop = 0;
    foreach(x; plotIICR) {
        if(!(x > 0 && x < double.max))
			++frontDrop;
		else
			break;
    }

	auto backDrop = 0;
    foreach_reverse(x; plotIICR) {
        if(!(x > 0 && x < double.max))
			++backDrop;
		else
			break;
    }

	assert(frontDrop + backDrop < len);
    plotTime = plotTime[frontDrop .. $ - backDrop];
    plotIICR = plotIICR[frontDrop .. $ - backDrop];
	assert(plotIICR.length + frontDrop + backDrop == len);
}

auto getT2Samples(string[] msCommand, string basename = null) {
	write("simulating T2 samples...");
	auto ms = pipeProcess(msCommand, Redirect.stdout);
	scope(exit) 
		wait(ms.pid);

	File ms_output;
	if(basename) {
		auto ms_output_fn = basename ~ ".msout";
		ms_output = File(ms_output_fn, "w");
		foreach(line; ms.stdout.byLine) {
			ms_output.writeln(line);
		}
		ms_output.reopen(ms_output_fn, "r");
	}
	else {
		ms_output = ms.stdout;
	}
	
	auto chunks = ms_output.byLine.drop(3).chunks(4);
	auto samples = chunks.map!(chunk => chunk.drop(2).front.splitter('\t').drop(1).front.to!double * 2).array;

	write(" done.\n");
	return samples;
}

auto getHistogram(R1, R2)(R1 bins, R2 data) /*if(isInputRange!R)*/ {
	assert(bins.length > 1);
	auto hist = new double[] (cast(int) bins.length - 1);
	
	hist[] = 0.0;
	foreach(immutable x; data) {
		immutable i = countUntil!"a > b"(bins, x);
		if(i < 1 || i > hist.length) {
			continue;
		}
		
		++hist[i - 1];
	}

	return hist;
}