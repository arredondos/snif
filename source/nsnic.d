import std.stdio;
import std.conv : to;
import ldc.attributes : fastmath;

import mir.ndslice;
import mir.math.common;

alias matrix = Slice!(Contiguous, [2], double*);
alias vector = Slice!(Contiguous, [1], double*);

///
enum ErrorStrategy {
    Visual,
    ExactPDF,
    ApproximatePDF
}

///
struct Scenario {
	int islands;
	double[] eventTimes;
	double[] migrationRates;
	double[] demeSizes;
    double referenceSize;
}

///
struct NSNIC {
    ///
    this(int n, double[] times, double[] migrations, double[] sizes, 
            double N0 = 0, bool report = false, bool sameDeme = true) {

        import std.range : iota;

        this.n = n;
        this.sizes = sizes;
        this.migrations = migrations;
        this.N0 = N0;
        
        if(N0 == 0) {
            this.times = times;
        } else {
            auto unscaled_times = new double[] (times.length);
            unscaled_times[] = times[] / (2 * N0);
            this.times = unscaled_times;
        }

        c = cast(int) times.length;
        assert(c == cast(int) migrations.length);
        assert(c == cast(int) sizes.length);
        assert(c > 0);
        
		qMatrices = new matrix[](c);
        qMatrices[0] = slice!double([3, 3], 0.0);
        createNIslandQMatrix(0, qMatrices[0]);

   		cumulativeExps = new matrix[](c);
        cumulativeExps[0] = slice!double([3, 3], 0.0);
        cumulativeExps[0].diagonal[] = 1.0;

        foreach(immutable i; iota!int(1, c)) {
            qMatrices[i] = slice!double([3, 3], 0.0);
            createNIslandQMatrix(i, qMatrices[i]);

			immutable dt = this.times[i] - this.times[i - 1];
            auto expm = slice!double([3, 3], 0.0);
            createNIslandExpMatrix(i - 1, dt, expm);
            cumulativeExps[i] = slice!double([3, 3], 0.0);
            mtimes33(cumulativeExps[i - 1], expm, cumulativeExps[i]);
			//cumulativeExps[i] = cumulativeExps[i - 1].mtimes(expm);
        }

        // writeln("qmat 0");
        // writeln(qMatrices[0]);

        // writeln("cumul exps 0");
        // writeln(cumulativeExps[0]);

        // foreach(immutable i; iota!int(1, c)) {
        //     writeln("\n\nqmat ", i);
        //     writeln(qMatrices[i]);

        //     writeln("cumul exps ", i);
        //     writeln(cumulativeExps[i]);
        // }
        
    }

    ///
    void allocateDistanceComputation(inout double[] times, bool computeIICR = false) {
        immutable k = cast(int) times.length;
        distanceArg = new double[] (k);
        
        fitTime = new double[] (k);
        fitTime[] = times[];

        fitIICR = new double[] (k);
        if(computeIICR)
            iicr(fitTime, fitIICR);
    }

    ///
    @fastmath @nogc
    void rebuild(double* parameters, bool size_changes, bool scaled) {
        import std.algorithm.iteration : map;
        import std.algorithm.mutation : copy;
        import std.algorithm.sorting : sort;
        import std.range : iota;
        import std.math : round;
        
        n = cast(int) parameters[0].round;        
        times[0] = 0.0;
        copy(parameters[1..c].map!(x => 10.0^^x), times[1..$]);
        times.sort();

        copy(parameters[c..2*c], migrations);
        if(size_changes)
            copy(parameters[2*c..3*c], sizes); 
        else
            sizes[] = 1.0;
        
        if(scaled) {
            immutable last_index = size_changes ? 3 * c : 2 * c;
            N0 = parameters[last_index];
            times[] /= 2 * N0;
        }
        else {
            N0 = 0.0;
        }

        createNIslandQMatrix(0, qMatrices[0]);
        cumulativeExps[0][] = 0.0;
        cumulativeExps[0].diagonal[] = 1.0;
        foreach(immutable i; iota!int(1, c)) {
            createNIslandQMatrix(i, qMatrices[i]);
			immutable dt = times[i] - times[i - 1];
            createNIslandExpMatrix(i - 1, dt, currentExp);
            mtimes33(cumulativeExps[i - 1], currentExp, cumulativeExps[i]);
        }
        iicr(fitTime, fitIICR);
    }
    
    ///
    this(string jsonFilename) {
        import std.range : array;
        import std.file : readText;
        import std.algorithm.iteration : map;
        import std.json : JSONType, JSONValue, parseJSON;

        auto jstring = readText(jsonFilename);
        auto jtree = parseJSON(jstring)["PSNIC"];
        
        auto parseDouble(ref JSONValue val) {
            double result;
            switch(val.type) {
                case JSONType.integer:
                    result = val.integer.to!double;
                    break;
                case JSONType.uinteger:
                    result = val.uinteger.to!double;
                    break;
                case JSONType.float_:
                    result = val.floating;
                    break;
                default:
                    assert(0); 
            }
            return result;            
        }
        
        auto n = jtree["n"].integer.to!int;
        auto t = jtree["t"].arrayNoRef.map!parseDouble.array;
        auto m = jtree["M"].arrayNoRef.map!parseDouble.array;
        auto s = jtree["s"].arrayNoRef.map!parseDouble.array;
        auto n0 = ("Nref" in jtree) ? parseDouble(jtree["Nref"]) : 0.0;
        
        this(n, t, m, s, n0);
    }

    ///
    this(Scenario scenario) {
        this(scenario.islands, 
             scenario.eventTimes, 
             scenario.migrationRates, 
             scenario.demeSizes, 
             scenario.referenceSize);
    }
    
    /// this constructor just pre-allocates the space
    this(int components) {
        c = components;
        n = int.init;
        times = new double[] (c);
        migrations = new double[] (c);
        sizes = new double[] (c);
        N0 = double.init;
        
		qMatrices = new matrix[](c);
   		cumulativeExps = new matrix[](c);
        currentExp = slice!double([3, 3], 0.0); 
        foreach(i; iota(c)) {
            qMatrices[i] = slice!double([3, 3], 0.0);   
            cumulativeExps[i] = slice!double([3, 3], 0.0);
        }
    }

    ///
	@fastmath
	auto cdf(inout double t, int row = 0) {
        auto st = (N0 == 0) ? t : t / (2 * N0);
		auto j = getComponentIndex(st); 
		auto hvec = cumulativeExps[j][row, 0 .. $];
		// auto vvec = createNIslandExpColumn(j, st - times[j]);
		// return hvec.mtimes(vvec);


        // immutable hvec = cumulativeExps[j][row, 0 .. $];

        double f11, f12, f21, f22;
        computeExpMatCoefs(j, st - times[j], f11, f12, f21, f22);

        immutable exp1 = 1.0 - f11 - f12;
        immutable exp2 = 1.0 - f21 - f22;
        return exp1 * hvec[0] + exp2 * hvec[1] + hvec[2];
	}

    ///
	@fastmath @nogc
	auto pdf(inout double t, int row = 0) {
        immutable st = (N0 == 0) ? t : t / (2 * N0);
        immutable j = getComponentIndex(st);
        immutable dt = st - times[j];

        double pdf, cdf;
        iicr_internal(j, dt, row, cdf, pdf);
        return (N0 == 0) ? pdf : pdf / N0;
        // auto st = (N0 == 0) ? t : t / (2 * N0);
		// auto j = getComponentIndex(st); 

		// auto hvec = cumulativeExps[j][row, 0 .. $];
		// auto vvec = createNIslandExpColumn(j, st - times[j]);
        // auto result = hvec.mtimes(qMatrices[j]).mtimes(vvec);
        // return (N0 == 0) ? result : result / N0;
	}

    ///
    @fastmath @nogc
    void iicr(inout double[] domain, ref double[] result, inout int row = 0) {
        import numeric : sameSign;
        import std.math : isNaN;

        immutable len = cast(int) domain.length;
        assert(result.length == len);
        
        // debug writeln("\ndomain = ", domain);

        auto i = 0;
        while(i < len) {
            double st = (N0 == 0) ? domain[i] : domain[i] / (2 * N0);
            immutable j = getComponentIndex(st);
            immutable limit = computeLimitIICR(j);

            // debug writeln("j = ", j);
            // debug writeln("st = ", st);
            // debug writeln("limit = ", limit);
            // debug writeln("event times = ", times);
            
            auto unstable = false;
            while((j == c - 1) || (st < times[j + 1])) {
                if(unstable) {
                    // debug writeln("unstable :(");
                    result[i] = limit;
                }
                else {
                    // immutable vvec = createNIslandExpColumn(j, st - times[j]);
                    // immutable cdf0  = hvec.mtimes(vvec);
                    // immutable pdf0  = hvec.mtimes(qMatrices[j]).mtimes(vvec);
                    // immutable iicr0 = ((N0 == 0.0) ? 1.0 : N0) * (1.0 - cdf0) / pdf0;                    

                    double cdf, pdf;
                    iicr_internal(j, st - times[j], row, cdf, pdf);
                    //debug writeln("states:\t", states);
                    immutable iicr = (N0 == 0.0 ? 1.0 : N0) * (1.0 - cdf) / pdf;
                    // debug writeln("pdf = ", pdf);
                    // debug writeln("cdf = ", cdf);
                    // debug writeln("iicr = ", iicr);
                    unstable =  isNaN(iicr)
                            //  || 1.0 - fabs(iicr/limit) < 1e-5
                            //  || fabs(cdf - 1.0) < 1e-10
                             || fabs(pdf) < 1e-14;
                    
                    result[i] = unstable ? limit : iicr;
                }
                
                ++i;
                if(!(i < len)) 
                    break;

                st = (N0 == 0) ? domain[i] : domain[i] / (2 * N0);
            }
        }

        //debug writeln("done");
    }

    ///
	@fastmath @nogc
	double distance(inout double[] iicr, ErrorStrategy errStrategy, 
                    double exponent = double.nan, NSNIC description = NSNIC.init) {

        import std.math : fabs;

		immutable k = fitTime.length;
		assert(iicr.length == k);

        final switch(errStrategy) {
            case ErrorStrategy.Visual:
                auto result = 0.0;
                foreach(immutable i; iota!int(k))
                    result += fabs(fitIICR[i] - iicr[i]);

                return result;

            case ErrorStrategy.ExactPDF:
                import numeric : trapezoidCuadrature;
                foreach(immutable i; iota!int(k))
                    distanceArg[i] = fabs(fitIICR[i] - iicr[i]) * description.pdf(fitTime[i]);

                return trapezoidCuadrature(fitTime, distanceArg);

            case ErrorStrategy.ApproximatePDF:
                import numeric : trapezoidCuadrature;
                getApproximatePDF(fitTime, iicr, distanceArg);
                foreach(immutable i; iota!int(k))
                    distanceArg[i] = fabs(fitIICR[i] - iicr[i]) * (distanceArg[i] ^^ exponent);

                return trapezoidCuadrature(fitTime, distanceArg);
        }
        assert(0);
	}

    ///
    @nogc
    void copyScenario(ref Scenario target, bool timesInGenerations) {
        import std.algorithm.mutation : copy;
        target.islands = n;
        target.referenceSize = N0;
        copy(times, target.eventTimes);
        copy(migrations, target.migrationRates);
        copy(sizes, target.demeSizes);
        if (timesInGenerations && !(N0 == 0))
            target.eventTimes[] *= 2 * N0;
    }

    ///
    auto getScenario(bool timesInGenerations) {
        auto scaled_times = times.dup;
        if (timesInGenerations && !(N0 == 0))
            scaled_times[] *= 2 * N0;

        return Scenario(n, scaled_times, migrations.dup, sizes.dup, N0);
    }

    ///
    auto timeScale() {
        return (N0 == 0.0) ? 1.0 : 2 * N0;
    }

private:
	int n, c;
    double N0;
    double[] times, migrations, sizes;
    double[] fitTime, fitIICR, distanceArg;
	matrix[] cumulativeExps, qMatrices;
    matrix currentExp;

	@fastmath @nogc
	auto getComponentIndex(inout double t) {
		import std.algorithm.searching : countUntil;
		import std.range : assumeSorted;

		assert(t >= 0.0);
		immutable j = cast(int) countUntil!"a > b"(assumeSorted(times), t);
		return (j == -1) ? (c - 1) : (j - 1);
	}

    @fastmath @nogc
    void createNIslandQMatrix(inout int i, ref matrix qmatrix) {
        immutable m = migrations[i];
        immutable z = 1.0 / sizes[i];

        qmatrix[0, 0] = -m - z;
        qmatrix[0, 1] = m;
        qmatrix[0, 2] = z;
        qmatrix[1, 0] = m / (n - 1);
        qmatrix[1, 1] = -m / (n - 1);
        
        qmatrix[1, 2] = 0.0;
        qmatrix[2, 0] = 0.0;
        qmatrix[2, 1] = 0.0;
        qmatrix[2, 2] = 0.0;
    }

    @fastmath @nogc
    void createNIslandExpMatrix(inout int i, inout double t, ref matrix expm) {
        double f11, f12, f21, f22;
        computeExpMatCoefs(i, t, f11, f12, f21, f22);
        expm[0, 0] = f11;
        expm[0, 1] = f12;
        expm[0, 2] = 1.0 - (f11 + f12);
        expm[1, 0] = f21;
        expm[1, 1] = f22;
        expm[1, 2] = 1.0 - (f21 + f22);
        expm[2, 0] = 0.0;
        expm[2, 1] = 0.0;
        expm[2, 2] = 1.0;
    }

    @fastmath
    auto createNIslandExpColumn(inout int i, inout double t) {
        double f11, f12, f21, f22;
        computeExpMatCoefs(i, t, f11, f12, f21, f22);
        vector col = [1.0 - (f11 + f12), 
                      1.0 - (f21 + f22), 
                      1.0].sliced;

        return col; 
    }

    @fastmath @nogc
    void computeExpMatCoefs(inout int i, inout double time, 
            out double f11, out double f12, out double f21, out double f22) {

        immutable m = migrations[i] * sizes[i];
        immutable t = time / sizes[i];
        immutable alpha = m * n + n - 1;
        immutable delta = sqrt(alpha * alpha - 4 * m * (n - 1));
        immutable e1 = exp(0.5 * t * ( delta - alpha) / (n - 1));
        immutable e2 = exp(0.5 * t * (-delta - alpha) / (n - 1));

        f21 = m * (e1 - e2) / delta;
        f12 = (n - 1) * f21;
        f11 = (e2 * (delta + alpha - 2 * m) + e1 * (delta - alpha + 2 * m)) / (2 * delta);
        f22 = (e1 * (delta + alpha - 2 * m) + e2 * (delta - alpha + 2 * m)) / (2 * delta);
    }

    @fastmath @nogc
    auto iicr_internal(inout int j, inout double dt, inout int row, out double cdf, out double pdf) {
        immutable hvec = cumulativeExps[j][row, 0 .. $];

        double f11, f12, f21, f22;
        computeExpMatCoefs(j, dt, f11, f12, f21, f22);

        immutable exp1 = 1.0 - f11 - f12;
        immutable exp2 = 1.0 - f21 - f22;
        cdf = exp1 * hvec[0] + exp2 * hvec[1] + hvec[2];
        immutable double[3] temp = [
            hvec[0]*qMatrices[j][0,0] + hvec[1]*qMatrices[j][1,0] + hvec[2]*qMatrices[j][2,0],
            hvec[0]*qMatrices[j][0,1] + hvec[1]*qMatrices[j][1,1] + hvec[2]*qMatrices[j][2,1],
            hvec[0]*qMatrices[j][0,2] + hvec[1]*qMatrices[j][1,2] + hvec[2]*qMatrices[j][2,2],
        ];
        pdf = temp[0] * exp1 + temp[1] * exp2 + temp[2];

        //double[3] states = [f11 * hvec[0] + f21 * hvec[1], f12 * hvec[0] + f22 * hvec[1], cdf];
        //return states;
    }

    @fastmath @nogc
    double computeLimitIICR(inout int i) {
        assert(i < migrations.length);

        immutable gamma = migrations[i] * sizes[i] / (n - 1);
        immutable delta = sqrt((1 + n * gamma)^^2 - 4 * gamma);
        immutable alpha = 0.5 * (1 + n * gamma + delta);
        immutable beta  = 0.5 * (1 + n * gamma - delta);

        immutable limit = sizes[i] * (alpha - 1.0) / (gamma - beta);
        return (N0 == 0) ? limit : N0 * limit;
    }

    @fastmath @nogc
    void mtimes33(ref matrix left, ref matrix right, ref matrix result) {
        result[0, 2] = 1.0;
        result[0, 2] -= (result[0, 0] = left[0,0] * right[0,0] + left[0,1] * right[1,0]);
        result[0, 2] -= (result[0, 1] = left[0,0] * right[0,1] + left[0,1] * right[1,1]);

        result[1, 2] = 1.0;
        result[1, 2] -= (result[1, 0] = left[1,0] * right[0,0] + left[1,1] * right[1,0]);
        result[1, 2] -= (result[1, 1] = left[1,0] * right[0,1] + left[1,1] * right[1,1]);

        result[2, 2] = 1.0;
        result[2, 0] = (result[2, 1] = 0.0);
    }
}

@fastmath @nogc
auto getExpIntMinusLambda(const double[] times, const double[] iicr, ref double[] result) {
    import std.algorithm.iteration : map;
    import numeric : trapezoidCuadrature;

    immutable len = cast(int) iicr.length;
    assert(times.length == len);
    assert(len > 0);

    // auto arg = new double[] (len);
    // arg[] = 1.0 / iicr[];

    result[0] = 1.0 - (1.0 / iicr[0]);
    for(int i = 1; i < len; i++)
        result[i] = exp(-trapezoidCuadrature(
            times[0 .. i + 1], 
            iicr[0 .. i + 1].map!(x => 1.0 / x)
        ));

    return result;
}

@fastmath @nogc
auto getApproximateCDF(inout double[] times, inout double[] iicr, ref double[] result) {
    getExpIntMinusLambda(times, iicr, result);
    result[] = 1.0 - result[];
    return result;
}

@fastmath @nogc
auto getApproximatePDF(inout double[] times, inout double[] iicr, ref double[] result) {
    getExpIntMinusLambda(times, iicr, result);
    result[] = result[] / iicr[];
    return result;
}
