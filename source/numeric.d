import std.range;

///
@nogc
auto trapezoidCuadrature(T, R)(T x, R y) {
    auto k = cast(int) x.length;
    assert(k > 1);
    assert(y.length == k);

    auto result = 0.0;
    foreach(immutable i; iota!int(1, k))
        result += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]);

    return result;
}

///
auto GetLogTimeVector(inout double min, inout double max, inout int intervals) {
	import std.math : log, exp;

	auto x = new double[](intervals + 1);
	x[0] = 0;
	x[1] = min;
	
	for(int i = 1; i < intervals - 1; i++) 
		x[i + 1] = min * exp((i * log(max/min)) / (intervals - 1));
	
	x[intervals] = max;
	return x;
}

///
auto sameSign(T)(T x, T y) {
	return (x > 0 && y > 0) || (x < 0 && y < 0);
}