/++
$(H2 Level 1)

$(SCRIPT inhibitQuickIndex = 1;)

This is a submodule of $(MREF mir,glas).

The Level 1 GLAS perform vector and vector-vector operations.

$(BOOKTABLE $(H2 Vector-vector operations),
$(T2 rot, apply Givens rotation)
$(T2 axpy, constant times a vector plus a vector)
$(T2 dot, dot product)
$(T2 dotc, dot product, conjugating the first vector)
)

$(BOOKTABLE $(H2 Vector operations),
$(TR $(TH Function Name) $(TH Description))
$(T2 nrm2, Euclidean norm)
$(T2 sqnrm2, square of Euclidean norm)
$(T2 asum, sum of absolute values)
$(T2 iamax, index of max abs value)
$(T2 amax, max abs value)
)

All functions except $(LREF iamax) work with multidimensional tensors.

GLAS does not provide `swap`, `scal`, and `copy`  functions.
This functionality is part of $(MREF_ALTTEXT ndslice, mir, ndslice) package. Examples can be found below.

License:   $(HTTP boost.org/LICENSE_1_0.txt, Boost License 1.0).
Copyright: Copyright © 2016-, Ilya Yaroshenko
Authors:   Ilya Yaroshenko

Macros:
T2=$(TR $(TDNW $(LREF $1)) $(TD $+))
SUBMODULE = $(MREF_ALTTEXT $1, mir, glas, $1)
SUBREF = $(REF_ALTTEXT $(TT $2), $2, mir, glas, $1)$(NBSP)
NDSLICEREF = $(REF_ALTTEXT $(TT $2), $2, mir, ndslice, $1)$(NBSP)
+/
module mir.glas.l1;

/// SWAP
unittest
{
    import std.algorithm.mutation: swap;
    import mir.ndslice.allocation: slice;
    import mir.ndslice.algorithm: each;
    import std.typecons: Yes;
    auto x = slice!double(4);
    auto y = slice!double(4);
    x[] = [0, 1, 2, 3];
    y[] = [4, 5, 6, 7];
    each!(swap)(x, y);
    assert(x == [4, 5, 6, 7]);
    assert(y == [0, 1, 2, 3]);
}

/// SCAL
unittest
{
    import mir.ndslice.allocation: slice;
    import std.typecons: Yes;
    auto x = slice!double(4);
    x[] = [0, 1, 2, 3];
    x[] *= 2.0;
    assert(x == [0, 2, 4, 6]);
}

/// COPY
unittest
{
    import mir.ndslice.allocation: slice;
    auto x = slice!double(4);
    auto y = slice!double(4);
    x[] = [0, 1, 2, 3];
    y[] = x;
    assert(y == [0, 1, 2, 3]);
}

import std.traits;
import std.meta;

import mir.math.common;
import mir.internal.utility;
import mir.ndslice.slice;
import mir.ndslice.algorithm : reduce, each;

version(LDC)
    import ldc.attributes : fastmath;
else
    enum fastmath;

@fastmath:

template _rot(alias c, alias s)
{
    @fastmath
    void _rot(X, Y)(ref X xr, ref Y yr)
    {
        auto x = xr;
        auto y = yr;
        auto t1 = c * x + s * y;
        static if (isComplex!(typeof(c)))
        {
            auto t2 = (c.re - c.im * 1fi) * y;
        }
        else
            auto t2 = c * y;
        static if (isComplex!(typeof(s)))
        {
            t2 -= (s.re - s.im * 1fi) * x;
        }
        else
            t2 -= s * x;
        xr = t1;
        yr = t2;
    }
}

template _axpy(alias a)
{
    @fastmath
    void _axpy(X, Y)(ref X x, ref Y y)
    {
        y += a * x;
    }
}

A _fmuladd(A, B, C)(A a, in B b, in C c)
{
    return a + b * c;
}

A _fmuladdc(A, B, C)(A a, in B b, in C c)
{
    static if (isComplex!B)
    {
        return a + (b.re - b.im * 1fi) * c;
    }
    else
        return a + b * c;
}

A _nrm2(A, B)(A a, in B b)
{
    static if (isComplex!B)
        return a + b.re * b.re + b.im * b.im;
    else
        return a + b * b;
}

A _asum(A, B)(A a, in B b)
{
    static if (isComplex!B)
    {
        return a + (b.re.fabs + b.im.fabs);
    }
    else
    static if (isFloatingPoint!B)
    {
        return a + b.fabs;
    }
    else
    {
        static if (isUnsigned!B)
            return a + b;
        else
            return a + (b >= 0 ? b : -b);
    }
}

A _amax(A, B)(A a, in B b)
{
    static if (isComplex!B)
    {
        return a.fmax(b.re.fabs + b.im.fabs);
    }
    else
    static if (isFloatingPoint!B)
    {
        return a.fmax(b.fabs);
    }
    else
    {
        static if (!isUnsigned!B)
            b = (b >= 0 ? b : -b);
        return a >= b ? a : b;
    }
}

private enum _shouldBeCastedToUnqual(T) = isPointer!T && !is(Unqual!T == T);

/++
Applies a plane rotation, where the  `c` (cos) and `s` (sin) are scalars.
Uses unrolled loops for strides equal to one.
Params:
    c = cos scalar
    s = sin scalar
    x = first n-dimensional tensor
    y = second n-dimensional tensor
BLAS: SROT, DROT, CROT, ZROT, CSROT, ZDROTF
+/
void rot(C, S, SliceKind kind1, SliceKind kind2, size_t[] packs, Iterator1, Iterator2)(in C c, in S s, Slice!(kind1, packs, Iterator1) x, Slice!(kind2, packs, Iterator2) y)
{
    assert(x.shape == y.shape, "constraints: x and y must have equal shapes");
    pragma(inline, false);
    each!(_rot!(c, s))(x, y);
}

///
unittest
{
    import mir.ndslice.allocation: slice;
    auto x = slice!double(4);
    auto y = slice!double(4);
    auto a = slice!double(4);
    auto b = slice!double(4);
    double cos = 3.0 / 5;
    double sin = 4.0 / 5;
    x[] = [0, 1, 2, 3];
    y[] = [4, 5, 6, 7];
    foreach (i; 0 .. 4)
    {
        a[i] = cos * x[i] + sin * y[i];
        b[i] = cos * y[i] - sin * x[i];
    }
    rot(cos, sin, x, y);
    assert(x == a);
    assert(y == b);
}

/++
Constant times a vector plus a vector.
Uses unrolled loops for strides equal to one.
Params:
    a = scale parameter
    x = first n-dimensional tensor
    y = second n-dimensional tensor
BLAS: SAXPY, DAXPY, CAXPY, ZAXPY
+/
void axpy(A, SliceKind kind1, SliceKind kind2, size_t[] packs, Iterator1, Iterator2)(in A a, Slice!(kind1, packs, Iterator1) x, Slice!(kind2, packs, Iterator2) y)
{
    static if (_shouldBeCastedToUnqual!Iterator2)
    {
        .axpy(a, cast(Slice!(N, Unqual!Iterator1))x, cast(Slice!(N, Unqual!Iterator2))y);
    }
    else
    {
        assert(x.shape == y.shape, "constraints: x and y must have equal shapes");
        pragma(inline, false);
        each!(_axpy!a)(x, y);
    }
}

/// SAXPY, DAXPY
unittest
{
    import mir.ndslice.allocation: slice;
    auto x = slice!double(4);
    auto y = slice!double(4);
    x[] = [0, 1, 2, 3];
    y[] = [4, 5, 6, 7];
    axpy(2.0, x, y);
    assert(y == [4, 7, 10, 13]);
}

/// SAXPY, DAXPY
unittest
{
    import mir.ndslice.allocation: slice;

    auto a = 3 + 4i;
    auto x = slice!cdouble(2);
    auto y = slice!cdouble(2);
    x[] = [0 + 1i, 2 + 3i];
    y[] = [4 + 5i, 6 + 7i];
    axpy(a, x, y);
    assert(y == [a * (0 + 1i) + (4 + 5i), a * (2 + 3i) + (6 + 7i)]);
}

/++
Forms the dot product of two vectors.
Uses unrolled loops for strides equal to one.
Returns: dot product `conj(xᐪ) × y`
Params:
    F = type for summation (optional template parameter)
    x = first n-dimensional tensor
    y = second n-dimensional tensor
BLAS: SDOT, DDOT, SDSDOT, DSDOT, CDOTC, ZDOTC
+/
F dot(F, SliceKind kind1, SliceKind kind2, size_t[] packs, Iterator1, Iterator2)(Slice!(kind1, packs, Iterator1) x, Slice!(kind2, packs, Iterator2) y)
{
    static if (allSatisfy!(_shouldBeCastedToUnqual, Iterator1, Iterator2))
    {
        return .dot!F(cast(Slice!(kind1, packs, Unqual!Iterator1))x, cast(Slice!(kind2, packs, Unqual!Iterator2))y);
    }
    else
    {
        assert(x.shape == y.shape, "constraints: x and y must have equal shapes");
        pragma(inline, false);
        return reduce!(_fmuladd)(cast(F)(0), x, y);
    }
}

/// SDOT, DDOT
unittest
{
    import mir.ndslice.allocation: slice;
    auto x = slice!double(4);
    auto y = slice!double(4);
    x[] = [0, 1, 2, 3];
    y[] = [4, 5, 6, 7];
    assert(dot(x, y) == 5 + 12 + 21);
}

/// ditto
auto dot(SliceKind kind1, SliceKind kind2, size_t[] packs, Iterator1, Iterator2)(Slice!(kind1, packs, Iterator1) x, Slice!(kind2, packs, Iterator2) y)
{
    return .dot!(Unqual!(typeof(x[0] * y[0])))(x, y);
}

/// SDOT, DDOT
unittest
{
    import mir.ndslice.allocation: slice;
    auto x = slice!double(4);
    auto y = slice!double(4);
    x[] = [0, 1, 2, 3];
    y[] = [4, 5, 6, 7];
    assert(dot(x, y) == 5 + 12 + 21);
}

/// SDSDOT, DSDOT
unittest
{
    import mir.ndslice.allocation: slice;
    auto x = slice!float(4);
    auto y = slice!float(4);
    x[] = [0, 1, 2, 3];
    y[] = [4, 5, 6, 7];
    assert(dot!real(x, y) == 5 + 12 + 21); // 80-bit FP for x86 CPUs
}

/// CDOTU, ZDOTU
unittest
{
    import mir.ndslice.allocation: slice;

    auto x = slice!cdouble(2);
    auto y = slice!cdouble(2);
    x[] = [0 + 1i, 2 + 3i];
    y[] = [4 + 5i, 6 + 7i];
    version(LDC) // DMD Internal error: backend/cgxmm.c 628
    assert(dot(x, y) == (0 + 1i) * (4 + 5i) + (2 + 3i) * (6 + 7i));
}

/++
Forms the dot product of two complex vectors.
Uses unrolled loops for strides equal to one.
Returns: dot product `xᐪ × y`
Params:
    F = type for summation (optional template parameter)
    x = first n-dimensional tensor
    y = second n-dimensional tensor
BLAS: CDOTU, ZDOTU
+/
F dotc(F, SliceKind kind1, SliceKind kind2, size_t[] packs, Iterator1, Iterator2)(Slice!(kind1, packs, Iterator1) x, Slice!(kind2, packs, Iterator2) y)
    if (isComplex!(DeepElementType!(typeof(x))) && isComplex!(DeepElementType!(typeof(y))))
{
    static if (allSatisfy!(_shouldBeCastedToUnqual, Iterator1, Iterator2))
    {
        return .dotc!F(cast(Slice!(N, Unqual!Iterator1))x, cast(Slice!(N, Unqual!Iterator2))y);
    }
    else
    {
        assert(x.shape == y.shape, "constraints: x and y must have equal shapes");
        pragma(inline, false);
        return reduce!(_fmuladdc)(cast(F)(0), x, y);
    }
}

/// ditto
auto dotc(SliceKind kind1, SliceKind kind2, size_t[] packs, Iterator1, Iterator2)(Slice!(kind1, packs, Iterator1) x, Slice!(kind2, packs, Iterator2) y)
{
    return .dotc!(Unqual!(typeof(x[x.shape.init] * y[y.shape.init])))(x, y);
}

/// CDOTC, ZDOTC
unittest
{
    import mir.ndslice.allocation: slice;

    auto x = slice!cdouble(2);
    auto y = slice!cdouble(2);
    x[] = [0 + 1i, 2 + 3i];
    y[] = [4 + 5i, 6 + 7i];
    version(LDC) // DMD Internal error: backend/cgxmm.c 628
    assert(dotc(x, y) == (0 + -1i) * (4 + 5i) + (2 + -3i) * (6 + 7i));
}

/++
Returns the euclidean norm of a vector.
Uses unrolled loops for stride equal to one.
Returns: euclidean norm `sqrt(conj(xᐪ)  × x)`
Params:
    F = type for summation (optional template parameter)
    x = n-dimensional tensor
BLAS: SNRM2, DNRM2, SCNRM2, DZNRM2
+/
F nrm2(F, SliceKind kind, size_t[] packs, Iterator)(Slice!(kind, packs, Iterator) x)
{
    static if (_shouldBeCastedToUnqual!Iterator)
        return .sqnrm2!F(cast(Slice!(N, Unqual!R))x).sqrt;
    else
        return .sqnrm2!F(x).sqrt;
}

/// ditto
auto nrm2(SliceKind kind, size_t[] packs, Iterator)(Slice!(kind, packs, Iterator) x)
{
    return .nrm2!(realType!(typeof(x[x.shape.init] * x[x.shape.init])))(x);
}

/// SNRM2, DNRM2
unittest
{
    import mir.ndslice.allocation: slice;
    import std.math: sqrt, approxEqual;
    auto x = slice!double(4);
    x[] = [0, 1, 2, 3];
    assert(nrm2(x).approxEqual(sqrt(1.0 + 4 + 9)));
}

/// SCNRM2, DZNRM2
unittest
{
    import mir.ndslice.allocation: slice;
    import std.math: sqrt, approxEqual;

    auto x = slice!cdouble(2);
    x[] = [0 + 1i, 2 + 3i];

    assert(nrm2(x).approxEqual(sqrt(1.0 + 4 + 9)));
}

/++
Forms the square of the euclidean norm.
Uses unrolled loops for stride equal to one.
Returns: `conj(xᐪ) × x`
Params:
    F = type for summation (optional template parameter)
    x = n-dimensional tensor
+/
F sqnrm2(F, SliceKind kind, size_t[] packs, Iterator)(Slice!(kind, packs, Iterator) x)
{
    static if (_shouldBeCastedToUnqual!Iterator)
    {
        return .sqnrm2!F(cast(Slice!(N, Unqual!R))x);
    }
    else
    {
        pragma(inline, false);
        return reduce!(_nrm2)(F(0), x);
    }
}

/// ditto
auto sqnrm2(SliceKind kind, size_t[] packs, Iterator)(Slice!(kind, packs, Iterator) x)
{
    return .sqnrm2!(realType!(typeof(x[x.shape.init] * x[x.shape.init])))(x);
}

///
unittest
{
    import mir.ndslice.allocation: slice;
    auto x = slice!double(4);
    x[] = [0, 1, 2, 3];
    assert(sqnrm2(x) == 1.0 + 4 + 9);
}

///
unittest
{
    import mir.ndslice.allocation: slice;

    auto x = slice!cdouble(2);
    x[] = [0 + 1i, 2 + 3i];

    assert(sqnrm2(x) == 1.0 + 4 + 9);
}

/++
Takes the sum of the `|Re(.)| + |Im(.)|`'s of a vector and
    returns a single precision result.
Returns: sum of the `|Re(.)| + |Im(.)|`'s
Params:
    F = type for summation (optional template parameter)
    x = n-dimensional tensor
BLAS: SASUM, DASUM, SCASUM, DZASUM
+/
F asum(F, SliceKind kind, size_t[] packs, Iterator)(Slice!(kind, packs, Iterator) x)
{
    static if (_shouldBeCastedToUnqual!Iterator)
    {
        return .asum!F(cast(Slice!(N, Unqual!R))x);
    }
    else
    {
        pragma(inline, false);
        return reduce!(_asum)(F(0), x);
    }
}

/// ditto
auto asum(SliceKind kind, size_t[] packs, Iterator)(Slice!(kind, packs, Iterator) x)
{
    alias T = DeepElementType!(typeof(x));
    return .asum!(realType!T)(x);
}

/// SASUM, DASUM
unittest
{
    import mir.ndslice.allocation: slice;
    auto x = slice!double(4);
    x[] = [0, -1, -2, 3];
    assert(asum(x) == 1 + 2 + 3);
}

/// SCASUM, DZASUM
unittest
{
    import mir.ndslice.allocation: slice;

    auto x = slice!cdouble(2);
    x[] = [0 - 1i, -2 + 3i];

    assert(asum(x) == 1 + 2 + 3);
}

/++
Finds the index of the first element having maximum `|Re(.)| + |Im(.)|`.
Return: index of the first element having maximum `|Re(.)| + |Im(.)|`
Params: x = 1-dimensional tensor
BLAS: ISAMAX, IDAMAX, ICAMAX, IZAMAX
+/
sizediff_t iamax(SliceKind kind, Iterator)(Slice!(kind, [1], Iterator) x)
{
    static if (_shouldBeCastedToUnqual!Iterator)
    {
        return .iamax(cast(Slice!(1, Unqual!R))x);
    }
    else
    {
        pragma(inline, false);
        if (x.length == 0)
            return -1;
        if (x._stride == 0)
            return 0;
        alias T = Unqual!(DeepElementType!(typeof(x)));
        alias F = realType!T;
        static if (isFloatingPoint!F)
            auto m = -double.infinity;
        else
            auto m = F.min;
        sizediff_t l = x.length;
        sizediff_t r = x.length;
        do
        {
            auto f = x.front;
            static if (isComplex!T)
            {
                auto e = f.re.fabs + f.im.fabs;
            }
            else
            static if (isFloatingPoint!T)
            {
                auto e = f.fabs;
            }
            else
            {
                static if (isUnsigned!T)
                    auto e = f;
                else
                    auto e = (f >= 0 ? f : -f);
            }

            if (e > m)
            {
                m = e;
                r = x.length;
            }
            x.popFront;
        }
        while (x.length);
        return l - r;
    }
}

/// ISAMAX, IDAMAX
unittest
{
    import mir.ndslice.allocation: slice;
    auto x = slice!double(6);
    //     0  1   2   3   4  5
    x[] = [0, -1, -2, -3, 3, 2];
    assert(iamax(x) == 3);
    // -1 for empty vectors
    assert(iamax(x[0 .. 0]) == -1);
}

/// ICAMAX, IZAMAX
unittest
{
    import mir.ndslice.allocation: slice;

    auto x = slice!cdouble(4);
    //        0          1          2         3
    x[] = [0 + -1i, -2 + 3i, 2 + 3i, 2 + 2i];

    assert(iamax(x) == 1);
    // -1 for empty vectors
    assert(iamax(x[$ .. $]) == -1);
}

/++
Takes the sum of the `|Re(.)| + |Im(.)|`'s of a vector and
    returns a single precision result.
Returns: sum of the `|Re(.)| + |Im(.)|`'s
Params:
    F = type for summation (optional template parameter)
    x = n-dimensional tensor
BLAS: SASUM, DASUM, SCASUM, DZASUM
+/
auto amax(SliceKind kind, size_t[] packs, Iterator)(Slice!(kind, packs, Iterator) x)
{
    static if (_shouldBeCastedToUnqual!Iterator)
    {
        return .amax(cast(Slice!(N, Unqual!R))x);
    }
    else
    {
        pragma(inline, false);
        alias T = DeepElementType!(typeof(x));
        alias F = realType!T;
        return reduce!(_amax)(F(0), x);
    }
}

///
unittest
{
    import mir.ndslice.allocation: slice;
    auto x = slice!double(6);
    x[] = [0, -1, -2, -7, 6, 2];
    assert(amax(x) == 7);
    // 0 for empty vectors
    assert(amax(x[0 .. 0]) == 0);
}

///
unittest
{
    import mir.ndslice.allocation: slice;

    auto x = slice!cdouble(4);
    x[] = [0 + -1i, -7 + 3i, 2 + 3i, 2 + 2i];

    assert(amax(x) == 10);
    // 0 for empty vectors
    assert(amax(x[$ .. $]) == 0);
}
