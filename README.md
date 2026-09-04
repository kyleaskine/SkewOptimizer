Skew Optimizer
==============

Command line util to optimize the Murphy E value by adjusting the skew value. <br />

Contents
--------

* msieve_poly.cpp -- code from msieve.
* skewopt.cpp -- main program for skewopt.

Building on Linux
-----------------
Have the gmp library already installed, plus g++.


g++ -c msieve_poly.cpp <br />
g++ skewopt.cpp msieve_poly.o -lgmp -o skewopt <br />

Run the source-built regression test with `./test-skewopt.sh`.

Running
-------
./skewopt <br />
Gives the command line parameters required <br />
<br />

By default, skewopt searches for the skew with the best Murphy E value. To
evaluate a specific skew instead, add it after the final coefficient:<br />
`./skewopt y0 y1 c0 c1 c2 c3 c4 c5 c6 c7 c8 skew`<br />
The optimizing search expands by decades in either direction and refines the
maximum in logarithmic skew space. Automatic optimization is limited to the
reciprocal range `1e-12` through `1e12`; an explicitly requested skew is not
subject to that search limit. If the nominal skew is already outside the
automatic range, it is scored separately so optimization can never return a
worse score.<br />
<br />
For a degree-9 polynomial, use the explicit `-deg9` mode. The existing
degree-8-and-lower command line remains unchanged:<br />
`./skewopt -deg9 y0 y1 c0 c1 c2 c3 c4 c5 c6 c7 c8 c9 [skew]`<br />
<br />

Example
-------
2*x^6 + 1, x - 118571099379011784113736688648896417641748464297615937576404566024103044751294464 <br />
<br />
./skewopt  -118571099379011784113736688648896417641748464297615937576404566024103044751294464 1 1 0 0 0 0 0 2 0 0 <br />
Best Skew: 1.29779341 <br />
