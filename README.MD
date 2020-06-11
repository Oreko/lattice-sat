<pre>
This is a class project for CS541 at oregon state university.

The purpose of this tool is to reduce the DCVP (decisional closest vector problem) to the SMT
(satifiable modulo theories) problem and solve it using an SMT solver.


Necessary tools:
    z3
    m4

Build with "cargo build --release"
Executable name: lattice-sat

Test input and general input form can be found in test_data.json
Output format is json as well.

Example usage:
lattice-sat -l -r 10
            ^ Use LLL to reduce the basis
               ^ Generate a random instance of DCVP
                  ^ Of degree 10
lattice-sat -n test_data.json
            ^ Solve over the L2 norm
               ^ On the instance encoded in this file

--help output
USAGE:
    lattice-sat [FLAGS] [OPTIONS] [INPUT]

FLAGS:
    -h, --help       Prints help information
    -l               Use LLL to simplify the basis before running
    -n               Use the l2 norm instead of the linf norm
    -V, --version    Prints version information

OPTIONS:
    -r <DEGREE>        Generate a random degree DEGREE instance of CVP

ARGS:
    <INPUT>    Sets the input file to use



An Aside:

This is my first project using rust, and it was (due to lack of foresight) done in a rushed timeframe. 
As such, be aware that there will be sins and antipatterns, lack of generics, and unsafe code. 
Use and read at your own risk.
</pre>
