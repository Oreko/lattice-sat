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
               ^ Generate a random instance of CVP
                   ^ Of degree 10

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