use std::env;

fn main() 
{
    let args: Vec<String> = env::args().collect();
    // Vectors in the basis are (101, 0, 0), (0, 208, 0), and (0, 0, 97)
    // Base vector is (543, 600, 220)
    // Distance is 52 (closest integer above) and d^2 = 2696
    // Closest vector should be (505, 624, 194) = (5*x + 3*y + 2*z).

    // We get (x-543)^2 + (y-600)^2 + (z-194)^2 <= 2025

    let base     = vec![543, 600, 220];
    let x0       = vec![101, 0, 0];
    let x1       = vec![0, 208, 0];
    let x2       = vec![0, 0, 97];
    let xs       = vec![&x0, &x1, &x2];
    let distance = 52;
    // Runs in 15s. Small dim, orthogonal and small numbers

    // let base     = vec![543, 600, 220, -420, 62, 234];
    // let x0       = vec![101, 0, 0, 0, 0, 0];
    // let x1       = vec![0, 208, 0, 0, 0, 0];
    // let x2       = vec![0, 0, 97, 0, 0, 0];
    // let x3       = vec![0, 0, 0, -122, 0, 0];
    // let x4       = vec![0, 0, 0, 0, 32, 0];
    // let x5       = vec![0, 0, 0, 0, 0, -42];
    // let xs       = vec![&x0, &x1, &x2, &x3, &x4, &x5];
    // let distance = 110;
    // Runs in 17m43.476s. medium dim, orthogonal and small numbers

    // let base     = vec![543, 600, 220, -420, 62, 234];
    // let x0       = vec![101, -30, 27, 79, 64, 11];
    // let x1       = vec![58, 208, 73, 23, 82, -32];
    // let x2       = vec![35, -56, 97, -83, -78, 57];
    // let x3       = vec![36, 46, 0, -122, -23, -28];
    // let x4       = vec![0, 34, 89, 23, 32, 47];
    // let x5       = vec![-23, -82, 35, 68, 123, -42];
    // let xs       = vec![&x0, &x1, &x2, &x3, &x4, &x5];
    // let distance = 200;
    // Runs in 17m43.476s. medium dim, not orthogonal and small numbers

    // let base     = vec![5432, 6000];
    // let x0       = vec![1001, 0];
    // let x1       = vec![0, 2008];
    // let xs       = vec![&x0, &x1];
    // let distance = 428;
    // 5005, 6024 is a lattice point. So the distance is 428
    // Runs in 39m52.889s. Small dim, orthogonal and larger numbers.

    // let base     = vec![100, 100, 100];
    // let x0       = vec![35, 72, -100];
    // let x1       = vec![-10, 0, -25];
    // let x2       = vec![-20, -279, 678];
    // let xs       = vec![&x0, &x1, &x2];
    // let distance = 1;
    // 100,99,100 is a lattice point. So the distance is 1
    // Runs in 52m55.033s. Small dim, not orthogonal and medium numbers (on average)



    let mut constraint: String = format!("(>= (sq {}) (+ ", distance);

    for idx in 0..x0.len() // fixed dimension
    {
        constraint.push_str("(sq (- (+ ");
        for jdx in 0..xs.len() // fixed coefficient
        {
            let tmp = format!("(* x{} {}) ", jdx, xs[jdx][idx]);
            constraint.push_str(&tmp);
        }
        let tmp = format!(") {})) ", base[idx]);
        constraint.push_str(&tmp);
    }
    constraint.push_str("))");

    println!("testing {}",constraint);
    decision(&constraint, x0.len());
    search(&constraint, x0.len());
}


// Dear rust devs. Please make a do notation. Sincerely, everyone.
// Also, although these functions seem to return a result type, the compiler seems to think they don't. Therefore, the .unwrap()s.
fn decision(constraint: &String, num_vectors: usize)
{
    use rsmt2::Solver;
    let mut solver = Solver::default(()).unwrap();
    solver.set_option(":parallel.enable", true).unwrap();

    solver.define_fun(
        "sq", & [ ("n", "Int") ], "Int", "(* n n)"
    ).unwrap();

    for i in 0..num_vectors
    {
        let tmp = format!("x{}", i);
        solver.declare_const(&tmp, "Int").unwrap();
    }

    solver.assert(constraint).unwrap();

    let is_sat = solver.check_sat().unwrap();
    assert_eq! { is_sat, true }
}

fn search(constraint: &String, num_vectors: usize)
{
    use rsmt2::{ Solver, SmtRes };
    use rsmt2::parse::{ IdentParser, ModelParser };

    #[derive(Clone, Copy)]
    struct Parser;

    impl<'a> IdentParser<String, String, & 'a str> for Parser
    {
        fn parse_ident(self, input: & 'a str) -> SmtRes<String>
        {
            Ok(input.into())
        }
        fn parse_type(self, input: & 'a str) -> SmtRes<String>
        {
            Ok(input.into())
        }
    }

    impl<'a> ModelParser<String, String, String, & 'a str> for Parser 
    {
        fn parse_value(
                self, input: & 'a str,
                _ident: & String, _params: & [ (String, String) ], _typ: & String,
            ) -> SmtRes<String>
        {
            Ok(input.into())
        }
    }

    let mut solver = Solver::default(Parser).unwrap();
    solver.set_option(":parallel.enable", true).unwrap();

    solver.define_fun(
        "sq", & [ ("n", "Int") ], "Int", "(* n n)"
    ).unwrap();

    for i in 0..num_vectors
    {
        let tmp = format!("x{}", i);
        solver.declare_const(&tmp, "Int").unwrap();
    }

    solver.assert(constraint).unwrap();
    

    let is_sat = solver.check_sat().unwrap();
    assert_eq! { is_sat, true }

    let model = solver.get_model().unwrap();
    for (ident, _, kind, value) in model
    {
        println!("{}, {}, {}", ident, kind, value);
    }
}