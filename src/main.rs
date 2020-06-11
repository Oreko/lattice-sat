use std::cmp::Ordering;
use std::fs;
use std::time::{SystemTime, UNIX_EPOCH};

use clap::{Arg, App};
use lll_rs::{
    lll::biglll,
    matrix::Matrix,
    vector::BigVector,
};
use rug::{
    Integer,
    rand::RandState
};
use serde::{Deserialize, Serialize};
use serde_json;


#[derive(Deserialize)]
struct CVPInstance {
    basis:    Vec<Vec<i64>>,
    target:   Vec<i64>,
    distance: usize,
}
impl Default for CVPInstance {
    fn default () -> CVPInstance {
        CVPInstance{basis: vec![], target: vec![], distance: 0}
    }
}

#[derive(Serialize)]
struct CVPResult {
    formula:      String,
    decision:     bool,
    basis:        Vec<Vec<i64>>,
    target:       Vec<i64>,
    distance:     usize,
    witness:      Vec<i64>,
    coefficients: Vec<i64>,
}
impl Default for CVPResult {
    fn default () -> CVPResult {
        CVPResult{formula: String::new(), decision: false, basis: vec![], target: vec![], distance: 0, 
                  witness: vec![], coefficients:vec![]}
    }
}


// I'm not sure how to unify the errors from all sources. I could design my own error type, but that seems heavy handed. 
fn main() -> std::io::Result<()>
{
    // CLI arguments
    let matches = App::new("Lattice-Sat")
                          .version("0.0.1")
                          .author("Oreko")
                          .about("Solve small CVP instances using an SMT solver")
                          .arg(Arg::with_name("INPUT")
                               .help("Sets the input file to use")
                               .required(false)
                               .index(1))
                          .arg(Arg::with_name("lll")
                               .short("l")
                               .help("Use LLL to simplify the basis before running"))
                          .arg(Arg::with_name("n2")
                               .short("n")
                               .help("Use the l2 norm instead of the linf norm"))
                          .arg(Arg::with_name("rand")
                               .short("r")
                               .value_name("DEGREE")
                               .help("Generate a random degree DEGREE instance of CVP"))
                          .get_matches();

    let mut data: CVPInstance = CVPInstance::default(); // I'm not sure what the accepted way of doing this is
    if let Some(filename) = matches.value_of("INPUT")
    {
        if matches.is_present("rand")
        {
            println!("Please specify either an input configuration or a random instance.");
            return Ok(());
        }
        let json_string: String = fs::read_to_string(filename)?;
        data = serde_json::from_str(&json_string)?;
    } else if matches.is_present("rand")
    {
        let degree_str = matches.value_of("rand").unwrap(); // Unwrap is safe since we checked for rand's presence
        let degree = degree_str.parse::<usize>().unwrap_or(1);
        if degree > 1 // This felt cleaner than the corresponding match and check.
        {
            let l2_norm = matches.is_present("n2");
            data = generate_random_instance(degree, l2_norm);
        } else
        {
            println!("Please specify a positive integer greater than 1 for DEGREE.");
            return Ok(());
        }
        
    } else
    {
        println!("Please specify either an input configuration or a random instance.");
        return Ok(());
    }

    let m: usize = data.basis.len();
    let n: usize = data.basis[0].len();
    let mut basis: Matrix<BigVector> = Matrix::init(m, n);

    assert_eq!(n, data.target.len());

    for i in 0..m // For some reason, map() based conversion gave errors
    {
        let mut vec: Vec<Integer> = Vec::new();
        for j in 0..n
        {
            vec.push(Integer::from(data.basis[i][j])); 
        }
        basis[i] = BigVector::from_vector(vec);
    }

    if matches.is_present("lll")
    {
        // This sometimes fails with a divide by zero error?
        // I wonder if this is a problem with the library or a problem with my code.
        // If I come back to this I'll look into it more.
        biglll::lattice_reduce(&mut basis);
    }

    let mut constraint: String = String::new();
    if matches.is_present("n2")
    {
        constraint.push_str(&l2_constraint_builder(&basis, &data.target, data.distance));
    } else
    {
        constraint.push_str(&linf_constraint_builder(&basis, &data.target, data.distance));
    }
    match search(&constraint, basis.dimensions().0)
    {
        Ok(mut result) =>
        {
            result.witness = vec_matrix_mult(&result.coefficients, &data.basis).expect("Caught integer over/underflow");
            result.target = data.target;
            result.distance = data.distance;
            result.basis = data.basis;
            result.formula.push_str(&constraint);
            let j = serde_json::to_string_pretty(&result).unwrap(); // Safe unwrap
            println!("{}",j);
        }
        Err(e) =>
        {
            // I'd like to propagate this error, but the types don't unify.
            println!("{}", e);
        }
    }
    return Ok(());
}

fn triangle(n: usize) -> usize
{
    let mut res = 0;
    for i in 1..n+1
    {
        res += i;
    }
    return res;
}

fn generate_random_instance(degree: usize, l2_norm: bool) -> CVPInstance
{
    let mut instance = CVPInstance::default();

    // Set up our random number generator
    let seed = Integer::from(SystemTime::now().duration_since(UNIX_EPOCH).expect("clock error").as_millis());
    let mut rand = RandState::new_mersenne_twister();
    rand.seed(&seed);

    // Number of non 0/1 values in the two matrices for degree "degree"
    let per_matrix = triangle(degree-1);
    let rands_needed = 2*per_matrix;
    let mut unimodular: Vec<i64> = vec![0; rands_needed];
    let mut base: Vec<i64>       = vec![];
    for i in 0..(rands_needed)
    {
        // This will generate numbers between -5 and 5
        // This follows the ideas proposed in the GGH paper
        let mut elem = rand.below(5) as i64;
        let sign = rand.bits(1) as i64;
        if sign == 1
        {
            elem = elem * -1;
        }
        unimodular[i] = elem;
    }
    for _ in 0..degree
    {
        // Random numbers for the "good" basis
        let mut basis_elem = rand.below(50) as i64;
        let basis_sign = rand.bits(1);
        if basis_sign == 1
        {
            basis_elem = basis_elem * -1;
        }
        base.push(basis_elem);

        // Random numbers for the target vector
        let mut target_elem = rand.below(50) as i64;
        let target_sign = rand.bits(1);
        if target_sign == 1
        {
            target_elem = target_elem * -1;
        }
        instance.target.push(target_elem);
    }
    // Generate the basis

    // Generate the two unimodular matrices
    let mut u1: Vec<Vec<i64>> = vec![];
    for i in 0..degree
    {
        let mut inner: Vec<i64> = vec![0; degree];
        let sign = rand.bits(1);
        if sign == 1
        {
            inner[i] = -1;
        } else
        {
            inner[i] = 1;
        }
        
        for j in 0..i
        {
            // Loop skips for i = 0, so triangle(i-1) is defined.
            inner[j] = unimodular[triangle(i-1) + j];
        }
        u1.push(inner);
    }

    let mut u2: Vec<Vec<i64>> = vec![];
    let mut offset: usize = 0;
    for i in 0..degree
    {
        let mut inner: Vec<i64> = vec![0; degree];
        let sign = rand.bits(1);
        if sign == 1
        {
            inner[i] = -1;
        } else
        {
            inner[i] = 1;
        }
        
        
        for j in 1..degree-i
        {
            inner[i+j] = unimodular[per_matrix - 1 + offset + j];
        }
        offset += degree-1-i;
        u2.push(inner);
    }
    
    // Generate the original basis
    let mut b: Vec<Vec<i64>> = vec![];
    for i in 0..degree
    {
        let mut inner: Vec<i64> = vec![0; degree];
        inner[i] = base[i];
        b.push(inner);
    }

    // Calculate the norm using the "good" basis
    if l2_norm
    {
        // Calculate the l2 norm
        let mut sum_of_squares = 0;
        for i in 0..degree
        {
            // This isn't the tightest bound, but given the basis obfuscation it should be good enough.
            sum_of_squares += (b[i][i] - instance.target[i]).pow(2);
        }
        instance.distance = integer_sqrt(sum_of_squares).unwrap() as usize; // this coversion should be ok. We know the size and sign of the data.
    } else
    {
        // Calculate the linf norm
        let mut norms: Vec<usize> = vec![0; degree];
        for i in 0..degree
        {
            // This isn't the tightest bound, but given the basis obfuscation it should be good enough.
            let tmp = b[i][i] - instance.target[i];
            norms[i] = tmp.abs() as usize;
        }
        instance.distance = *norms.iter().max().unwrap(); // Both the unwrap and dereference are safe here.
    }
    
    // Basis = U1 U2 B
    let u = matrix_mult(&u1, &u2).expect("Caught integer over/underflow");
    instance.basis = matrix_mult(&u, &b).expect("Caught integer over/underflow");
    
    return instance;
}

// This is a quick and dirty algorithm and can be optimized using euclidian division for instance.
fn integer_sqrt(a: i64) -> Option<i64>
{
    if a < 0 {return None}

    let mut i: i64 = 1;
    let mut result: i64 = 1; 
    while result <= a
    { 
      i += 1; 
      result = i * i; 
    } 

    return Some(i+1);
}

// We trust the input data. I may fix this if I come back to this code.
fn matrix_mult(a: &Vec<Vec<i64>>, b: &Vec<Vec<i64>>) -> Option<Vec<Vec<i64>>>
{
    let mut result: Vec<Vec<i64>> = vec![];
    let n: usize = a.len();
    let m: usize = a[0].len();
    let l: usize = b[0].len();
    for i in 0..n
    {
        result.push(vec![0; l]);
        for j in 0..l
        {
            for k in 0..m
            {
                let product = a[i][k].checked_mul(b[k][j])?;
                let sum = result[i][j].checked_add(product)?;
                result[i][j] = sum;
            }
            
        }
    }
    return Some(result);
}

// We trust the input data. I may fix this if I come back to this code.
// 1xM times MxN
fn vec_matrix_mult(a: &Vec<i64>, b: &Vec<Vec<i64>>) -> Option<Vec<i64>>
{
    let m: usize = a.len();
    let n: usize = b[0].len();
    let mut result: Vec<i64> = vec![0; m];
    for i in 0..n
    {
        for j in 0..m
        {
            let product = a[j].checked_mul(b[j][i])?;
            let sum = result[i].checked_add(product)?;
            result[i] = sum;
        }
    }
    return Some(result);
}

fn l2_constraint_builder(basis: &Matrix<BigVector>, target: &Vec<i64>, distance: usize) -> String
{
    let dims = basis.dimensions();
    let mut l2_constraint: String = format!("(>= (sq {}) (+ ", distance);

    for idx in 0..dims.1 // fixed dimension
    {
        l2_constraint.push_str("(sq (- (+ ");
        for jdx in 0..dims.0 // fixed coefficient
        {
            let tmp = format!("(* x{} {}) ", jdx, basis[jdx][idx]);
            l2_constraint.push_str(&tmp);
        }
        let tmp = format!(") {})) ", target[idx]);
        l2_constraint.push_str(&tmp);
    }
    l2_constraint.push_str("))");

    return l2_constraint;
}

fn linf_constraint_builder(basis: &Matrix<BigVector>, target: &Vec<i64>, distance: usize) -> String
{
    let dims = basis.dimensions();
    let mut linf_constraint: String = "(and ".to_string();

    for idx in 0..dims.1 // fixed dimension
    {
        let tmp = format!("(>= {} (abs (- {} (+ ", distance, target[idx]);
        linf_constraint.push_str(&tmp);
        for jdx in 0..dims.0 // fixed coefficient
        {
            let tmp = format!("(* x{} {}) ", jdx, basis[jdx][idx]);
            linf_constraint.push_str(&tmp);
        }
        linf_constraint.push_str("))))");
    }
    linf_constraint.push_str(")");

    return linf_constraint;
}


// Dear rust devs. Please make a do notation. Sincerely, everyone.
fn search(constraint: &String, num_vectors: usize) -> Result<CVPResult, rsmt2::errors::Error>
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

    let mut solver = Solver::default(Parser)?;
    solver.set_option(":parallel.enable", true)?;

    solver.define_fun(
        "sq", & [ ("n", "Int") ], "Int", "(* n n)"
    )?;

    solver.define_fun(
        "abs", & [ ("n", "Int") ], "Int", "(ite (< n 0) (- 0 n) n)"
    )?;

    for i in 0..num_vectors
    {
        let tmp = format!("x{}", i);
        solver.declare_const(&tmp, "Int")?;
    }

    solver.assert(constraint)?;
    
    let mut result: CVPResult = CVPResult::default();

    match solver.check_sat()
    {
        Ok(res) => result.decision = res,
        Err(err) => return Err(err),
    }
    if result.decision
    {
        let mut model = solver.get_model()?;
        // I don't trust rsmt2 to provide the elements of the model in any particular order.
        // We know that each identifier will be "x#" and unique, so lexicographic sort will work.
        model.sort_by(|a, b| cmp_coef(&a.0, &b.0));
        for (_, _, _, value) in model
        {
            match value.parse::<i64>()// I wonder if there's a cleaner way of doing this without writing my own parser.
            {
                Ok(number) => result.coefficients.push(number),
                Err(_) => // Happens in the case of negative values. Then we need to parse it ourselves.
                {
                    // Negative strings come in the form (- #####)
                    // So we'll simply strip whitespace and paranthesis.
                    let value = value.replace(&['(', ')', ' '][..], "");
                    let res = value.parse::<i64>().expect("Caught integer over/underflow, coefficient was too large.");
                    result.coefficients.push(res);
                }
            }
        }
    }
    return Ok(result);
}

fn cmp_coef(a: &str, b:&str) -> Ordering
{
    let first  = &a[1..a.len()].parse::<i64>().unwrap(); // Potentially unsafe due to overflows.
    let second = &b[1..b.len()].parse::<i64>().unwrap(); // Potentially unsafe due to overflows.

    return first.cmp(second);
}