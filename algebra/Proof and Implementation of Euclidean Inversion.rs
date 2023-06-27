/*
    The Modular Inversion by Means of the Extended Euclidean Algorithm: 
                    Implementation in Rust and Proof

                            Aleksei Vambol
                              June 2023
*/
// Computes the multiplicative inverse of x modulo n by applying the Extended 
// Euclidean Algorithm; panics in the case of n < 2. If x and n are not coprime, 
// the aforementioned inverse does not exist, so None is returned
fn mod_inv(x: i64, n: i64) -> Option<i64> {
    if n < 2 { panic!("The modulus must be greater than 1!"); }
    // Working not with x, but with such x' that 0 <= x' < n and x = x' (mod n)
    let (mut s, mut x_s, mut b, mut x_b) = (((x % n) + n) % n, 1, n, 0);
    // Now s = x', b = n; "s" and "b" stand for "small" and "big", respectively.
    // From now on we have s = x_s * x' + n_s * n and b = x_b * x' + n_b * n,
    // where x' and n are immutable. In each iteration until s = 0 we use   
    // the formula "GCD(s, b) = GCD(b mod s, s)" and update the variables 
    // accordingly. We do not need to store the values of n_s and n_b
    while s > 0 {
        let q = b / s;
        (s, x_s, b, x_b) = (b - q * s, x_b - q * x_s, s, x_s);
    }
    // Now we have b = GCD(0, b) = GCD(x', n). If b > 1, then x' is  
    // not invertible modulo n. If b = 1, then 1 = x_b * x' + n_b * n, 
    // so x_b * x' = 1 (mod n); since it is proven that |x_b| does not 
    // exceed n, we return either x_b or x_b + n
    if b == 1 { Some(if x_b < 0 { x_b + n } else { x_b }) } else { None }
}

fn main() {
    let (x, n) = (3, 10);
    match mod_inv(x, n) {
        Some(r) => println!("{}", r),
        None => println!("{} and {} are not coprime!", x, n),
    }
}