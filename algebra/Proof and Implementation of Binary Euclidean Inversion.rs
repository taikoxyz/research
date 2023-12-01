/*
The Modular Inversion by Means of the Binary Extended Euclidean Algorithm: 
                    Implementation in Rust and Proof

                             Aleksei Vambol
                               June 2023
*/
// Computes the multiplicative inverse of x modulo n by applying the binary  
// Extended Euclidean Algorithm. For applying this method n must be odd,   
// x and n must be coprime (because if they are not coprime, the inverse  
// does not exist), both x and n must be positive
fn mod_inv(x: i64, n: i64) -> i64 {
    let (mut a, mut b, mut u, mut v) = (x, n, 1, 0);
    // Now a = x, b = n;
    // (1) b is odd; (2) u < n and v < n; (3) a = u * x (mod n);  
    // (4) b = v * x (mod n); (5) GCD(a, b) = GCD(x, n) = 1; 
    // (6) a, b, u and v are non-negative. 
    // In each iteration we perform the transformation of a and b as well as  
    // their accompanying coefficients u and v, which preserves (1)-(6) and 
    // decreases a + b. Thus, when a is 0, b is 1 due to (5), so v is the 
    // inverse of x modulo n due to (4). Also, 0 < v < n due to (2) and (6). 
    // In the case of the classic Extended Euclidean Algorithm we would have 
    // a = u * x + i * n and b = v * x + j * n instead of (3) and (4), only (5) 
    // would still hold true. Since we do not seek all the Bezout coefficients 
    // and have "(mod n)" in both (3) and (4), we discard "i * n" and "j * n" 
    while a > 0 {
        if (a & 1) > 0 {
            // Both a and b are odd here. We decrease the greatest by the   
            // smallest and satisfy (5), because GCD(p, q) = GCD(p - q, q),  
            // update the greatest's accompanying coefficient to satisfy   
            // (3) and (4), swap the values for a and b as well as for their  
            // accompanying coefficients, if this is required, to satisfy (1)  
            // without breaking (3)-(5). Thus, (1)-(5) are satisfied
            if a >= b {
                (a, u) = (a - b, u - v);
            } else {
                (a, b, u, v) = (b - a, a, v - u, u);
            }
            // We conditionally update u to satisfy (6) without breaking (1)-(5)
            if u < 0 { u += n; }
        }
        // Here a is even and (1)-(6) are satisfied. We divide a by 2 and still 
        // satisfy (5), since b is odd due to (1) and GCD(p, q) = GCD(p / 2, q) 
        // for even p and odd q. As a result, only (3) is not satisfied 
        a >>= 1;
        // In order to satisfy (3) without breaking (1)-(2) and (4)-(6),   
        // u should be set to u * 2^(-1) mod n. If u is even, it is done by  
        // dividing u by 2. For odd u we set u to (u + n) / 2, since n is odd, 
        // u < n due to (2) and u is non-negative due to (6)
        if u & 1 > 0 { u += n; }
        u = u >> 1;
    }
    v
}

fn main() {
    let (x, n) = (13, 97);
    let i = mod_inv(x, n);
    assert!((i * x) % n == 1, "Incorrect inverse!");
    println!("{}", i);
}
