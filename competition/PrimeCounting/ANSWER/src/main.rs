use std::collections::HashMap;
use std::time::Instant;

fn main() {
    let start = Instant::now();
    let n: u64 = 1_000_000_000_000;
    let count = prime_pi(n);
    println!("{:?}", start.elapsed());
    println!("{}", count);
}

// Lehmer の公式:
//   π(n) = φ(n, a) + a - 1 - P₂ - P₃
//   a = π(n^{1/4})
//
// P₂ = Σ_{i=a+1}^{b}  ( π(n/pᵢ) - (i-1) )          b = π(√n)
// P₃ = Σ_{i=a+1}^{c}  Σ_{j=i}^{π(√(n/pᵢ))}  ( π(n/(pᵢpⱼ)) - (j-1) )   c = π(n^{1/3})
//
// π(v) は Lucy DP で O(√n) 個の floor(n/k) 全てについて O(n^{3/4}/log n) で前計算。
// φ(x, a) は再帰 + メモ化。pₐ² > x のとき π(x) - a + 1 で打ち切り。

fn prime_pi(n: u64) -> i64 {
    if n < 2 { return 0; }

    let sqn = (n as f64).sqrt() as usize;

    // ── Lucy DP ──────────────────────────────────────────────────────────
    let mut small = vec![0i64; sqn + 2]; // small[v] = π(v)
    let mut large = vec![0i64; sqn + 2]; // large[k] = π(⌊n/k⌋)

    for i in 1..=sqn { small[i] = i as i64 - 1; }
    for k in 1..=sqn { large[k] = (n / k as u64) as i64 - 1; }

    let mut is_prime = vec![true; sqn + 2];
    let mut primes: Vec<usize> = vec![];

    for p in 2..=sqn {
        if !is_prime[p] { continue; }
        primes.push(p);
        let (pu, p2, pcnt) = (p as u64, (p * p) as u64, small[p - 1]);

        for k in 1..=sqn {
            let v = n / k as u64;
            if v < p2 { break; }
            let kp = k * p;
            let svp = if kp <= sqn { large[kp] } else { small[(v / pu) as usize] };
            large[k] -= svp - pcnt;
        }
        if p2 <= sqn as u64 {
            for i in (p2 as usize..=sqn).rev() { small[i] -= small[i / p] - pcnt; }
            for j in (p2 as usize..=sqn).step_by(p) { is_prime[j] = false; }
        }
    }

    // ── π ルックアップ ───────────────────────────────────────────────────
    let pi = |v: u64| -> i64 {
        if v <= sqn as u64 { small[v as usize] } else { large[(n / v) as usize] }
    };

    // ── 閾値 ─────────────────────────────────────────────────────────────
    let iroot = |k: u32| -> u64 {
        let mut r = (n as f64).powf(1.0 / k as f64) as u64;
        while r.saturating_pow(k) > n { r -= 1; }
        while (r + 1).saturating_pow(k) <= n { r += 1; }
        r
    };
    let a = pi(iroot(4)) as usize; // π(n^{1/4})
    let b = pi(sqn as u64) as usize; // π(√n)
    let c = pi(iroot(3)) as usize; // π(n^{1/3})

    // ── P₂ ───────────────────────────────────────────────────────────────
    let p2: i64 = (a + 1..=b)
        .map(|i| pi(n / primes[i - 1] as u64) - (i as i64 - 1))
        .sum();

    // ── P₃ ───────────────────────────────────────────────────────────────
    let mut p3: i64 = 0;
    for i in (a + 1)..=c {
        let xi = n / primes[i - 1] as u64;
        let bi = pi(isqrt(xi)) as usize;
        for j in i..=bi {
            p3 += pi(xi / primes[j - 1] as u64) - (j as i64 - 1);
        }
    }

    // ── φ(n, a) ──────────────────────────────────────────────────────────
    let mut cache = HashMap::new();
    let phi_na = phi(n, a, &primes, &small, &large, n, sqn, &mut cache);

    phi_na + a as i64 - 1 - p2 - p3
}

fn isqrt(n: u64) -> u64 {
    let mut r = (n as f64).sqrt() as u64;
    while r * r > n { r -= 1; }
    while (r + 1) * (r + 1) <= n { r += 1; }
    r
}

fn phi(
    x: u64, a: usize,
    primes: &[usize],
    small: &[i64], large: &[i64],
    n: u64, sqn: usize,
    cache: &mut HashMap<(u64, usize), i64>,
) -> i64 {
    if a == 0 { return x as i64; }
    if x == 0 { return 0; }

    let pa = primes[a - 1] as u64;
    if pa * pa > x {
        // この区間に残る数は 1 とpₐより大きい素数のみ → π(x) - a + 1
        let pi_x = if x <= sqn as u64 { small[x as usize] } else { large[(n / x) as usize] };
        return pi_x - a as i64 + 1;
    }

    if let Some(&v) = cache.get(&(x, a)) { return v; }
    let result = phi(x, a - 1, primes, small, large, n, sqn, cache)
               - phi(x / pa, a - 1, primes, small, large, n, sqn, cache);
    cache.insert((x, a), result);
    result
}
