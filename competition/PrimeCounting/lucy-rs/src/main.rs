use std::time::Instant;

fn main() {
    let start = Instant::now();
    let n: u64 = 1_000_000_000_000;
    let count = prime_count(n);
    let elapsed = start.elapsed();
    println!("{:?}", elapsed);
    println!("{}", count);
}

// Lucy_Hedgehog法による素数計数
// floor(N/k) の形をとる全ての値に対して S(v) = π(v) を同時に求める
//
// 初期化: S(v) = v - 1
// 各素数 p について:
//   S(v) -= S(v/p) - S(p-1)  ← p が最小素因数の合成数を除去
fn prime_count(n: u64) -> i64 {
    if n < 2 {
        return 0;
    }

    let sqn = (n as f64).sqrt() as usize;

    // small[i] = S(i),          i <= sqn
    // large[k] = S(floor(n/k)), k <= sqn
    let mut small = vec![0i64; sqn + 1];
    let mut large = vec![0i64; sqn + 1];

    for i in 1..=sqn {
        small[i] = i as i64 - 1;
    }
    for k in 1..=sqn {
        large[k] = (n / k as u64) as i64 - 1;
    }

    let mut is_prime = vec![true; sqn + 1];

    for p in 2..=sqn {
        if !is_prime[p] {
            continue;
        }

        let p_u64 = p as u64;
        let p2 = p_u64 * p_u64;
        let pcnt = small[p - 1];

        // large[k] の更新: floor(n/k) >= p^2 の範囲のみ
        for k in 1..=sqn {
            let v = n / k as u64;
            if v < p2 {
                break;
            }
            let kp = k * p;
            let svp = if kp <= sqn {
                large[kp]
            } else {
                small[(v / p_u64) as usize]
            };
            large[k] -= svp - pcnt;
        }

        // small[i] の更新: i in [p^2, sqn]、大きい方から処理
        if p2 <= sqn as u64 {
            for i in (p2 as usize..=sqn).rev() {
                small[i] -= small[i / p] - pcnt;
            }
            for j in (p2 as usize..=sqn).step_by(p) {
                is_prime[j] = false;
            }
        }
    }

    large[1]
}
