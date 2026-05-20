use rayon::prelude::*;
use std::time::Instant;

fn main() {
    let start = Instant::now();

    let size: usize = 1_000_000_000_000;
    let size_root: usize = (size as f64).sqrt().ceil() as usize;
    let block_size: usize = (size_root as f64 * 0.9) as usize;

    let prime_list_as_usize = sieve(size_root);

    let count = devision_sieve_method(size, size_root, block_size, prime_list_as_usize);

    let time = start.elapsed();

    println!("{:?}", time);
    println!("{}", count);
}

fn sieve(size: usize) -> Vec<usize> {
    let mut last_prime_as_index: usize = 0;

    let mut prime_list_as_bool: Vec<bool> = make_list(size);

    // O(√n)
    // indexから値への変換公式
    //
    // 現在見ている素数pの二乗以下の間だけ篩を続ける。
    while {
        let p: usize = index_to_value(last_prime_as_index);
        p * p < size
    } {
        prime_list_as_bool = sieving(last_prime_as_index, prime_list_as_bool);

        last_prime_as_index = find_next_prime(last_prime_as_index, &prime_list_as_bool);
    }

    let prime_list_as_usize: Vec<usize> = make_primelist_from_bool(&prime_list_as_bool);

    return prime_list_as_usize;
}

// O(n)
fn make_list(size: usize) -> Vec<bool> {
    return vec![true; (size - 1) / 2];
}

fn devision_sieve_method(
    size: usize,
    size_root: usize,
    block_size: usize,
    prime_list: Vec<usize>,
) -> usize {
    if size < 2 {
        return 0;
    }

    let base_count = prime_list.len();

    // size_root 以下は prime_list.len() で数え済みなので、
    // 次の区間は size_root + 1 から始める
    let start = size_root + 1;
    let end = size + 1;

    // [start, end) を block_size ごとに分ける
    let block_count = end.saturating_sub(start).div_ceil(block_size);

    let segment_count: usize = (0..block_count)
        .into_par_iter()
        .map(|block_index| {
            let low = start + block_index * block_size;
            let high = (low + block_size).min(end);

            count_segment_primes(low, high, &prime_list)
        })
        .sum();

    base_count + segment_count
}

fn range_def_sieving(
    base_1: usize, base_5: usize, high: usize,
    prime_list: &[usize],
    mut sieve_1: Vec<bool>, mut sieve_5: Vec<bool>,
) -> usize {
    // p * p <= high - 1 を満たす素数の上限インデックスを二分探索で1回だけ求める。
    // 旧コードはループ内で毎素数ごとに除算していたが、これで除算が1ブロックに1回になる。
    let prime_end = prime_list.partition_point(|&p| p * p <= high - 1);
    let len_1 = sieve_1.len();
    let len_5 = sieve_5.len();

    for &p in &prime_list[2..prime_end] {

        if let Some(idx) = find_start(p, base_1, 1, high) {
            let mut i = idx;
            while i < len_1 { sieve_1[i] = false; i += p; }
        }
        if let Some(idx) = find_start(p, base_5, 5, high) {
            let mut i = idx;
            while i < len_5 { sieve_5[i] = false; i += p; }
        }
    }

    return counting(&sieve_1) + counting(&sieve_5);
}

fn find_start(p: usize, base: usize, target: usize, high: usize) -> Option<usize> {
    let first_multiple = base.div_ceil(p) * p;
    let mut start = first_multiple.max(p * p);

    let rem = start % 6;
    if rem != target {
        let diff = (target + 6 - rem) % 6;
        // p % 6 が 1 なら inv=1、5 なら inv=5（mod 6 の逆数）
        let inv_p_mod6 = if p % 6 == 1 { 1 } else { 5 };
        let steps = (diff * inv_p_mod6) % 6;
        start += steps * p;
    }

    if start < high { Some((start - base) / 6) } else { None }
}

fn sieving(prime_index: usize, mut v: Vec<bool>) -> Vec<bool> {
    let p = index_to_value(prime_index);
    let start = p * p;
    let start_index = (start - 3) / 2;

    for i in (start_index..v.len()).step_by(p) {
        v[i] = false;
    }

    return v;
}

fn find_next_prime(last_prime_as_index: usize, v: &[bool]) -> usize {
    for i in last_prime_as_index + 1..v.len() {
        if v[i] {
            return i;
        }
    }

    return v.len();
}

// fn find_composite(prime_as_num: usize, low: usize, high: usize) -> Option<usize> {
//     let p = prime_as_num;

//     // p * p <= high - 1 は呼び出し側の partition_point で保証済みなので
//     // 除算によるガードチェックは不要。
//     let first_multiple = low.div_ceil(p) * p;
//     let mut start = first_multiple.max(p * p);

//     if start % 2 == 0 {
//         start += p;
//     }

//     if start < high { Some(start) } else { None }
// }

fn make_primelist_from_bool(prime_list_as_bool: &Vec<bool>) -> Vec<usize> {
    let mut prime_list_as_usize = vec![2];

    for i in 0..prime_list_as_bool.len() {
        if prime_list_as_bool[i] {
            prime_list_as_usize.push(index_to_value(i));
        }
    }

    return prime_list_as_usize;
}

fn counting(v: &[bool]) -> usize {
    return v.iter().filter(|&&is_prime| is_prime).count();
}

fn index_to_value(num: usize) -> usize {
    return 2 * num + 3;
}

fn count_segment_primes(low: usize, high: usize, prime_list: &[usize]) -> usize {
    // [low, high) の中で最初の =1 mod 6と=5 mod 6の値を求める
    let base_1 = first_with_residue(low, 1);
    let base_5 = first_with_residue(low, 5);

    // 各配列の要素数
    let count_1 = if base_1 < high { (high - 1 - base_1) / 6 + 1 } else { 0 };
    let count_5 = if base_5 < high { (high - 1 - base_5) / 6 + 1 } else { 0 };

    let sieve_1 = vec![true; count_1];
    let sieve_5 = vec![true; count_5];

    range_def_sieving(base_1, base_5, high, prime_list, sieve_1, sieve_5)


}

fn first_with_residue(low: usize, target: usize) -> usize {
    let r = low % 6;
    if r == target {
        low
    } else if r < target {
        low + (target - r)
    } else {
        low + (6 - r + target)
    }
}
