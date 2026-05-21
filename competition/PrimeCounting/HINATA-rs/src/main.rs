use rayon::prelude::*;
use std::time::Instant;

#[global_allocator]
static ALLOC: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

fn main() {
    let start = Instant::now();

    let size: usize = 1_000_000_000;
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
    odd_low: usize,
    high: usize,
    prime_list_as_usize: &[usize],
    mut sieve_list: Vec<bool>,
) -> usize {
    // p * p <= high - 1 を満たす素数の上限インデックスを二分探索で1回だけ求める。
    // 旧コードはループ内で毎素数ごとに除算していたが、これで除算が1ブロックに1回になる。
    let prime_end = prime_list_as_usize.partition_point(|&p| p * p <= high - 1);
    let len = sieve_list.len();

    for &p in &prime_list_as_usize[1..prime_end] {
        let Some(x) = find_composite(p, odd_low, high) else {
            continue;
        };

        // (x - odd_low) / 2 をループ外で1回だけ計算し、
        // 内側ループは index への加算のみにして命令数を減らす。
        let mut index = (x - odd_low) / 2;
        while index < len {
            sieve_list[index] = false;
            index += p;
        }
    }

    return counting(&sieve_list);
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

fn find_composite(prime_as_num: usize, low: usize, high: usize) -> Option<usize> {
    let p = prime_as_num;

    // p * p <= high - 1 は呼び出し側の partition_point で保証済みなので
    // 除算によるガードチェックは不要。
    let first_multiple = low.div_ceil(p) * p;
    let mut start = first_multiple.max(p * p);

    if start % 2 == 0 {
        start += p;
    }

    if start < high { Some(start) } else { None }
}

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
    let odd_low = if low % 2 == 0 { low + 1 } else { low };

    let odd_count = if odd_low >= high {
        0
    } else {
        (high - odd_low + 1) / 2
    };

    let sieve_list = vec![true; odd_count];

    range_def_sieving(odd_low, high, prime_list, sieve_list)
}
