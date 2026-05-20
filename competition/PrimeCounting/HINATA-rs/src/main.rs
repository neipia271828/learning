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

    let mut prime_list_as_bool = Bitset::new_all_true((size - 1) / 2);

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
    mut sieve_list: Bitset,
) -> usize {
    for &p_usize in &prime_list_as_usize[1..] {
        if p_usize > (high - 1) / p_usize {
            break;
        }
        let p = p_usize;
        let Some(mut x) = find_composite(p, odd_low, high) else {
            continue;
        };

        while x < high {
            let index = (x - odd_low) / 2;
            sieve_list.set_false(index);
            x += 2 * p;
        }
    }

    return sieve_list.count_ones();
}

fn sieving(prime_index: usize, mut v: Bitset) -> Bitset {
    let p = index_to_value(prime_index);
    let start = p * p;
    let start_index = (start - 3) / 2;

    for i in (start_index..v.len()).step_by(p) {
        v.set_false(i);
    }

    return v;
}

fn find_next_prime(last_prime_as_index: usize, v: &Bitset) -> usize {
    for i in last_prime_as_index + 1..v.len() {
        if v.get(i) {
            return i;
        }
    }

    return v.len();
}

fn find_composite(prime_as_num: usize, low: usize, high: usize) -> Option<usize> {
    let p = prime_as_num;

    if p == 0 || high == 0 || p > (high - 1) / p {
        return None;
    }

    let first_multiple = low.div_ceil(p) * p;
    let mut start = first_multiple.max(p * p);

    if start % 2 == 0 {
        start += p;
    }

    if start < high { Some(start) } else { None }
}

fn make_primelist_from_bool(prime_list_as_bool: &Bitset) -> Vec<usize> {
    let mut prime_list_as_usize = vec![2];

    for i in 0..prime_list_as_bool.len() {
        if prime_list_as_bool.get(i) {
            prime_list_as_usize.push(index_to_value(i));
        }
    }

    return prime_list_as_usize;
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

    let sieve_list = Bitset::new_all_true(odd_count);

    range_def_sieving(odd_low, high, prime_list, sieve_list)
}

struct Bitset {
    words: Vec<u64>,
    len: usize,
}

impl Bitset {
    // 全ビットをtrueで初期化
    fn new_all_true(len: usize) -> Self {
        let num_words = len.div_ceil(64);
        let mut words = vec![u64::MAX; num_words];

        // 最後のwordの余りビットを0にする
        let remainder = len % 64;
        if remainder != 0 {
            words[num_words - 1] = (1u64 << remainder) - 1;
        }

        Self { words, len }
    }

    // ビットをfalse(合成数)にする
    #[inline]
    fn set_false(&mut self, index: usize) {
        self.words[index / 64] &= !(1u64 << (index % 64));
    }

    // trueのビット数を数える (popcntが使われる)
    fn count_ones(&self) -> usize {
        self.words.iter().map(|w| w.count_ones() as usize).sum()
    }

    fn get(&self, index: usize) -> bool {
        (self.words[index / 64] >> (index % 64)) & 1 == 1
    }

    fn len(&self) -> usize {
        self.len
    }
}