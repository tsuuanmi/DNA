//! A small deterministic generator for the renderer: xoshiro256** seeded by
//! `SplitMix64`, normal deviates by Box-Muller, and a hash of a text to a
//! standard normal value for the sequence-determined patterns.

/// A deterministic pseudo-random generator (xoshiro256**).
#[derive(Debug, Clone)]
pub struct Rng {
    state: [u64; 4],
}

/// One `SplitMix64` step: the seeding and hashing mixer.
const fn splitmix(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// A uniform value in (0, 1) from 53 random bits, never 0.
fn open_unit(bits: u64) -> f64 {
    ((bits >> 11) as f64 + 0.5) / (1_u64 << 53) as f64
}

impl Rng {
    /// A generator from `seed`.
    #[must_use]
    pub fn new(seed: u64) -> Self {
        let mut mix = seed;
        Self {
            state: std::array::from_fn(|_| splitmix(&mut mix)),
        }
    }

    /// The next 64 random bits.
    pub const fn next_u64(&mut self) -> u64 {
        let result = self.state[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = self.state[1] << 17;
        self.state[2] ^= self.state[0];
        self.state[3] ^= self.state[1];
        self.state[1] ^= self.state[2];
        self.state[0] ^= self.state[3];
        self.state[2] ^= t;
        self.state[3] = self.state[3].rotate_left(45);
        result
    }

    /// A uniform value in (0, 1).
    pub fn uniform(&mut self) -> f64 {
        open_unit(self.next_u64())
    }

    /// A standard normal value.
    pub fn normal(&mut self) -> f64 {
        let u1 = self.uniform();
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }
}

/// `BLAKE2b`'s initialisation vector (RFC 7693).
const IV: [u64; 8] = [
    0x6A09_E667_F3BC_C908,
    0xBB67_AE85_84CA_A73B,
    0x3C6E_F372_FE94_F82B,
    0xA54F_F53A_5F1D_36F1,
    0x510E_527F_ADE6_82D1,
    0x9B05_688C_2B3E_6C1F,
    0x1F83_D9AB_FB41_BD6B,
    0x5BE0_CD19_137E_2179,
];

/// `BLAKE2b`'s message schedule (RFC 7693).
const SIGMA: [[usize; 16]; 12] = [
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
    [14, 10, 4, 8, 9, 15, 13, 6, 1, 12, 0, 2, 11, 7, 5, 3],
    [11, 8, 12, 0, 5, 2, 15, 13, 10, 14, 3, 6, 7, 1, 9, 4],
    [7, 9, 3, 1, 13, 12, 11, 14, 2, 6, 5, 10, 4, 0, 15, 8],
    [9, 0, 5, 7, 2, 4, 10, 15, 14, 1, 11, 12, 6, 8, 3, 13],
    [2, 12, 6, 10, 0, 11, 8, 3, 4, 13, 7, 5, 15, 14, 1, 9],
    [12, 5, 1, 15, 14, 13, 4, 10, 0, 7, 6, 3, 9, 2, 8, 11],
    [13, 11, 7, 14, 12, 1, 3, 9, 5, 0, 15, 4, 8, 6, 2, 10],
    [6, 15, 14, 9, 11, 3, 0, 8, 12, 2, 13, 7, 1, 4, 10, 5],
    [10, 2, 8, 4, 7, 6, 1, 5, 15, 11, 9, 14, 3, 12, 13, 0],
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
    [14, 10, 4, 8, 9, 15, 13, 6, 1, 12, 0, 2, 11, 7, 5, 3],
];

/// One `BLAKE2b` compression of `block` into `h`, `count` bytes in, `last`
/// for the final block.
fn compress(h: &mut [u64; 8], block: &[u8; 128], count: u128, last: bool) {
    let m: [u64; 16] = std::array::from_fn(|i| {
        let mut word = [0_u8; 8];
        word.copy_from_slice(&block[8 * i..8 * i + 8]);
        u64::from_le_bytes(word)
    });
    let mut v = [0_u64; 16];
    v[..8].copy_from_slice(h);
    v[8..].copy_from_slice(&IV);
    let [low, high] = [count % (1 << 64), count >> 64].map(|part| u64::try_from(part).unwrap_or(0));
    v[12] ^= low;
    v[13] ^= high;
    if last {
        v[14] = !v[14];
    }
    let mut g = |a: usize, b: usize, c: usize, d: usize, x: u64, y: u64| {
        v[a] = v[a].wrapping_add(v[b]).wrapping_add(x);
        v[d] = (v[d] ^ v[a]).rotate_right(32);
        v[c] = v[c].wrapping_add(v[d]);
        v[b] = (v[b] ^ v[c]).rotate_right(24);
        v[a] = v[a].wrapping_add(v[b]).wrapping_add(y);
        v[d] = (v[d] ^ v[a]).rotate_right(16);
        v[c] = v[c].wrapping_add(v[d]);
        v[b] = (v[b] ^ v[c]).rotate_right(63);
    };
    for s in &SIGMA {
        g(0, 4, 8, 12, m[s[0]], m[s[1]]);
        g(1, 5, 9, 13, m[s[2]], m[s[3]]);
        g(2, 6, 10, 14, m[s[4]], m[s[5]]);
        g(3, 7, 11, 15, m[s[6]], m[s[7]]);
        g(0, 5, 10, 15, m[s[8]], m[s[9]]);
        g(1, 6, 11, 12, m[s[10]], m[s[11]]);
        g(2, 7, 8, 13, m[s[12]], m[s[13]]);
        g(3, 4, 9, 14, m[s[14]], m[s[15]]);
    }
    for i in 0..8 {
        h[i] ^= v[i] ^ v[i + 8];
    }
}

/// The unkeyed `BLAKE2b` digest of `data`, 16 bytes (RFC 7693), as Python's
/// `hashlib.blake2b(data, digest_size=16)` gives it.
fn blake2b_16(data: &[u8]) -> [u8; 16] {
    let mut h = IV;
    h[0] ^= 0x0101_0000 ^ 16;
    let mut count = 0_u128;
    let mut chunks = data.chunks(128).peekable();
    if chunks.peek().is_none() {
        compress(&mut h, &[0; 128], 0, true);
    }
    while let Some(chunk) = chunks.next() {
        let mut block = [0_u8; 128];
        block[..chunk.len()].copy_from_slice(chunk);
        count += chunk.len() as u128;
        compress(&mut h, &block, count, chunks.peek().is_none());
    }
    let mut out = [0_u8; 16];
    out[..8].copy_from_slice(&h[0].to_le_bytes());
    out[8..].copy_from_slice(&h[1].to_le_bytes());
    out
}

/// A standard normal value determined by `text` and `seed` alone, so the
/// same context gives the same value in every read: Box-Muller on the
/// `BLAKE2b` digest of `"{seed}:{text}"`, exactly as the research simulator
/// computes it, so both carry one realisation of the sequence-determined
/// patterns.
#[must_use]
pub fn hash_normal(text: &[u8], seed: u64) -> f64 {
    let mut message = format!("{seed}:").into_bytes();
    message.extend_from_slice(text);
    let digest = blake2b_16(&message);
    let word = |bytes: &[u8]| {
        let mut word = [0_u8; 8];
        word.copy_from_slice(bytes);
        u64::from_be_bytes(word)
    };
    let scale = 2.0_f64.powi(64);
    let u1 = (word(&digest[..8]) as f64 + 0.5) / scale;
    let u2 = (word(&digest[8..]) as f64 + 0.5) / scale;
    (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `model._hash_normal("pNACGTAC", 7)` in the research simulator.
    const PY_PNACGTAC: f64 = 1.087_831_887_427_764_6;

    #[test]
    fn deviates_are_deterministic_and_standard() {
        let mut first = Rng::new(7);
        let mut second = Rng::new(7);
        let values = (0..20_000).map(|_| first.normal()).collect::<Vec<_>>();
        assert!(values.iter().all(|&value| value == second.normal()));
        let mean = values.iter().sum::<f64>() / values.len() as f64;
        let variance = values
            .iter()
            .map(|value| (value - mean).powi(2))
            .sum::<f64>()
            / values.len() as f64;
        assert!(mean.abs() < 0.03 && (variance - 1.0).abs() < 0.05);
        assert_ne!(hash_normal(b"ACG", 7), hash_normal(b"ACT", 7));
        // The research simulator's values (model._hash_normal).
        assert_eq!(blake2b_16(b"abc")[..4], [0xCF, 0x4A, 0xB7, 0x91]);
        assert!((hash_normal(b"pNACGTAC", 7) - PY_PNACGTAC).abs() < 1e-12);
    }
}
