#![no_main]

use libfuzzer_sys::fuzz_target;

// SRS-IN-002 SRS-IN-006
fuzz_target!(|data: &[u8]| {
    dna::fuzzing::parse_abif(data);
});
