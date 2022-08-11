use std::str::FromStr;

use num_bigint::BigUint;
use num_traits::identities::One;
use num_traits::ToPrimitive;
use num_traits::Zero;

// TEST VECTORS
/// 0xe0a77c19a07df2f666ea36f7879462e36fc76959f60cd29ac96341c4ffffffb
/// 6350874878119819312338956282401532410528162663560392320966563075034087161851
///
// const R: Fr = Fr([
//     0xac96341c4ffffffb,
//     12436184717236109307
//
//     0x36fc76959f60cd29,
//     3962172157175319849
//
//     0x666ea36f7879462e,
//     7381016538464732718
//
//     0x0e0a77c19a07df2f,
//     1011752739694698287
// ]);

fn main() {
    let v = BigUint::from_str(
        "16950150798460657717958625567821834550301663161624707787222815936182638968203",
    )
    .unwrap();
    let res = biguint_to_real_u64_vec(v, 4);
    println!("{:?}", res);
}

/// Convert BigUint into a vector of 64-bit limbs.
fn biguint_to_real_u64_vec(mut v: BigUint, limbs: usize) -> Vec<u64> {
    let m = BigUint::one() << 64;
    let mut ret = vec![];

    while v > BigUint::zero() {
        let limb: BigUint = &v % &m;
        ret.push(limb.to_u64().unwrap());
        v >>= 64;
    }

    while ret.len() < limbs {
        ret.push(0);
    }

    assert!(ret.len() == limbs);

    ret
}
