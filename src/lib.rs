// BabyJubJub elliptic curve implementation in Rust.
// For LICENSE check https://github.com/arnaucube/babyjubjub-rs

use halo2curves::bn256::Fr;
use halo2curves::group::ff::Field;
use num_bigint::BigUint;

// D = 168696
pub const D: Fr = Fr::from_raw([0x292F8, 0x00, 0x00, 0x00]);

// A = 168700
pub const A: Fr = Fr::from_raw([0x292FC, 0x00, 0x00, 0x00]);

// 5299619240641551281634865583518297030282874472190772894086521144482721001553
pub const B8_X: Fr = Fr::from_raw([
    0x2893F3F6BB957051,
    0x2AB8D8010534E0B6,
    0x4EACB2E09D6277C1,
    0xBB77A6AD63E739B,
]);

// 16950150798460657717958625567821834550301663161624707787222815936182638968203
pub const B8_Y: Fr = Fr::from_raw([
    0x4B3C257A872D7D8B,
    0xFCE0051FB9E13377,
    0x25572E1CD16BF9ED,
    0x25797203F7A0B249,
]);

pub const B8: Point = Point { x: B8_X, y: B8_Y };

#[derive(Clone, Debug)]
pub struct PointProjective {
    pub x: Fr,
    pub y: Fr,
    pub z: Fr,
}

impl PointProjective {
    pub fn affine(&self) -> Point {
        if bool::from(self.z.is_zero()) {
            return Point {
                x: Fr::zero(),
                y: Fr::zero(),
            };
        }

        let zinv = self.z.invert().unwrap();
        let x = self.x.mul(&zinv);
        let y = self.y.mul(&zinv);

        Point { x, y }
    }

    #[allow(clippy::many_single_char_names)]
    pub fn add(&self, q: &PointProjective) -> PointProjective {
        // add-2008-bbjlp https://hyperelliptic.org/EFD/g1p/auto-twisted-projective.html#doubling-dbl-2008-bbjlp

        // A = Z1*Z2
        let a = self.z.mul(&q.z);
        // B = A2
        let b = a.square();
        // C = X1*X2
        let c = self.x.mul(&q.x);
        // D = Y1*Y2
        let d = self.y.mul(&q.y);
        // E = d*C*D
        let e = D.mul(&c).mul(&d);
        // F = B-E
        let f = b.sub(&e);
        // G = B+E
        let g = b.add(&e);
        // X3 = A*F*((X1+Y1)*(X2+Y2)-C-D)
        let x3 = a
            .mul(&f)
            .mul(&self.x.add(&self.y).mul(&q.x.add(&q.y)).sub(&c).sub(&d));
        // Y3 = A*G*(D-a*C)
        let y3 = a.mul(&g).mul(&d.sub(&A.mul(&c)));
        // Z3 = F*G
        let z3 = f.mul(&g);

        PointProjective {
            x: x3,
            y: y3,
            z: z3,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Point {
    pub x: Fr,
    pub y: Fr,
}

impl Point {
    pub fn projective(&self) -> PointProjective {
        PointProjective {
            x: self.x,
            y: self.y,
            z: Fr::one(),
        }
    }

    pub fn mul_scalar(&self, f: &Fr) -> Point {
        let mut r: PointProjective = PointProjective {
            x: Fr::zero(),
            y: Fr::one(),
            z: Fr::one(),
        };
        let mut exp: PointProjective = self.projective();
        let b = f.to_bytes();
        let n = BigUint::from_bytes_le(&b);
        for i in 0..n.bits() {
            if test_bit(&b, i) {
                r = r.add(&exp);
            }
            exp = exp.add(&exp);
        }
        r.affine()
    }

    pub fn equals(&self, p: Point) -> bool {
        self.x == p.x && self.y == p.y
    }
}

pub fn test_bit(b: &[u8], i: usize) -> bool {
    b[i / 8] & (1 << (i % 8)) != 0
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use super::*;
    use num_bigint::BigUint;
    use num_traits::identities::One;
    use num_traits::ToPrimitive;
    use num_traits::Zero;

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

    #[test]
    fn test_add_same_point() {
        let p_x = BigUint::from_str(
            "17777552123799933955779906779655732241715742912184938656739573121738514868268",
        )
        .unwrap();
        let p_y = BigUint::from_str(
            "2626589144620713026669568689430873010625803728049924121243784502389097019475",
        )
        .unwrap();
        let p_x_limbs = biguint_to_real_u64_vec(p_x, 4);
        let p_y_limbs = biguint_to_real_u64_vec(p_y, 4);

        let p: PointProjective = PointProjective {
            x: Fr::from_raw([p_x_limbs[0], p_x_limbs[1], p_x_limbs[2], p_x_limbs[3]]),
            y: Fr::from_raw([p_y_limbs[0], p_y_limbs[1], p_y_limbs[2], p_y_limbs[3]]),
            z: Fr::one(),
        };

        let q_x = BigUint::from_str(
            "17777552123799933955779906779655732241715742912184938656739573121738514868268",
        )
        .unwrap();
        let q_y = BigUint::from_str(
            "2626589144620713026669568689430873010625803728049924121243784502389097019475",
        )
        .unwrap();
        let q_x_limbs = biguint_to_real_u64_vec(q_x, 4);
        let q_y_limbs = biguint_to_real_u64_vec(q_y, 4);

        let q: PointProjective = PointProjective {
            x: Fr::from_raw([q_x_limbs[0], q_x_limbs[1], q_x_limbs[2], q_x_limbs[3]]),
            y: Fr::from_raw([q_y_limbs[0], q_y_limbs[1], q_y_limbs[2], q_y_limbs[3]]),
            z: Fr::one(),
        };

        let res = p.add(&q).affine();
        let res_x = BigUint::from_str(
            "6890855772600357754907169075114257697580319025794532037257385534741338397365",
        )
        .unwrap();
        let res_y = BigUint::from_str(
            "4338620300185947561074059802482547481416142213883829469920100239455078257889",
        )
        .unwrap();
        let res_x_limbs = biguint_to_real_u64_vec(res_x, 4);
        let res_y_limbs = biguint_to_real_u64_vec(res_y, 4);

        assert_eq!(
            res.x,
            Fr::from_raw([
                res_x_limbs[0],
                res_x_limbs[1],
                res_x_limbs[2],
                res_x_limbs[3]
            ])
        );
        assert_eq!(
            res.y,
            Fr::from_raw([
                res_y_limbs[0],
                res_y_limbs[1],
                res_y_limbs[2],
                res_y_limbs[3]
            ])
        );
    }

    #[test]
    fn test_add_different_points() {
        let p_x = BigUint::from_str(
            "17777552123799933955779906779655732241715742912184938656739573121738514868268",
        )
        .unwrap();
        let p_y = BigUint::from_str(
            "2626589144620713026669568689430873010625803728049924121243784502389097019475",
        )
        .unwrap();
        let p_x_limbs = biguint_to_real_u64_vec(p_x, 4);
        let p_y_limbs = biguint_to_real_u64_vec(p_y, 4);

        let p: PointProjective = PointProjective {
            x: Fr::from_raw([p_x_limbs[0], p_x_limbs[1], p_x_limbs[2], p_x_limbs[3]]),
            y: Fr::from_raw([p_y_limbs[0], p_y_limbs[1], p_y_limbs[2], p_y_limbs[3]]),
            z: Fr::one(),
        };

        let q_x = BigUint::from_str(
            "16540640123574156134436876038791482806971768689494387082833631921987005038935",
        )
        .unwrap();
        let q_y = BigUint::from_str(
            "20819045374670962167435360035096875258406992893633759881276124905556507972311",
        )
        .unwrap();
        let q_x_limbs = biguint_to_real_u64_vec(q_x, 4);
        let q_y_limbs = biguint_to_real_u64_vec(q_y, 4);

        let q: PointProjective = PointProjective {
            x: Fr::from_raw([q_x_limbs[0], q_x_limbs[1], q_x_limbs[2], q_x_limbs[3]]),
            y: Fr::from_raw([q_y_limbs[0], q_y_limbs[1], q_y_limbs[2], q_y_limbs[3]]),
            z: Fr::one(),
        };

        let res = p.add(&q).affine();
        let res_x = BigUint::from_str(
            "7916061937171219682591368294088513039687205273691143098332585753343424131937",
        )
        .unwrap();
        let res_y = BigUint::from_str(
            "14035240266687799601661095864649209771790948434046947201833777492504781204499",
        )
        .unwrap();
        let res_x_limbs = biguint_to_real_u64_vec(res_x, 4);
        let res_y_limbs = biguint_to_real_u64_vec(res_y, 4);
        assert_eq!(
            res.x,
            Fr::from_raw([
                res_x_limbs[0],
                res_x_limbs[1],
                res_x_limbs[2],
                res_x_limbs[3]
            ])
        );
        assert_eq!(
            res.y,
            Fr::from_raw([
                res_y_limbs[0],
                res_y_limbs[1],
                res_y_limbs[2],
                res_y_limbs[3]
            ])
        );
    }

    #[test]
    fn test_mul_scalar() {
        let p_x = BigUint::from_str(
            "17777552123799933955779906779655732241715742912184938656739573121738514868268",
        )
        .unwrap();
        let p_y = BigUint::from_str(
            "2626589144620713026669568689430873010625803728049924121243784502389097019475",
        )
        .unwrap();
        let p_x_limbs = biguint_to_real_u64_vec(p_x, 4);
        let p_y_limbs = biguint_to_real_u64_vec(p_y, 4);
        let p: Point = Point {
            x: Fr::from_raw([p_x_limbs[0], p_x_limbs[1], p_x_limbs[2], p_x_limbs[3]]),
            y: Fr::from_raw([p_y_limbs[0], p_y_limbs[1], p_y_limbs[2], p_y_limbs[3]]),
        };
        let res_m = p.mul_scalar(&Fr::from(3));
        let res_a = p.projective().add(&p.projective());
        let res_a = res_a.add(&p.projective()).affine();
        assert_eq!(res_m.x, res_a.x);

        let res_m_x = BigUint::from_str(
            "19372461775513343691590086534037741906533799473648040012278229434133483800898",
        )
        .unwrap();
        let res_m_y = BigUint::from_str(
            "9458658722007214007257525444427903161243386465067105737478306991484593958249",
        )
        .unwrap();
        let res_m_x_limbs = biguint_to_real_u64_vec(res_m_x, 4);
        let res_m_y_limbs = biguint_to_real_u64_vec(res_m_y, 4);

        assert_eq!(
            res_m.x,
            Fr::from_raw([
                res_m_x_limbs[0],
                res_m_x_limbs[1],
                res_m_x_limbs[2],
                res_m_x_limbs[3]
            ])
        );
        assert_eq!(
            res_m.y,
            Fr::from_raw([
                res_m_y_limbs[0],
                res_m_y_limbs[1],
                res_m_y_limbs[2],
                res_m_y_limbs[3]
            ])
        );

        let n = BigUint::from_str(
            "14035240266687799601661095864649209771790948434046947201833777492504781204499",
        )
        .unwrap();
        let n_limbs = biguint_to_real_u64_vec(n, 4);
        let n_f = Fr::from_raw([n_limbs[0], n_limbs[1], n_limbs[2], n_limbs[3]]);
        let res2 = p.mul_scalar(&n_f);

        let res2_x = BigUint::from_str(
            "17070357974431721403481313912716834497662307308519659060910483826664480189605",
        )
        .unwrap();
        let res2_y = BigUint::from_str(
            "4014745322800118607127020275658861516666525056516280575712425373174125159339",
        )
        .unwrap();
        let res2_x_limbs = biguint_to_real_u64_vec(res2_x, 4);
        let res2_y_limbs = biguint_to_real_u64_vec(res2_y, 4);
        assert_eq!(
            res2.x,
            Fr::from_raw([
                res2_x_limbs[0],
                res2_x_limbs[1],
                res2_x_limbs[2],
                res2_x_limbs[3]
            ])
        );
        assert_eq!(
            res2.y,
            Fr::from_raw([
                res2_y_limbs[0],
                res2_y_limbs[1],
                res2_y_limbs[2],
                res2_y_limbs[3]
            ])
        );
    }
}
