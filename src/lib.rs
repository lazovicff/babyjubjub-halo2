// BabyJubJub elliptic curve implementation in Rust.
// For LICENSE check https://github.com/arnaucube/babyjubjub-rs

use halo2curves::group::ff::Field;
use halo2curves::{bn256::Fr, group::ff::PrimeField};

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

#[derive(Clone, Copy, Debug)]
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
        for i in 0..Fr::NUM_BITS {
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

pub fn test_bit(b: &[u8], i: u32) -> bool {
    b[(i / 8) as usize] & (1 << (i % 8)) != 0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add_same_point() {
        let p: PointProjective = PointProjective {
            x: Fr::from_str_vartime(
                "17777552123799933955779906779655732241715742912184938656739573121738514868268",
            )
            .unwrap(),
            y: Fr::from_str_vartime(
                "2626589144620713026669568689430873010625803728049924121243784502389097019475",
            )
            .unwrap(),
            z: Fr::one(),
        };

        let q: PointProjective = PointProjective {
            x: Fr::from_str_vartime(
                "17777552123799933955779906779655732241715742912184938656739573121738514868268",
            )
            .unwrap(),
            y: Fr::from_str_vartime(
                "2626589144620713026669568689430873010625803728049924121243784502389097019475",
            )
            .unwrap(),
            z: Fr::one(),
        };

        let res = p.add(&q).affine();
        assert_eq!(
            res.x,
            Fr::from_str_vartime(
                "6890855772600357754907169075114257697580319025794532037257385534741338397365"
            )
            .unwrap(),
        );
        assert_eq!(
            res.y,
            Fr::from_str_vartime(
                "4338620300185947561074059802482547481416142213883829469920100239455078257889"
            )
            .unwrap(),
        );
    }

    #[test]
    fn test_add_different_points() {
        let p: PointProjective = PointProjective {
            x: Fr::from_str_vartime(
                "17777552123799933955779906779655732241715742912184938656739573121738514868268",
            )
            .unwrap(),
            y: Fr::from_str_vartime(
                "2626589144620713026669568689430873010625803728049924121243784502389097019475",
            )
            .unwrap(),
            z: Fr::one(),
        };

        let q: PointProjective = PointProjective {
            x: Fr::from_str_vartime(
                "16540640123574156134436876038791482806971768689494387082833631921987005038935",
            )
            .unwrap(),
            y: Fr::from_str_vartime(
                "20819045374670962167435360035096875258406992893633759881276124905556507972311",
            )
            .unwrap(),
            z: Fr::one(),
        };

        let res = p.add(&q).affine();
        assert_eq!(
            res.x,
            Fr::from_str_vartime(
                "7916061937171219682591368294088513039687205273691143098332585753343424131937"
            )
            .unwrap(),
        );
        assert_eq!(
            res.y,
            Fr::from_str_vartime(
                "14035240266687799601661095864649209771790948434046947201833777492504781204499"
            )
            .unwrap(),
        );
    }

    #[test]
    fn test_mul_scalar() {
        let p: Point = Point {
            x: Fr::from_str_vartime(
                "17777552123799933955779906779655732241715742912184938656739573121738514868268",
            )
            .unwrap(),
            y: Fr::from_str_vartime(
                "2626589144620713026669568689430873010625803728049924121243784502389097019475",
            )
            .unwrap(),
        };
        let res_m = p.mul_scalar(&Fr::from(3));
        let res_a = p.projective().add(&p.projective());
        let res_a = res_a.add(&p.projective()).affine();
        assert_eq!(res_m.x, res_a.x);

        assert_eq!(
            res_m.x,
            Fr::from_str_vartime(
                "19372461775513343691590086534037741906533799473648040012278229434133483800898"
            )
            .unwrap(),
        );
        assert_eq!(
            res_m.y,
            Fr::from_str_vartime(
                "9458658722007214007257525444427903161243386465067105737478306991484593958249"
            )
            .unwrap(),
        );

        let n = Fr::from_str_vartime(
            "14035240266687799601661095864649209771790948434046947201833777492504781204499",
        )
        .unwrap();
        let res2 = p.mul_scalar(&n);
        assert_eq!(
            res2.x,
            Fr::from_str_vartime(
                "17070357974431721403481313912716834497662307308519659060910483826664480189605"
            )
            .unwrap(),
        );
        assert_eq!(
            res2.y,
            Fr::from_str_vartime(
                "4014745322800118607127020275658861516666525056516280575712425373174125159339"
            )
            .unwrap(),
        );
    }
}
