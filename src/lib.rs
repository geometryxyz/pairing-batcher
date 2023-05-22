use std::{collections::HashMap, vec};

use ark_ec::pairing::Pairing;
use ark_ff::One;

pub struct PairingBatcher<E: Pairing> {
    /// Mapping of all G2 points serialized with correlated G1 points
    g2_to_g1: HashMap<E::G2, E::G1>,
    /// challenge
    challenge: E::ScalarField,
    /// running challenge
    running_challenge: E::ScalarField,
}


impl<E: Pairing> PairingBatcher<E> {
    pub fn new(challenge: E::ScalarField) -> Self {
        Self {
            g2_to_g1: HashMap::default(),
            challenge,
            running_challenge: E::ScalarField::one(),
        }
    }

    /// Adds new pairing equation that needs to be checked
    pub fn add_pairing(&mut self, pairs: &[(E::G1Affine, E::G2Affine)]) {
        let g2_points: Vec<E::G2> = pairs.iter().map(|&(_, g2)| g2.into()).collect();

        let mut is_present: bool = false;
        for g2 in g2_points.iter() {
            if self.g2_to_g1.get(g2).is_some() {
                is_present = true;
                break;
            }
        }

        let g1_points: Vec<E::G1> = if is_present {
            self.running_challenge *= self.challenge;
            pairs
                .iter()
                .map(|&(g1, _)| g1 * self.running_challenge)
                .collect()
        } else {
            pairs.iter().map(|pair| pair.0.into()).collect()
        };

        self.update_mapping(&g2_points, &g1_points);
    }

    /// Updates mapping based on pairs that are added
    fn update_mapping(&mut self, g2_points: &[E::G2], g1_points: &[E::G1]) {
        g2_points.iter().zip(
        g1_points
            .iter())
            .for_each(|(&g2, g1)| {
                self.g2_to_g1
                    .entry(g2)
                    .and_modify(|g1_point| *g1_point += g1)
                    .or_insert(*g1);
            });
    }

    /// Returns output that is ready to be called on MultiMillerLoop
    pub fn finalize(&self) -> (Vec<E::G1Prepared>, Vec<E::G2Prepared>) {
        let mut g1_prepared_points = vec![]; 
        let mut g2_prepared_points = vec![];
        self.g2_to_g1
            .iter()
            .for_each(|(g2, g1)| {
                g1_prepared_points.push(g1.into());
                g2_prepared_points.push(g2.into());
            });

        (g1_prepared_points, g2_prepared_points)
    }
}

#[cfg(test)]
mod test {
    use ark_ec::{bls12::{self, G1Prepared, G2Prepared, Bls12, Bls12Config}, Group, pairing::Pairing};
    use ark_std::{test_rng, UniformRand, One};
    use ark_bls12_381::{Fr, G1Affine, G1Projective, G2Affine, G2Projective, Fq12, Bls12_381};
    use ark_ff::Field;
    use ark_ec::pairing::{prepare_g1, prepare_g2};
    use std::ops::Neg;

    use crate::PairingBatcher;

    #[test]
    fn test() {
        /*
            e(a, b) = e(c, d)
            e(j, b) = e(f, g)
            e(e, d) = e(h, b)
        */

        let mut rng = test_rng();

        let a = Fr::rand(&mut rng);
        let b = Fr::rand(&mut rng);
        let c = Fr::rand(&mut rng);
        let d = a * b * c.inverse().unwrap();
        let f = Fr::rand(&mut rng);
        let j = Fr::rand(&mut rng);
        let g = j * b * f.inverse().unwrap();
        let e = Fr::rand(&mut rng);
        let h = e * d * b.inverse().unwrap();

        let a: G1Affine = (G1Projective::generator() * a).into();
        let b: G2Affine = (G2Projective::generator() * b).into();
        let c: G1Affine = (G1Projective::generator() * c).into();
        let d: G2Affine = (G2Projective::generator() * d).into();
        let j: G1Affine = (G1Projective::generator() * j).into();
        let f: G1Affine = (G1Projective::generator() * f).into();
        let g: G2Affine = (G2Projective::generator() * g).into();
        let e: G1Affine = (G1Projective::generator() * e).into();
        let h: G1Affine = (G1Projective::generator() * h).into();

        // Manual Miller loop
        {
            let mlo = {
                Bls12_381::multi_miller_loop(
                    vec![a, c.neg(), j, f.neg(), e, h.neg()],
                    vec![b, d, b, g, d, b]
                )
            };

            let pairing_result =  Bls12_381::final_exponentiation(mlo);
            assert_eq!(pairing_result.unwrap().0, Fq12::one());
        }

        {
            // Batched test
            let mut pairing_batcher = PairingBatcher::<Bls12_381>::new(Fr::rand(&mut rng));

            pairing_batcher.add_pairing(&[(a, b), ((-c), d)]);
            pairing_batcher.add_pairing(&[(j, b), ((-f), g)]);
            pairing_batcher.add_pairing(&[(e, d), ((-h), b)]);

            let batched_tuples = pairing_batcher.finalize();
            /*
                e(a, b) = e(c, d)
                e(j, b) = e(f, g)
                e(e, d) = e(h, b)

                ==>

                e(a + [R]j + [R^2]h, b).e(c + [R^2]e, d).e([R]f, g)
            */
            assert_eq!(3, batched_tuples.0.len());

            let mlo = {
                Bls12_381::multi_miller_loop(
                    batched_tuples.0,
                    batched_tuples.1
                )
            };
            let pairing_result =  Bls12_381::final_exponentiation(mlo);
            assert_eq!(pairing_result.unwrap().0, Fq12::one());
        }
    }
}
