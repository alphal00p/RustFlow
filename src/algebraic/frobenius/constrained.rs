//! Exact affine support and numerical matching share the existing Frobenius owner.
use super::*;
use crate::asymptotic::ExactAsymptoticConstraints;
use crate::frobenius::{ExactFrobeniusLimits, exact::Gaussian};

/// Affine maps use column zero for the constant and the remaining columns for
/// the supplied regular boundary. Reconstruction is checked component by component.
pub(crate) struct ConstrainedEndpointMap {
    pub endpoint: EndpointLinearMap,
    pub reconstruction: EndpointLinearMap,
}

impl PreparedAlgebraicFrobenius {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn constrained_endpoint_map(
        &self,
        point: &C,
        matching_parameter: &Atom,
        seeds: &BTreeMap<Symbol, RootSeed>,
        p: Precision,
        order: usize,
        winding: i32,
        constraints: &ExactAsymptoticConstraints,
        limits: &ExactFrobeniusLimits,
        context: &RunContext,
    ) -> Result<Option<ConstrainedEndpointMap>> {
        context.cancellation.check()?;
        let physical = self.lifted.size * self.lifted.count;
        let exact = self.prepared.exact_coefficients(order, limits, context)?;
        constraints.validate(physical, limits, context)?;
        // An incomplete uniformly known prefix asks for more order; actual
        // contradictions and divergent affine directions remain terminal.
        for column in &exact.columns {
            let exponent = crate::frobenius::exact::gaussian(&column.exponent)?;
            if &exponent.re + &Rational::from(order + 1) <= Rational::zero() {
                return Ok(None);
            }
            for relation in &constraints.relations {
                context.cancellation.check()?;
                for (selector, weight) in &relation.terms {
                    if crate::frobenius::exact::gaussian(weight)?.is_zero() {
                        continue;
                    }
                    let difference =
                        crate::frobenius::exact::gaussian(&selector.power)?.re - &exponent.re;
                    if difference.is_integer() && difference > order {
                        return Ok(None);
                    }
                }
            }
        }
        let equations = if self.lifted.monomials.len() > 1 {
            let Some(equations) = self.sheet_equations(
                &exact,
                matching_parameter,
                seeds,
                winding,
                p,
                limits,
                context,
            )?
            else {
                return Ok(None);
            };
            equations
        } else {
            Vec::new()
        };
        let space = exact.constrain_with_rows(constraints, equations, limits, context)?;
        let n = exact.dimension();
        let (finite_offset, finite_directions) =
            exact.finite_projection(&space, physical, context)?;
        let basis = exact.numerical(p);
        let matching = Matrix::from_nested_vec(
            basis.evaluate_with_winding(point, &Default::default(), winding)?,
            FloatField::from_rep(p.zero()),
        )
        .map_err(Error::Numerical)?;
        let (inverse, inverse_error) = conditioned_inverse(&matching, p, context)?;
        let field = FloatField::from_rep(stored_ball(&p.zero(), p));
        let a = matching.map(|value| stored_ball(value, p), field.clone());
        let zero = EpsilonBoundary {
            point: point.clone(),
            leading: 0,
            coefficients: vec![vec![p.zero(); self.lifted.size]; self.lifted.count],
        };
        self.lifted.lift_boundary_data(&zero, seeds, p)?;
        let roots = super::sheet_germ::matching_root_balls(
            &self.lifted,
            matching_parameter,
            seeds,
            p,
            limits,
            context,
        )?;
        let mut lift = Matrix::new(n as u32, physical as u32, field.clone());
        for (block, &mask) in self.lifted.monomials.iter().enumerate() {
            let factor = roots
                .iter()
                .enumerate()
                .filter(|(i, _)| mask & (1 << i) != 0)
                .fold(stored_ball(&p.i(1), p), |a, (_, root)| &a * root);
            for j in 0..physical {
                lift[((block * physical + j) as u32, j as u32)] = factor.clone();
            }
        }
        let lifted_inverse = &inverse * &lift;
        let affine = |offset: &[Gaussian], directions: &[Vec<Gaussian>]| {
            Matrix::from_nested_vec(
                offset
                    .iter()
                    .zip(directions)
                    .map(|(value, row)| {
                        std::iter::once(value)
                            .chain(row)
                            .map(|a| crate::ode::source::exact_ball(a, p))
                            .collect()
                    })
                    .collect(),
                field.clone(),
            )
            .map_err(Error::Numerical)
        };
        let free = space.free_amplitudes();
        let mut amplitudes = Matrix::new((free + 1) as u32, (physical + 1) as u32, field.clone());
        amplitudes[(0, 0)] = stored_ball(&p.i(1), p);
        for (i, &coordinate) in space.free_coordinates.iter().enumerate() {
            for j in 0..physical {
                amplitudes[((i + 1) as u32, (j + 1) as u32)] =
                    lifted_inverse[(coordinate as u32, j as u32)].clone();
            }
        }
        let finite = affine(&finite_offset, &finite_directions)?;
        let reconstructed = &a * &affine(&space.offset, &space.directions)?;
        let make_map = |left: &Matrix<FloatField<ComplexBall>>| -> Result<EndpointLinearMap> {
            context.cancellation.check()?;
            let product = left * &amplitudes;
            let mut values = vec![vec![p.zero(); physical + 1]; physical];
            let mut errors = vec![vec![p.real(0); physical + 1]; physical];
            for i in 0..physical {
                context.cancellation.check()?;
                let row_norm = (1..=free).fold(p.real(0), |sum, j| {
                    sum.add_round(
                        &ball_upper(&left[(i as u32, j as u32)], p),
                        p.bits,
                        RoundingDirection::Up,
                    )
                });
                for j in 0..=physical {
                    let ball = &product[(i as u32, j as u32)];
                    values[i][j] = C::new(ball.re.center.clone(), ball.im.center.clone());
                    errors[i][j] =
                        ball.re
                            .radius
                            .add_round(&ball.im.radius, p.bits, RoundingDirection::Up);
                    if j > 0 {
                        let lift_norm = (0..n).fold(p.real(0), |sum, k| {
                            sum.add_round(
                                &ball_upper(&lift[(k as u32, (j - 1) as u32)], p),
                                p.bits,
                                RoundingDirection::Up,
                            )
                        });
                        errors[i][j] = errors[i][j].add_round(
                            &row_norm
                                .mul_round(&inverse_error, p.bits, RoundingDirection::Up)
                                .mul_round(&lift_norm, p.bits, RoundingDirection::Up),
                            p.bits,
                            RoundingDirection::Up,
                        );
                    }
                }
            }
            Ok(EndpointLinearMap {
                values,
                arithmetic_errors: errors,
            })
        };
        Ok(Some(ConstrainedEndpointMap {
            endpoint: make_map(&finite)?,
            reconstruction: make_map(&reconstructed)?,
        }))
    }
}
