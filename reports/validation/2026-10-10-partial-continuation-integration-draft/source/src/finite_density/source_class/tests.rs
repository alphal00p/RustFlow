use super::super::guarded::{GuardedIdentity, GuardedIdentityTerm, GuardedMeasureIdentity};
use super::*;
use rustred::solver::guarded::GuardedSource;
fn data() -> (
    [Symbol; 6],
    [IndexRole; 6],
    IndexDomain<6>,
    SourceImageBinding,
) {
    let indices = std::array::from_fn(|i| symbol!(format!("source_class_test::n{i}")));
    let roles = [
        IndexRole::Ordinary,
        IndexRole::RequiredCut,
        IndexRole::Ordinary,
        IndexRole::Occupation,
        IndexRole::Occupation,
        IndexRole::Ordinary,
    ];
    let a = IndexDomain::new([
        IndexBounds::unbounded(),
        IndexBounds::unbounded(),
        IndexBounds::new(None, Some(0)).unwrap(),
        IndexBounds::new(Some(0), None).unwrap(),
        IndexBounds::new(Some(0), None).unwrap(),
        IndexBounds::fixed(0),
    ])
    .unwrap();
    let b = SourceImageBinding {
        input_identity: "synthetic-class-only".into(),
        family_signature: vec!["synthetic-role-layout".into()],
        origin_identity: "synthetic-domain-check-not-physical-permit".into(),
        shifted_slots: vec![0],
        source_options: Default::default(),
        expected_roles: roles[..5].to_vec(),
        physical_arity: 5,
        physical_slots: 2,
        input_slots: 3,
    };
    (indices, roles, a, b)
}
fn context(
    coefficient: Atom,
    shift: [i16; 6],
    domain: IndexDomain<6>,
    conditions: Vec<Atom>,
) -> GuardedContext<6> {
    let (indices, roles, _, _) = data();
    GuardedContext::new_with_physical_arity(
        GuardedMeasureIdentity {
            measure: "source-class-test".into(),
            support: "only synthetic class regression".into(),
            orientation: "positive".into(),
            normalization: "unit".into(),
            branch: "formal".into(),
            deformation: "slot0".into(),
        },
        roles,
        indices,
        vec![symbol!("source_class_test::x")],
        vec![GuardedIdentity {
            id: "row".into(),
            domain,
            terms: vec![GuardedIdentityTerm { shift, coefficient }],
            nonzero_conditions: conditions,
        }],
        5,
    )
    .unwrap()
}
fn cert(c: &GuardedContext<6>) -> Result<SourceImageClassCertificate> {
    let (_, _, a, b) = data();
    certify_source_image_class(c, &a, &b, Default::default(), &RunContext::default())
}
fn copy_sources(
    c: &GuardedContext<6>,
    change: impl FnOnce(&mut Vec<GuardedSource<6>>),
) -> Arc<GuardedSourceSystem<6>> {
    let original = c.sources();
    let mut rows = original
        .sources()
        .iter()
        .zip(original.native_sources().rows())
        .map(|(i, r)| {
            GuardedSource::new(i.id.clone(), r.clone(), i.domain.clone())
                .with_nonzero_conditions(i.nonzero_conditions.clone())
        })
        .collect::<Vec<_>>();
    change(&mut rows);
    Arc::new(
        GuardedSourceSystem::new(
            original.measure_id(),
            *original.roles(),
            *original.native_sources().index_variables(),
            rows,
        )
        .unwrap(),
    )
}
#[test]
fn source_class_completion_and_whole_shift_boundary_zero_are_exact() {
    let (indices, _, a, _) = data();
    let n = Atom::var(indices[2]);
    let c = context(-&n, [0, 0, 1, 0, 0, 0], a.clone(), vec![]);
    assert_eq!(cert(&c).unwrap().counts.checked_faces, 1);
    let mut b = *a.bounds();
    b[2] = IndexBounds::new(None, Some(-1)).unwrap();
    let c = context(
        -(&n + Atom::one()),
        [0, 0, 2, 0, 0, 0],
        IndexDomain::new(b).unwrap(),
        vec![],
    );
    assert_eq!(cert(&c).unwrap().counts.checked_faces, 1);
    let c = context(
        &n * (&n + Atom::one()),
        [0, 0, 2, 0, 0, 0],
        a.clone(),
        vec![],
    );
    assert_eq!(cert(&c).unwrap().counts.checked_faces, 2);
    for shift in [1, 2] {
        assert!(
            cert(&context(
                Atom::one(),
                [0, 0, shift, 0, 0, 0],
                a.clone(),
                vec![]
            ))
            .is_err()
        );
    }
}
#[test]
fn source_class_all_fixed_coordinates_are_specialized_before_image_admission() {
    let (indices, _, a, _) = data();
    let mut g = *a.bounds();
    g[0] = IndexBounds::fixed(7);
    let c = context(
        Atom::var(indices[0]) - Atom::num(7),
        [0, 0, 1, 0, 0, 0],
        IndexDomain::new(g).unwrap(),
        vec![],
    );
    let proof = cert(&c).unwrap();
    assert_eq!(proof.counts.merged_terms, 0);
}
#[test]
fn source_class_invalid_tuple_is_not_hidden_by_required_cut_zero() {
    let (_, _, a, b) = data();
    let mut g = *a.bounds();
    g[1] = IndexBounds::fixed(0);
    g[3] = IndexBounds::fixed(0);
    let g = IndexDomain::new(g).unwrap();
    assert!(cert(&context(Atom::one(), [0, 0, 1, 0, 0, 0], g.clone(), vec![])).is_err());
    assert!(
        cert(&context(
            Atom::one(),
            [0, 0, 0, -1, 0, 0],
            g.clone(),
            vec![]
        ))
        .is_err()
    );
    let c = context(Atom::one(), [0, -1, 0, 0, 0, 0], g, vec![]);
    cert(&c).unwrap();
    let escaped = copy_sources(&c, |rows| {
        rows[0].row[0].integral = NativeIntegral::symbolic([0, -1, 0, 0, 0, 1]).unwrap()
    });
    assert!(
        certify_native(
            &escaped,
            5,
            &a,
            &b,
            Default::default(),
            &RunContext::default()
        )
        .is_err()
    );
}
#[test]
fn source_class_conditions_and_complete_context_are_bound_without_erasure() {
    let (indices, _, a, b) = data();
    let n = Atom::var(indices[2]);
    let x = Atom::var(symbol!("source_class_test::x"));
    let c = context(-&n / &x, [0, 0, 1, 0, 0, 0], a.clone(), vec![x.clone()]);
    let before = c.sources().sources()[0].nonzero_conditions.clone();
    let proof = cert(&c).unwrap();
    assert!(proof.counts.source_conditions > 0);
    assert_eq!(before, c.sources().sources()[0].nonzero_conditions);
    proof
        .validate_context(&c, &a, &b, &RunContext::default())
        .unwrap();
    let changed = context(
        -&n / &x,
        [0, 0, 1, 0, 0, 0],
        a.clone(),
        vec![&x + Atom::one()],
    );
    assert!(
        proof
            .validate_context(&changed, &a, &b, &RunContext::default())
            .is_err()
    );
    for mode in 0..4 {
        let mut other = b.clone();
        match mode {
            0 => other.origin_identity.push('x'),
            1 => other.family_signature.push("other".into()),
            2 => other.input_identity.push('x'),
            _ => other.source_options.free_virtual_zero_sectors = true,
        };
        assert!(
            proof
                .validate_context(&c, &a, &other, &RunContext::default())
                .is_err()
        );
    }
}
#[test]
fn source_class_limits_and_cancel_fail_before_a_certificate() {
    let (indices, _, a, b) = data();
    let c = context(
        -Atom::var(indices[2]),
        [0, 0, 1, 0, 0, 0],
        a.clone(),
        vec![],
    );
    let defaults = SourceImageBudget::default();
    for budget in [
        SourceImageBudget {
            max_sources: 0,
            ..defaults
        },
        SourceImageBudget {
            max_terms: 0,
            ..defaults
        },
        SourceImageBudget {
            max_face_specializations: 0,
            ..defaults
        },
        SourceImageBudget {
            max_polynomial_degree: 0,
            ..defaults
        },
        SourceImageBudget {
            max_binding_bytes: 0,
            ..defaults
        },
        SourceImageBudget {
            max_source_bytes: 0,
            ..defaults
        },
    ] {
        assert!(certify_source_image_class(&c, &a, &b, budget, &RunContext::default()).is_err());
    }
    let run = RunContext::default();
    run.cancellation.cancel();
    assert!(matches!(
        certify_source_image_class(&c, &a, &b, defaults, &run),
        Err(Error::Cancelled)
    ));
}
#[test]
fn source_class_zero_domains_and_foreign_fixed_patterns_are_conservative() {
    let (_, _, a, b) = data();
    let c = context(Atom::one(), [0; 6], a.clone(), vec![]);
    let mut bad = *a.bounds();
    bad[2] = IndexBounds::fixed(1);
    let sources = copy_sources(&c, |_| {});
    let sources = Arc::try_unwrap(sources)
        .unwrap()
        .with_zero_domains(vec![IndexDomain::new(bad).unwrap()])
        .unwrap();
    assert!(
        certify_native(
            &Arc::new(sources),
            5,
            &a,
            &b,
            Default::default(),
            &RunContext::default()
        )
        .is_err()
    );
    // Native source construction rejects fixed raw patterns before the
    // checker can receive one. Preserve that explicit rejection contract.
    let original = c.sources();
    let mut row = original.native_sources().rows()[0].clone();
    row[0].integral = NativeIntegral::numeric([0; 6]).unwrap();
    let fixed = GuardedSourceSystem::new(
        original.measure_id(),
        *original.roles(),
        *original.native_sources().index_variables(),
        vec![GuardedSource::new("fixed", row, a.clone())],
    );
    assert!(
        matches!(fixed,Err(rustred::solver::SolverError::InvalidInput(ref m)) if m.contains("must have a symbolic integral index"))
    );
    let mut widened = *a.bounds();
    widened[2] = IndexBounds::unbounded();
    assert!(
        certify_source_image_class(
            &c,
            &IndexDomain::new(widened).unwrap(),
            &b,
            Default::default(),
            &RunContext::default()
        )
        .is_err()
    );
}
#[test]
fn source_class_cancellation_of_identical_images_does_not_hide_nonzero_bad_images() {
    let (_, _, a, b) = data();
    let c = context(Atom::one(), [0; 6], a.clone(), vec![]);
    let sources = copy_sources(&c, |rows| {
        let mut t = rows[0].row[0].clone();
        t.integral = NativeIntegral::symbolic([0, 0, 1, 0, 0, 0]).unwrap();
        rows[0].row = vec![
            t.clone(),
            rustred::solver::Term {
                integral: t.integral,
                coefficient: -t.coefficient,
            },
        ];
    });
    let p = certify_native(
        &sources,
        5,
        &a,
        &b,
        Default::default(),
        &RunContext::default(),
    )
    .unwrap();
    assert_eq!(p.counts.merged_terms, 0);
    // Even if another well-behaved row would reduce the bad child later, this
    // row's nonzero invalid intermediate is rejected before native discovery.
    let sources = copy_sources(&c, |rows| {
        let mut extra = rows[0].clone();
        extra.id = "bad-intermediate".into();
        extra.row[0].integral = NativeIntegral::symbolic([0, 0, 1, 0, 0, 0]).unwrap();
        rows.push(extra);
    });
    assert!(
        certify_native(
            &sources,
            5,
            &a,
            &b,
            Default::default(),
            &RunContext::default()
        )
        .is_err()
    );
}
#[test]
fn source_class_actual_singleton_and_both_partial_double_placements_cover_all_emitted_rows() {
    use super::super::DensityInput;
    use super::super::partial_origin::PartialOriginCapability;
    let input = serde_json::from_str::<DensityInput>(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/examples/finite_density/massless_three_loop_chain.json"
    )))
    .unwrap()
    .prepare()
    .unwrap();
    let options = WeightedSourceOptions {
        free_virtual_zero_sectors: true,
        policy: super::super::preparation::WeightedSourcePolicy::PolynomialClosure,
        ..Default::default()
    };
    for (cuts, shifted) in [
        (vec![3], vec![0]),
        (vec![0, 3], vec![1, 2]),
        (vec![0, 3], vec![2, 4]),
    ] {
        let family = input
            .occupied_cut(&cuts, 1024)
            .unwrap()
            .at_physical_masses();
        let origin = PartialOriginCapability::new_with_context(
            &input,
            &family,
            &shifted,
            options,
            Default::default(),
            &RunContext::default(),
        )
        .unwrap();
        origin
            .validate_binding(&input, &family, &shifted, options)
            .unwrap();
        let identity = origin.identity().unwrap();
        let formal = family
            .guarded_sources_with_options::<16>(
                symbol!("source_class_actual_eps"),
                4,
                symbol!("source_class_actual_eta"),
                &shifted,
                4096,
                vec![],
                GuardedMeasureIdentity {
                    measure: format!("actual partial source class; origin={identity}"),
                    support: "bound polynomial partial-origin class".into(),
                    orientation: "existing occupied routing".into(),
                    normalization: "existing C/H".into(),
                    branch: "fixed positive eta; regulated continuation".into(),
                    deformation: format!("native eta shifted {shifted:?}"),
                },
                WeightedSourceOptions {
                    free_virtual_zero_sectors: false,
                    ..options
                },
            )
            .unwrap();
        let original = formal.context.sources();
        let rows = original
            .sources()
            .iter()
            .zip(original.native_sources().rows())
            .map(|(i, r)| {
                GuardedSource::new(i.id.clone(), r.clone(), i.domain.clone())
                    .with_nonzero_conditions(i.nonzero_conditions.clone())
            })
            .collect();
        let sources = Arc::new(
            GuardedSourceSystem::new(
                format!("{};bound-origin={identity}", original.measure_id()),
                *original.roles(),
                *original.native_sources().index_variables(),
                rows,
            )
            .unwrap()
            .with_zero_domains(
                origin
                    .zero_domains_with_context::<16>(
                        &family,
                        &shifted,
                        options,
                        Default::default(),
                        &RunContext::default(),
                    )
                    .unwrap(),
            )
            .unwrap(),
        );
        assert_eq!(
            sources.native_sources().rows(),
            original.native_sources().rows()
        );
        let family_signature = family
            .factors()
            .iter()
            .map(Atom::to_canonical_string)
            .chain(family.coordinates().iter().map(Atom::to_canonical_string))
            .collect();
        let binding = SourceImageBinding {
            input_identity: input.identity().into(),
            family_signature,
            origin_identity: identity,
            shifted_slots: shifted.clone(),
            source_options: options,
            expected_roles: family.roles().to_vec(),
            physical_arity: family.factors().len(),
            physical_slots: family.physical_slots(),
            input_slots: family.input_slots(),
        };
        let proof = certify_native(
            &sources,
            family.factors().len(),
            formal.deformation.admitted_domain(),
            &binding,
            Default::default(),
            &RunContext::default(),
        )
        .unwrap();
        assert_eq!(proof.counts.sources, sources.sources().len());
        assert!(proof.counts.checked_faces > 0);
        if let Ok(path) = std::env::var("PARTIAL_SOURCE_CLASS_REPORT") {
            std::fs::create_dir_all(&path).unwrap();
            std::fs::write(std::path::Path::new(&path).join(format!("shift-{}.json",shifted.iter().map(usize::to_string).collect::<Vec<_>>().join("-"))),serde_json::to_vec_pretty(&json!({"certificate":proof.report(),"native_zero_boxes":sources.zero_domains().len(),"all_formal_rows_unchanged":true,"imported_rules":0,"physical_capability_identity":binding.origin_identity})).unwrap()).unwrap();
        }
    }
}

#[test]
fn source_class_retained_arithmetic_size_and_fixed_face_caps_are_enforced() {
    let (indices, _, a, b) = data();
    let c = context(Atom::one(), [0; 6], a.clone(), vec![]);
    let twice = copy_sources(&c, |rows| {
        let t = rows[0].row[0].clone();
        rows[0].row.push(t);
    });
    let cap = SourceImageBudget {
        max_integer_bits: 1,
        ..Default::default()
    };
    assert!(
        matches!(certify_native(&twice,5,&a,&b,cap,&RunContext::default()),Err(Error::Limit(ref m)) if m.contains("integer coefficient"))
    );
    let mut guard = *a.bounds();
    guard[0] = IndexBounds::fixed(65537);
    let huge = context(
        Atom::var(indices[0]),
        [0; 6],
        IndexDomain::new(guard).unwrap(),
        vec![],
    );
    assert!(matches!(cert(&huge),Err(Error::Limit(ref m)) if m.contains("fixed guard coordinate")));
    let mut inverse = b.clone();
    inverse.source_options.positive_compact_energy_powers = true;
    assert!(
        matches!(certify_source_image_class(&c,&a,&inverse,Default::default(),&RunContext::default()),Err(Error::Unsupported(ref m)) if m.contains("inverse completions"))
    );
}
