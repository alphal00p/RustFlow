use super::*;

#[test]
fn explicit_occupied_placement_binds_parent_physical_slots() {
    let input: crate::finite_density::DensityInput = serde_json::from_str(include_str!(
        "../../examples/finite_density/massless_three_loop_chain.json"
    ))
    .unwrap();
    let input = input.prepare().unwrap();
    let family = input
        .occupied_cut(&[0, 3], 4096)
        .unwrap()
        .at_physical_masses();
    let mut options = FlowOptions::default();
    for (mode, expected) in [
        (MassMode::All, vec![1, 2, 4]),
        (MassMode::Auto, vec![1, 2, 4]),
        (MassMode::Propagators(vec![4, 2]), vec![2, 4]),
    ] {
        options.mass_mode = mode;
        assert_eq!(selected_flow_slots(&family, &options).unwrap(), expected);
    }
    for slots in [vec![], vec![1, 1], vec![0, 2], vec![3], vec![5], vec![999]] {
        options.mass_mode = MassMode::Propagators(slots);
        assert!(matches!(
            selected_flow_slots(&family, &options),
            Err(Error::InvalidInput(_))
        ));
    }
    options.mass_mode = MassMode::Loop;
    assert!(matches!(
        selected_flow_slots(&family, &options),
        Err(Error::Unsupported(_))
    ));
    // A complete amplitude has several cuts: its options cannot reinterpret a
    // parent's explicit indices by dropping whichever slots happen to be cut.
    options.mass_mode = MassMode::Propagators(vec![1, 2]);
    assert!(matches!(
        validate_options(&options),
        Err(Error::Unsupported(_))
    ));
}

#[test]
fn selected_placement_does_not_erase_contour_or_recursion_options() {
    let input: crate::finite_density::DensityInput = serde_json::from_str(include_str!(
        "../../examples/finite_density/massless_three_loop_chain.json"
    ))
    .unwrap();
    let family = input
        .prepare()
        .unwrap()
        .occupied_cut(&[0, 3], 4096)
        .unwrap();
    let options = FlowOptions {
        mass_mode: MassMode::Propagators(vec![1, 2]),
        prescription: Prescription::MinusI0,
        ..Default::default()
    };
    assert!(matches!(
        selected_flow_slots(&family, &options),
        Err(Error::Unsupported(_))
    ));
}
