pub(crate) fn imaginary_parameter() -> Symbol {
    symbol!("symbolica_amflow::imaginary_unit")
}

pub(crate) fn encode_complex(a: &Atom) -> Atom {
    a.replace_map(|view, _, out| {
        if let AtomView::Num(n) = view
            && let symbolica::coefficient::Coefficient::Complex(c) = n.get_coeff_view().to_owned()
            && !c.im.is_zero()
        {
            **out = Atom::num(c.re) + Atom::var(imaginary_parameter()) * Atom::num(c.im);
        }
    })
}

