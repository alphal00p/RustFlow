use symbolica::atom::{Atom, AtomCore, NamespacedSymbol, SymbolBuilder};

fn main() {
    let original = Atom::num(1);
    let builder = SymbolBuilder::new(NamespacedSymbol::parse("export_review::generated"))
        .with_generator(move |_, builder| {
            eprintln!("entered generator");
            let mut bytes = Vec::new();
            original.export(&mut bytes).unwrap();
            eprintln!("export completed: {} bytes", bytes.len());
            builder
        });
    SymbolBuilder::build_group(vec![builder]).unwrap();
    eprintln!("group completed");
}
