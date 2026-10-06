//! Symbolica c340 owner-level export race, no RustFlow dependency.
//! The writer callback deterministically schedules unrelated registry insertion
//! after its serialized count is written but before the native iterator starts.
use std::io::{self, Write};
use symbolica::domains::finite_field::{FiniteFieldCore, Zp64};
use symbolica::prelude::*;

struct RegisteringWriter {
    bytes: Vec<u8>,
    mode: String,
    parameter: Symbol,
    triggered: bool,
    retained: Vec<Atom>,
}

impl Write for RegisteringWriter {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        let offset = self.bytes.len();
        self.bytes.extend_from_slice(bytes);
        // Fresh process, Atom::num(1), no coefficient registries: magic/version/
        // flags/symbol-count occupy 15 bytes; field-count and list-count follow.
        let trigger = match self.mode.as_str() {
            "finite" => offset == 15,
            "polynomial" => offset == 23,
            "control" => false,
            _ => unreachable!(),
        };
        if trigger {
            assert!(!self.triggered);
            assert_eq!(bytes, [0_u8; 8], "expected an initially empty registry");
            self.triggered = true;
            if self.mode == "finite" {
                let field = Zp64::new(7);
                let element = field.to_element(3);
                self.retained
                    .push(Atom::num(Coefficient::from_finite_field(field, element)));
            } else {
                self.retained
                    .push((Atom::var(self.parameter) + 1).set_coefficient_ring(self.parameter));
            }
            eprintln!(
                "registered {} state after count at offset {}",
                self.mode, offset
            );
        }
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

fn main() {
    let mode = std::env::args()
        .nth(1)
        .expect("mode: control | finite | polynomial");
    assert!(["control", "finite", "polynomial"].contains(&mode.as_str()));
    let parameter = symbol!("export_registry_race::p");
    let original = Atom::num(1);
    let mut writer = RegisteringWriter {
        bytes: Vec::new(),
        mode: mode.clone(),
        parameter,
        triggered: false,
        retained: Vec::new(),
    };
    original
        .export(&mut writer)
        .expect("native export itself succeeds");
    assert_eq!(writer.triggered, mode != "control");
    if let Some(path) = std::env::args().nth(2) {
        std::fs::write(path, &writer.bytes).unwrap();
    }
    println!(
        "mode={mode} bytes={} export=Ok callback={}",
        writer.bytes.len(),
        writer.triggered
    );
    let imported = std::panic::catch_unwind(|| {
        let mut source = writer.bytes.as_slice();
        Atom::import(&mut source, None).map(|atom| (atom, source.len()))
    });
    match imported {
        Ok(Ok((atom, trailing))) => {
            println!(
                "import=Ok atom={atom} trailing={trailing} equal={}",
                atom == original
            );
            assert_eq!(atom, original);
            assert_eq!(trailing, 0);
            // Patched native exporter must succeed in every mode.
        }
        Ok(Err(error)) => {
            println!("import=Err kind={:?} message={error}", error.kind());
            panic!("patched export must be importable in {mode} mode");
        }
        Err(_) => {
            println!("import=panic");
            panic!("patched export must be importable in {mode} mode");
        }
    }
    // The same new registry entries are harmless once registration has ended.
    let mut stable = Vec::new();
    original.export(&mut stable).unwrap();
    let mut source = stable.as_slice();
    let stable_atom = Atom::import(&mut source, None).unwrap();
    assert_eq!(stable_atom, original);
    assert!(source.is_empty());
    println!(
        "post_registration_control=Ok bytes={} atom={stable_atom} trailing=0",
        stable.len()
    );
}
