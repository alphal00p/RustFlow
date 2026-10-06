//! Cumulative methodological and scientific credit for successful API preparation.
//! Native software dependencies retain their own usage tracking in the host.
use std::sync::atomic::{AtomicU8, Ordering};
use symbolica::api::python::Citation;

const TRANSPORT: u8 = 1;
const AUTOMATIC: u8 = 2;
const HIGGS_JET: u8 = 4;
static USED: AtomicU8 = AtomicU8::new(0);

pub(super) fn transport() {
    USED.fetch_or(TRANSPORT, Ordering::Relaxed);
}

#[cfg(feature = "automatic")]
pub(super) fn automatic() {
    USED.fetch_or(AUTOMATIC, Ordering::Relaxed);
}

pub(super) fn higgs_jet() {
    USED.fetch_or(HIGGS_JET, Ordering::Relaxed);
}

fn citation(id: &str, reference: &str, bibtex: &str, reasons: Vec<String>) -> Citation {
    Citation {
        id: id.into(),
        reference: reference.into(),
        bibtex: bibtex.into(),
        reasons,
        description:
            "Literature underlying the numerical integration methods and scientific inputs.".into(),
        relevance: None,
    }
}

/// Return cumulative citations for successfully prepared native integration APIs.
/// Module registration does not activate these references. Reading never resets
/// usage; the community host merges these native Citation objects by stable ID.
pub fn get_citations() -> Vec<Citation> {
    let used = USED.load(Ordering::Relaxed);
    if used == 0 {
        return Vec::new();
    }
    let mut continuation_reasons = Vec::new();
    if used & TRANSPORT != 0 {
        continuation_reasons.push("Prepared native differential-equation transport using local series and analytic continuation methods.".into());
    }
    if used & AUTOMATIC != 0 {
        continuation_reasons.push("Prepared native automatic integral evaluation combining auxiliary-mass flow and differential-equation continuation.".into());
    }
    if used & HIGGS_JET != 0 {
        continuation_reasons.push("Prepared the Higgs-jet workflow with native auxiliary-mass-flow starting values and physical series transport.".into());
    }
    let mut citations = vec![
        citation(
            "arXiv:2607.08477",
            "Rui-Jun Huang, Xiao Liu and Yan-Qing Ma. AMFlow 2.0: significant algorithmic and software improvements for Feynman integral evaluation (2026). arXiv:2607.08477.",
            r#"@article{Huang:2026AMFlow2,
  author = {Huang, Rui-Jun and Liu, Xiao and Ma, Yan-Qing},
  title = {{AMFlow 2.0}: significant algorithmic and software improvements for {Feynman} integral evaluation},
  eprint = {2607.08477},
  archivePrefix = {arXiv},
  primaryClass = {hep-ph},
  year = {2026},
  url = {https://arxiv.org/abs/2607.08477}
}"#,
            continuation_reasons.clone(),
        ),
        citation(
            "arXiv:2006.05510",
            "Martijn Hidding. DiffExp, a Mathematica package for computing Feynman integrals in terms of one-dimensional series expansions (2020). arXiv:2006.05510.",
            r#"@article{Hidding:2020DiffExp,
  author = {Hidding, Martijn},
  title = {{DiffExp}, a {Mathematica} package for computing {Feynman} integrals in terms of one-dimensional series expansions},
  eprint = {2006.05510},
  archivePrefix = {arXiv},
  primaryClass = {hep-ph},
  year = {2020},
  url = {https://arxiv.org/abs/2006.05510}
}"#,
            continuation_reasons,
        ),
    ];
    if used & (AUTOMATIC | HIGGS_JET) != 0 {
        let mut reasons = Vec::new();
        if used & AUTOMATIC != 0 {
            reasons.push("Auxiliary-mass-flow methodology for the prepared native automatic integral evaluator.".into());
        }
        if used & HIGGS_JET != 0 {
            reasons.push("Auxiliary-mass-flow methodology for native starting boundaries in the supplied-boundary Higgs-jet workflow.".into());
        }
        citations.push(citation(
            "arXiv:2201.11669",
            "Xiao Liu and Yan-Qing Ma. AMFlow: a Mathematica package for Feynman integrals computation via Auxiliary Mass Flow. Comput. Phys. Commun. 283 (2023) 108565. doi:10.1016/j.cpc.2022.108565.",
            r#"@article{Liu:2022AMFlow,
  author = {Liu, Xiao and Ma, Yan-Qing},
  title = {{AMFlow}: a {Mathematica} package for {Feynman} integrals computation via {Auxiliary Mass Flow}},
  eprint = {2201.11669},
  archivePrefix = {arXiv},
  primaryClass = {hep-ph},
  journal = {Computer Physics Communications},
  volume = {283},
  pages = {108565},
  year = {2023},
  doi = {10.1016/j.cpc.2022.108565},
  url = {https://arxiv.org/abs/2201.11669}
}"#,
            reasons,
        ));
    }
    if used & HIGGS_JET != 0 {
        citations.push(citation(
            "arXiv:2112.07578",
            "Matteo Becchetti, Francesco Moriello and Armin Schweitzer. Two-loop amplitude for mixed QCD-EW corrections to gg -> Hg. JHEP 04 (2022) 139. doi:10.1007/JHEP04(2022)139.",
            r#"@article{Becchetti:2021HiggsJet,
  author = {Becchetti, Matteo and Moriello, Francesco and Schweitzer, Armin},
  title = {Two-loop amplitude for mixed {QCD-EW} corrections to {$gg \to Hg$}},
  eprint = {2112.07578},
  archivePrefix = {arXiv},
  primaryClass = {hep-ph},
  journal = {JHEP},
  volume = {04},
  pages = {139},
  year = {2022},
  doi = {10.1007/JHEP04(2022)139},
  url = {https://arxiv.org/abs/2112.07578}
}"#,
            vec!["Scientific Higgs-jet differential systems, form-factor projectors and tensor structures used by the prepared native application.".into()],
        ));
    }
    citations
}
